"""Run one task of a run: build the chain, equilibrate it, then record a trajectory.

    env/bin/python src/simulate.py runs/<run name> <task id>

You normally never call this yourself: ./launch-run does it, once per task.

Everything a task needs comes from its run folder:
    tasks.csv      which sequence, physics set and random seed this task uses
    physics.yaml   the physics settings, frozen when the run was launched
and everything it produces goes back there:
    logs/task_0001.json               status and progress (read by ./status)
    trajectories/task_0001_<name>.gsd production frames
    checkpoints/task_0001.gsd         full state saved every hour, deleted when done

If a task is stopped (e.g. it ran out of cluster time), running it again resumes
from its last checkpoint instead of starting over.

HOOMD in brief: a `hoomd.Simulation` holds the `state` (box, particles, bonds) and the
`operations`: an integrator that moves the particles under the forces from
forcefield.py, and writers that save frames to a file. `sim.run(n)` advances n steps.
"""

import os
import socket
import sys
import time
from datetime import datetime
from pathlib import Path

import gsd.fl
import gsd.hoomd
import hoomd
import hoomd.md
import numpy as np
import yaml

from forcefield import (AA_PARAMS, AMINO_ACIDS, residue_pair_forces, thermal_energy,
                        wall_contact_distance, wall_forces)
from runs import read_log, read_tasks, run_lengths, update_log

REPORT_EVERY_HOURS = 1     # how often to print progress, update the log and save a checkpoint
FIRST_REPORT_SECONDS = 60  # one early report, so you can see the speed right away
CHUNK_STEPS = 10_000       # steps between checks of the clock (cheap; not a physics setting)


# ── Geometry and starting configuration ───────────────────────────────────────
def system_geometry(sequence, physics):
    """The simulation box and, for a grafted chain, where the surface, ceiling and anchor are.

    bulk:    one chain in a large periodic cube.
    grafted: the surface sits near the bottom of the box. A purely repulsive ceiling
             sits above it, higher than the fully stretched chain can reach (contour
             length plus wall_headroom), so beads can never wrap through the periodic
             boundary in z. Above the ceiling there is one interaction cutoff of empty
             space, so beads near the ceiling never interact with beads near the surface
             through that boundary either.
    """
    if physics["geometry"] == "bulk":
        return dict(box=[physics["box_size"]] * 3)
    if physics["geometry"] != "grafted":
        raise ValueError(f"unknown geometry '{physics['geometry']}' (use bulk or grafted)")

    if physics["graft_terminus"] not in ("N", "C"):
        raise ValueError("graft_terminus must be N or C")
    contour = (len(sequence) - 1) * physics["harmonic_r0"]
    cutoff = max(physics["lj_r_cut"], physics["yukawa_r_cut"], physics["wall_r_cut"])
    height = contour + physics["wall_headroom"]  # surface to ceiling
    box_z = height + cutoff + 1.0
    z_wall = -box_z / 2 + 0.5
    return dict(box=[physics["box_size"], physics["box_size"], box_z],
                z_wall=z_wall, z_ceiling=z_wall + height,
                # the grafted residue: the first one (N-terminus) or the last one (C-terminus)
                anchor=0 if physics["graft_terminus"] == "N" else len(sequence) - 1)


def grafted_start_positions(sequence, physics, geometry, rng):
    """A random, non-overlapping chain growing upward from the anchor.

    The anchor sits at its contact distance from the surface, where the wall exerts no
    force on it. Each following bead is one bond length away in a random upward
    direction, retried until it is at least 0.5 nm from every earlier bead.
    """
    r0 = physics["harmonic_r0"]
    anchor_sigma = AA_PARAMS[sequence[geometry["anchor"]]]["sigma"]
    z_anchor = geometry["z_wall"] + wall_contact_distance(anchor_sigma, physics)

    grown = [np.array([0.0, 0.0, z_anchor])]
    while len(grown) < len(sequence):
        for _ in range(10_000):
            direction = rng.normal(size=3)
            direction /= np.linalg.norm(direction)
            direction[2] = abs(direction[2])
            candidate = grown[-1] + r0 * direction
            earlier = np.array(grown[:-1]).reshape(-1, 3)
            if np.all(np.linalg.norm(earlier - candidate, axis=1) >= 0.5):
                break
        else:
            raise RuntimeError("could not grow a non-overlapping starting chain")
        grown.append(candidate)

    positions = np.array(grown)
    # beads were grown from the anchor outward; for C-terminal grafting that is backwards
    return positions if geometry["anchor"] == 0 else positions[::-1]


def starting_snapshot(sequence, physics, geometry, seed):
    """The starting configuration: one bead per residue, bonded in a chain.

    A HOOMD snapshot is a plain description of the system: the box, and per particle
    its type, mass, charge and position, plus the list of bonds.
    """
    n = len(sequence)
    if physics["geometry"] == "bulk":
        # started as a straight line along x; equilibration makes it a random coil
        r0 = physics["harmonic_r0"]
        positions = [[i * r0 - n * r0 / 2, 0, 0] for i in range(n)]
    else:
        positions = grafted_start_positions(sequence, physics, geometry, np.random.default_rng(seed))

    snapshot = hoomd.Snapshot()
    snapshot.configuration.box = geometry["box"] + [0, 0, 0]  # Lx, Ly, Lz, no tilt
    snapshot.particles.N = n
    snapshot.particles.types = AMINO_ACIDS
    snapshot.particles.typeid[:] = [AMINO_ACIDS.index(aa) for aa in sequence]
    snapshot.particles.mass[:] = [AA_PARAMS[aa]["mass"] for aa in sequence]
    snapshot.particles.charge[:] = [AA_PARAMS[aa]["charge"] for aa in sequence]
    snapshot.particles.diameter[:] = [AA_PARAMS[aa]["sigma"] for aa in sequence]  # for viewers only
    snapshot.particles.position[:] = positions
    snapshot.bonds.N = n - 1
    snapshot.bonds.types = ["AA-bond"]
    snapshot.bonds.typeid[:] = 0
    snapshot.bonds.group[:] = [[i, i + 1] for i in range(n - 1)]
    return snapshot


# ── Simulation and forces ─────────────────────────────────────────────────────
def create_simulation(seed):
    """An empty HOOMD simulation, on the GPU if this HOOMD was built for one.

    The seed drives the thermostat's random kicks."""
    device = hoomd.device.GPU() if hoomd.version.gpu_enabled else hoomd.device.CPU()
    return hoomd.Simulation(device=device, seed=seed)


def add_forces(sim, physics, geometry):
    """Give a simulation (whose particles already exist) its forces and its integrator."""
    # Tree neighbour list: efficient for one chain in a big, mostly empty box. The
    # buffer makes it rebuild only when some bead has moved more than 0.1 nm.
    nlist = hoomd.md.nlist.Tree(buffer=0.2, exclusions=["bond"])
    forces = residue_pair_forces(physics, nlist)
    moving = hoomd.filter.All()
    if physics["geometry"] == "grafted":
        forces += wall_forces(physics, geometry["z_wall"], geometry["z_ceiling"])
        # The anchor is left out of the integrator: it feels forces but never moves.
        moving = hoomd.filter.SetDifference(hoomd.filter.All(), hoomd.filter.Tags([geometry["anchor"]]))
    # Langevin thermostat: friction plus random kicks that keep the temperature.
    # Friction gamma = 1 for every residue (HOOMD's default), i.e. a velocity
    # relaxation time m/gamma of ~100 time units (~50 ps) for a typical residue.
    langevin = hoomd.md.methods.Langevin(filter=moving, kT=thermal_energy(physics["temperature_K"]))
    sim.operations.integrator = hoomd.md.Integrator(dt=physics["dt"], methods=[langevin], forces=forces)


# ── Checkpoints ───────────────────────────────────────────────────────────────
def save_checkpoint(sim, path):
    """Save the full state (positions and velocities, in double precision).

    Written to a temporary file first, so a task killed while saving still has its
    previous checkpoint."""
    tmp = path.with_suffix(".tmp")
    hoomd.write.GSD.write(state=sim.state, filename=str(tmp), mode="wb", precision="double")
    os.replace(tmp, path)


def drop_frames_after(trajectory, step):
    """When resuming: keep only frames up to the checkpoint; the rest is recomputed.

    The frames are copied chunk by chunk with the low-level gsd.fl module, so the new
    file is in exactly the format HOOMD wrote and HOOMD can keep appending to it.
    """
    tmp = trajectory.with_suffix(".tmp")
    with gsd.fl.open(trajectory, "r") as old, gsd.fl.open(
            tmp, "w", application=old.application, schema=old.schema,
            schema_version=old.schema_version) as new:
        names = old.find_matching_chunk_names("")
        for frame in range(old.nframes):
            if old.read_chunk(frame, "configuration/step")[0] > step:
                break
            for chunk in names:
                if old.chunk_exists(frame, chunk):
                    new.write_chunk(chunk, old.read_chunk(frame, chunk))
            new.end_frame()
    os.replace(tmp, trajectory)


def now():
    return datetime.now().isoformat(timespec="seconds")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    run_dir, task_id = Path(sys.argv[1]), int(sys.argv[2])
    task = {int(t["task_id"]): t for t in read_tasks(run_dir)}[task_id]
    physics = yaml.safe_load((run_dir / "physics.yaml").read_text())[task["physics"]]
    name, sequence, seed = task["sequence_name"], task["sequence"], int(task["seed"])

    stem = f"task_{task_id:04d}"
    trajectory = run_dir / "trajectories" / f"{stem}_{name}.gsd"
    checkpoint = run_dir / "checkpoints" / f"{stem}.gsd"
    trajectory.parent.mkdir(exist_ok=True)
    checkpoint.parent.mkdir(exist_ok=True)

    log = read_log(run_dir, task_id)
    if log.get("status") == "completed":
        print(f"Task {task_id} ({name}) is already completed.")
        return

    equilibration, production, save_every = run_lengths(physics, len(sequence))
    total = equilibration + production
    geometry = system_geometry(sequence, physics)

    # ── System: resume from the last checkpoint, or start fresh ──
    sim = create_simulation(seed)
    if checkpoint.exists():
        sim.create_state_from_gsd(filename=str(checkpoint))
        if trajectory.exists():
            drop_frames_after(trajectory, sim.timestep)
        print(f"Resuming from checkpoint at step {sim.timestep:,} of {total:,}")
    else:
        sim.create_state_from_snapshot(starting_snapshot(sequence, physics, geometry, seed))
        trajectory.unlink(missing_ok=True)  # leftovers from an attempt that never checkpointed

    add_forces(sim, physics, geometry)

    # Frames every save_every steps, counted from the start of production
    writer = hoomd.write.GSD(filename=str(trajectory), mode="ab",
                             trigger=hoomd.trigger.Periodic(save_every, phase=equilibration % save_every))

    update_log(
        run_dir, task_id, status="running", sequence_name=name, n_residues=len(sequence),
        physics=task["physics"], geometry=geometry, seed=seed,
        equilibration_steps=equilibration, production_steps=production, total_steps=total,
        step=sim.timestep, trajectory=str(trajectory.relative_to(run_dir)),
        started_at=log.get("started_at", now()), attempts=log.get("attempts", 0) + 1,
        host=socket.gethostname(), hoomd_version=hoomd.version.version, last_update=now(),
    )
    print(f"Task {task_id}: {name}, {len(sequence)} residues, physics '{task['physics']}' "
          f"({physics['geometry']}), seed {seed}")
    print(f"  {equilibration:,} equilibration + {production:,} production steps, "
          f"a frame every {save_every:,} steps -> {trajectory}", flush=True)

    # ── Run, saving progress every REPORT_EVERY_HOURS ──
    # time.monotonic() is a stopwatch: unlike the wall clock it never jumps
    previous_run_time = log.get("run_time_s", 0)  # from earlier attempts of this task
    session_start = time.monotonic()
    last_report_time, last_report_step = session_start, sim.timestep
    next_report_time = session_start + FIRST_REPORT_SECONDS

    while sim.timestep < total:
        if sim.timestep >= equilibration and writer not in sim.operations.writers:
            sim.operations.writers.append(writer)  # production starts: begin saving frames
        end_of_phase = equilibration if sim.timestep < equilibration else total
        sim.run(min(CHUNK_STEPS, end_of_phase - sim.timestep))

        if time.monotonic() >= next_report_time:
            writer.flush()
            save_checkpoint(sim, checkpoint)
            clock = time.monotonic()
            speed = (sim.timestep - last_report_step) / (clock - last_report_time)
            phase = "equilibration" if sim.timestep < equilibration else "production"
            update_log(run_dir, task_id, step=sim.timestep, phase=phase,
                       steps_per_second=round(speed, 1),
                       run_time_s=round(previous_run_time + clock - session_start),
                       last_update=now())
            print(f"[{now()}] {phase:13s} {100 * sim.timestep / total:5.1f}%  "
                  f"{speed:,.0f} steps/s  ~{(total - sim.timestep) / speed / 3600:.1f} h left", flush=True)
            last_report_time, last_report_step = clock, sim.timestep
            next_report_time = clock + REPORT_EVERY_HOURS * 3600
    writer.flush()

    # ── Done ──
    with gsd.hoomd.open(trajectory, "r") as f:
        frames = len(f)
    run_time = round(previous_run_time + time.monotonic() - session_start)
    update_log(run_dir, task_id, status="completed", step=sim.timestep, phase="done", frames=frames,
               run_time_s=run_time, completed_at=now(), last_update=now())
    checkpoint.unlink(missing_ok=True)
    print(f"Done: {frames} frames in {run_time / 3600:.2f} h of run time")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        # Record the error for ./status and ./launch-run, then stop with the full
        # message (it ends up in the task's slurm/ output file, or on screen).
        update_log(Path(sys.argv[1]), int(sys.argv[2]), status="failed",
                   error=f"{type(error).__name__}: {error}", last_update=now())
        raise
