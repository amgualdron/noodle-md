"""Checks that the physics is what the docs say it is. Run by setup.sh; takes ~1 minute.

    env/bin/python src/check_physics.py

Prints PASS/FAIL per check and exits with an error if anything failed.
"""

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import gsd.hoomd
import hoomd
import numpy as np
import yaml

from forcefield import AA_PARAMS, AMINO_ACIDS, debye_huckel, thermal_energy, wall_contact_distance, wall_forces
from launch_run import resolve_physics
from runs import ROOT

ALL_SETS = yaml.safe_load((ROOT / "physics.yaml").read_text())
failures = []


def check(description, ok):
    print(f"  {'PASS' if ok else 'FAIL'}  {description}")
    if not ok:
        failures.append(description)


def wall_energy_and_force(physics, residue, height):
    """Energy and z-force on one residue at `height` above the surface (at z = -15)."""
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    snap = hoomd.Snapshot()
    snap.configuration.box = [20, 20, 40, 0, 0, 0]
    snap.particles.N = 1
    snap.particles.types = AMINO_ACIDS
    snap.particles.typeid[0] = AMINO_ACIDS.index(residue)
    snap.particles.position[0] = [0, 0, -15 + height]
    sim.create_state_from_snapshot(snap)
    forces = wall_forces(physics, z_wall=-15, z_ceiling=15)
    sim.operations.integrator = hoomd.md.Integrator(dt=0.001, forces=forces)
    sim.run(0)
    return sum(f.energy for f in forces), sum(f.forces[0][2] for f in forces)


def main():
    grafted = resolve_physics(ALL_SETS, "grafted")

    print("Electrostatics")
    bjerrum, kappa = debye_huckel(298, 0.150)
    check(f"Bjerrum length at 298 K is {bjerrum:.3f} nm (water: ~0.71)", 0.69 < bjerrum < 0.73)
    check(f"Debye length at 150 mM is {1 / kappa:.3f} nm (expected ~0.79)", 0.76 < 1 / kappa < 0.81)

    print("Neutral, purely repulsive surface (the default)")
    for potential in ("lj93", "lj126"):
        physics = dict(grafted, wall_potential=potential)
        for residue in ("G", "W"):
            contact = wall_contact_distance(AA_PARAMS[residue]["sigma"], physics)
            heights = np.linspace(0.8 * contact, contact, 20)
            energies, forces = zip(*(wall_energy_and_force(physics, residue, h) for h in heights))
            u_far, f_far = wall_energy_and_force(physics, residue, contact * 1.01)
            check(f"{potential} {residue}: zero energy and force beyond contact ({contact:.3f} nm)",
                  u_far == 0 and f_far == 0)
            check(f"{potential} {residue}: energy only rises toward the surface, force pushes away",
                  all(np.diff(energies) <= 0) and all(f >= 0 for f in forces))
            check(f"{potential} {residue}: energy goes smoothly to zero at contact",
                  energies[-1] < 1e-6 and energies[-2] < 1e-3)
        _, f_below = wall_energy_and_force(physics, "G", -0.05)
        check(f"{potential}: a bead that slipped past the surface is pushed back up", f_below > 0)

    print("Attractive surface (wall_attraction > 0)")
    for scale in ("uniform", "hydropathy"):
        physics = dict(grafted, wall_attraction=0.5, wall_attraction_scale=scale)
        lam = 0.5 if scale == "uniform" else 0.5 * AA_PARAMS["F"]["hydropathy"]
        contact = wall_contact_distance(AA_PARAMS["F"]["sigma"], physics)
        u_min, f_min = wall_energy_and_force(physics, "F", contact)
        expected = -physics["wall_epsilon"] * lam
        check(f"{scale}: well depth at contact is {u_min:.4f} kcal/mol (expected {expected:.4f})",
              abs(u_min - expected) < 0.02 * abs(expected))
        check(f"{scale}: no force at the bottom of the well", abs(f_min) < 1e-3)
    physics = dict(grafted, wall_attraction=0.5, wall_attraction_scale="hydropathy")
    contact = wall_contact_distance(AA_PARAMS["R"]["sigma"], physics)
    u_r, _ = wall_energy_and_force(physics, "R", contact + 0.3)
    check("hydropathy: arginine (hydropathy lambda = 0) feels no attraction", u_r == 0)

    print("Charged surface (Debye-Huckel)")
    physics = dict(grafted, wall_charge_density=-0.1)
    bjerrum, kappa = debye_huckel(physics["temperature_K"], physics["ionic_concentration_M"])
    kT = thermal_energy(physics["temperature_K"])
    for h in (1.0, 2.0):
        u, f = wall_energy_and_force(physics, "K", h)
        q_lysine = AA_PARAMS["K"]["charge"]
        analytic = 4 * np.pi * bjerrum * kT * q_lysine * (-0.1) / kappa * np.exp(-kappa * h)
        u_up, _ = wall_energy_and_force(physics, "K", h + 1e-4)
        u_down, _ = wall_energy_and_force(physics, "K", h - 1e-4)
        check(f"lysine at {h} nm: U = {u:.5f}, analytic {analytic:.5f} kcal/mol", abs(u - analytic) < 1e-9)
        check(f"lysine at {h} nm: force = -dU/dh", abs(f + (u_up - u_down) / 2e-4) < 1e-5)
    check("negative surface attracts lysine and repels glutamate",
          wall_energy_and_force(physics, "K", 1.0)[1] < 0 < wall_energy_and_force(physics, "E", 1.0)[1])

    print("Short simulations (bulk; grafted with attraction and charge; the grafted_charged set)")
    with tempfile.TemporaryDirectory() as tmp:
        run_dir = Path(tmp)
        (run_dir / "logs").mkdir()
        (run_dir / "trajectories").mkdir()
        sets = {"bulk": resolve_physics(ALL_SETS, "test_bulk"),
                "grafted": dict(resolve_physics(ALL_SETS, "test_grafted"),
                                wall_attraction=0.3, wall_charge_density=-0.05)}
        sets["grafted_C_terminus"] = dict(sets["grafted"], graft_terminus="C")
        sets["grafted_charged"] = resolve_physics(ALL_SETS, "test_grafted_charged")
        (run_dir / "physics.yaml").write_text(yaml.safe_dump(sets))
        sequence = "MSEQKKPLGEGFDNRYSTWAVHIC"
        names = list(sets)
        with open(run_dir / "tasks.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["task_id", "sequence_name", "physics", "replicate", "seed", "length", "sequence"])
            for i, name in enumerate(names, 1):
                writer.writerow([i, "test", name, 1, 42, len(sequence), sequence])
        for i, name in enumerate(names, 1):
            result = subprocess.run([sys.executable, str(ROOT / "src" / "simulate.py"), str(run_dir), str(i)],
                                    capture_output=True, text=True)
            log = json.loads((run_dir / "logs" / f"task_{i:04d}.json").read_text())
            check(f"{name}: task completes", result.returncode == 0 and log["status"] == "completed")
            if result.returncode != 0:
                print(result.stderr[-2000:])
                continue
            with gsd.hoomd.open(run_dir / log["trajectory"]) as traj:
                positions = np.array([frame.particles.position for frame in traj])
            check(f"{name}: {len(positions)} frames, all positions finite", np.isfinite(positions).all())
            if name.startswith("grafted"):
                geometry = log["geometry"]
                anchor = positions[:, geometry["anchor"], :]
                check(f"{name}: anchor (residue {geometry['anchor'] + 1}) never moved",
                      np.all(anchor == anchor[0]))
                check(f"{name}: every bead stayed between surface and ceiling",
                      positions[:, :, 2].min() > geometry["z_wall"] and positions[:, :, 2].max() < geometry["z_ceiling"])

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED.")
        sys.exit(1)
    print("All physics checks passed.")


if __name__ == "__main__":
    main()
