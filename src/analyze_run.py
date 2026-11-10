"""./noodle analyze <run name>  --  the shape of every finished chain, frame by frame and on average.

For each protein and physics set, the frames of all its completed replicates are pooled
into one small compressed table:
    analysis/timeseries/<physics set>/<protein>.csv.gz   one row per frame (and chain) of every replicate
    analysis/profiles/<physics set>/<protein>.csv.gz     brushes only: density vs height
and from those the summary is rebuilt every time:
    analysis/summary.csv   one row per protein and physics set: averages and error bars

The timeseries keep everything a statistic needs, so the trajectories themselves can
go to long-term storage (./noodle clean) and the analysis can still be redone, at
~60 KB per replicate instead of megabytes. Columns (lengths in nm):
    replicate     which replicate (independent simulation) the frame comes from
    frame         index of the saved frame within that replicate
    chain         brushes only: which chain of the brush, 0 to chains-1
    rg            radius of gyration, mass-weighted
    ree_x/y/z     end-to-end vector, from the first residue to the last
    asphericity   0 for a sphere-like shape, up to 1 for a rod
    com_height    grafted only: height of the chain's centre of mass above the surface

Replicates already in a table are skipped, and a replicate that finishes later is added
to its table, so this is quick to rerun while a run is going.
"""

import gzip
import math
import os
import sys

import gsd.fl
import gsd.hoomd
import numpy as np
import pandas as pd
import yaml

from runs import brush_size, read_log, read_tasks, run_from_command_line

BLOCKS = 10  # a lone replicate's frames are cut into this many blocks to estimate its error bar
OBSERVABLES = ["rg", "rg_sq", "ree_sq", "asphericity", "ree_z", "com_height"]
PROFILE_BIN = 0.25  # nm, height resolution of a brush's density profile
FEW_SAMPLES = 20    # fewer independent samples than this per chain: the run was too short
CHUNK_FRAMES = 100  # frames read at a time (a 10 x 10 brush of N=200: ~50 MB per 100 frames)


# ── One trajectory ────────────────────────────────────────────────────────────
def trajectory_metrics(gsd_path, chains=1, z_wall=None):
    """(one row per saved frame and chain, the density profile) of one trajectory.

    The profile is None except for a brush. The file holds positions wrapped into the
    periodic box. Each chain is made whole again by walking along it: each bond is
    shorter than half the box, so its true vector is its minimum-image vector, and
    residue i sits at residue 0 plus the sum of the first i bonds.

    The frames are read CHUNK_FRAMES at a time: a 10 x 10 brush of 200-residue chains
    has 20,000 beads per frame, and all its frames at once would take gigabytes.
    """
    with gsd.hoomd.open(gsd_path, "r") as traj:
        first = traj[0]
        n = first.particles.N // chains  # residues per chain; chain k is particles k*n ... k*n + n-1
        masses = first.particles.mass[:n].astype(np.float64)
        charges = first.particles.charge[:n].astype(np.float64)
        box = np.asarray(first.configuration.box[:3], dtype=np.float64)
    weights = masses / masses.sum()
    # profile bins up to the top of the box (empty bins at the top are trimmed at the end)
    edges = np.arange(0, box[2] + PROFILE_BIN, PROFILE_BIN)
    count, charge = np.zeros(len(edges) - 1), np.zeros(len(edges) - 1)

    parts = []
    with gsd.fl.open(gsd_path, "r") as f:  # positions only: much faster than whole frames
        frames = f.nframes
        for start in range(0, frames, CHUNK_FRAMES):
            stop = min(start + CHUNK_FRAMES, frames)
            x = np.stack([f.read_chunk(i, "particles/position") for i in range(start, stop)]).astype(np.float64)
            x = x.reshape(len(x), chains, n, 3)  # frame, chain, residue, xyz

            bonds = np.diff(x, axis=2)
            bonds -= box * np.round(bonds / box)
            x = np.concatenate([x[:, :, :1], x[:, :, :1] + np.cumsum(bonds, axis=2)], axis=2)

            com = np.einsum("i,fcij->fcj", weights, x)
            rel = x - com[:, :, None, :]
            gyration = np.einsum("i,fcij,fcik->fcjk", weights, rel, rel)  # gyration tensor per frame and chain
            moments = np.linalg.eigvalsh(gyration)
            rg_sq = moments.sum(axis=-1)
            asphericity = 1.5 * np.sum((moments - moments.mean(axis=-1, keepdims=True)) ** 2, axis=-1) / rg_sq**2
            ree = x[:, :, -1] - x[:, :, 0]

            part = pd.DataFrame(dict(frame=np.repeat(np.arange(start, stop), chains),
                                     chain=np.tile(np.arange(chains), stop - start),
                                     rg=np.sqrt(rg_sq).ravel(), ree_x=ree[..., 0].ravel(),
                                     ree_y=ree[..., 1].ravel(), ree_z=ree[..., 2].ravel(),
                                     asphericity=asphericity.ravel()))
            if z_wall is not None:
                part["com_height"] = (com[..., 2] - z_wall).ravel()  # z is not periodic for grafted chains
            parts.append(part)
            if chains > 1:
                heights = x[..., 2] - z_wall
                count += np.histogram(heights, edges)[0]
                charge += np.histogram(heights, edges, weights=np.broadcast_to(charges, heights.shape))[0]

    table = pd.concat(parts, ignore_index=True)
    if chains == 1:
        return table.drop(columns="chain"), None
    return table, density_profile(edges, count, charge, frames * box[0] * box[1] * PROFILE_BIN)


def density_profile(edges, count, charge, volume):
    """Residues and charge per nm^3 at each height above the surface, averaged over the
    frames: the histogram of every residue's height, divided by the volume of each slab
    times the number of frames. Bins above the highest residue are left out."""
    top = np.flatnonzero(count)[-1] + 1 if count.any() else 0
    return pd.DataFrame(dict(z=edges[:top] + PROFILE_BIN / 2, density=count[:top] / volume,
                             charge=charge[:top] / volume))


def output_path(run_dir, kind, protein, physics):
    """analysis/<kind>/<physics set>/<protein>.csv.gz, kind being timeseries or profiles."""
    return run_dir / "analysis" / kind / physics / f"{protein}.csv.gz"


def write_atomically(table, path):
    """Write to a temporary file, then put it in place in a single step: a file that
    exists is always complete (./noodle clean relies on this before deleting anything)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    if path.suffix == ".gz":
        with gzip.open(tmp, "wt", compresslevel=6) as f:
            table.to_csv(f, index=False, float_format="%.5g")
    else:
        table.to_csv(tmp, index=False, float_format="%.6g")
    os.replace(tmp, path)


def add_rows(path, new, order):
    """Add the new replicates' rows to the table at path (created if missing), replacing
    any rows those replicates already had there."""
    tables = new
    if path.exists():
        old = pd.read_csv(path)
        tables = [old[~old["replicate"].isin(pd.concat(new)["replicate"])]] + new
    table = pd.concat(tables, ignore_index=True).sort_values(order, kind="stable")
    write_atomically(table, path)
    return table


# ── Statistics ────────────────────────────────────────────────────────────────
def block_error(values):
    """Standard error of the mean of correlated frames: the frames are cut into BLOCKS
    consecutive blocks, long enough to be roughly independent, and the error comes from
    how much the block averages scatter."""
    n = len(values) // BLOCKS * BLOCKS
    if n == 0:
        return np.nan
    means = values[:n].reshape(BLOCKS, -1).mean(axis=1)
    return means.std(ddof=1) / np.sqrt(BLOCKS)


def correlation_time(values):
    """Integrated autocorrelation time of a series, in frames: how many frames apart two
    values must be to be roughly independent (twice this). The autocorrelation function
    is summed up to where it first drops below zero, where noise starts to dominate."""
    x = values - values.mean()
    if not x.any():
        return np.nan
    spectrum = np.fft.rfft(x, 2 * len(x))
    acf = np.fft.irfft(spectrum * np.conj(spectrum))[:len(x)]
    acf /= acf[0]
    first_negative = np.argmax(acf < 0) or len(acf)
    return 0.5 + acf[1:first_negative].sum()


def summary_row(protein, physics_name, physics, n_residues, table):
    """Averages over all pooled frames (and chains). The error bar is the scatter between
    the replicates' own averages (independent by construction) when there are several,
    and the block error of the single replicate otherwise (from its chain-averaged frames).

    independent_samples: how many independent Rg values each chain contributed per
    replicate (its frames divided by twice the correlation time of its Rg; the median
    over chains and replicates). It shows whether the run was long enough: below
    FEW_SAMPLES, averages and error bars are not to be trusted.
    """
    table = table.assign(rg_sq=table["rg"] ** 2,  # theory uses <Rg^2>, not <Rg>^2
                         ree_sq=table["ree_x"] ** 2 + table["ree_y"] ** 2 + table["ree_z"] ** 2)
    chain = ["chain"] if "chain" in table else []
    per_frame = table.groupby(["replicate", "frame"]).mean(numeric_only=True)  # average over the chains
    replicates = per_frame.groupby("replicate")
    chains, _ = brush_size(physics)
    row = dict(sequence_name=protein, physics=physics_name, n_residues=n_residues, chains=chains,
               grafting_density=physics.get("grafting_density") or 0, reduced_density=np.nan,
               replicates=replicates.ngroups, frames=len(per_frame))
    for name in OBSERVABLES:
        if name in ("ree_z", "com_height") and "com_height" not in table:
            row[name] = row[f"{name}_err"] = np.nan  # no meaning for a free chain
            continue
        means = replicates[name].mean()
        row[name] = means.mean()
        row[f"{name}_err"] = (means.std(ddof=1) / np.sqrt(len(means)) if len(means) > 1
                              else block_error(per_frame[name].to_numpy()))
    samples = [len(rg) / (2 * correlation_time(rg.to_numpy()))
               for _, rg in table.groupby(["replicate"] + chain)["rg"]]
    row["independent_samples"] = np.nanmedian(samples)
    return row


def add_reduced_density(summary, physics):
    """Grafting density relative to overlap, sigma * pi * <Rg^2>, with the <Rg^2> of the
    same protein as a single chain in this run (a bulk set if there is one, else a single
    grafted chain). Around 1 the chains start to touch: below it they are "mushrooms",
    well above it a brush. Left empty when the run has no single-chain set."""
    single = {}
    for geometry in ("grafted", "bulk"):  # bulk last, so it wins
        for row in summary:
            if row["chains"] == 1 and physics[row["physics"]]["geometry"] == geometry:
                single[row["sequence_name"]] = row["rg_sq"]
    for row in summary:
        if row["chains"] > 1 and row["sequence_name"] in single:
            row["reduced_density"] = row["grafting_density"] * math.pi * single[row["sequence_name"]]


# ── The command ───────────────────────────────────────────────────────────────
def analyze(run_dir):
    """Add every completed replicate not yet in its protein's table, then rebuild the summary.

    Returns the ids of every task whose frames are in a table."""
    tasks = pd.DataFrame(read_tasks(run_dir))
    physics = yaml.safe_load((run_dir / "physics.yaml").read_text())
    logs = {int(i): read_log(run_dir, i) for i in tasks["task_id"]}
    completed = sum(log.get("status") == "completed" for log in logs.values())

    analysed, summary, new, missing, errors = set(), [], 0, [], []
    for (protein, physics_name), group in tasks.groupby(["sequence_name", "physics"], sort=False):
        path = output_path(run_dir, "timeseries", protein, physics_name)
        table = pd.read_csv(path) if path.exists() else None
        have = set(table["replicate"]) if table is not None else set()

        added, profiles = [], []
        for _, task in group.iterrows():
            task_id, replicate = int(task["task_id"]), int(task["replicate"])
            log = logs[task_id]
            if replicate in have:
                analysed.add(task_id)
                continue
            if log.get("status") != "completed":
                continue
            trajectory = run_dir / log["trajectory"]
            if not trajectory.exists():
                missing.append(task_id)
                continue
            try:
                frames, profile = trajectory_metrics(trajectory, log.get("chains", 1), log["geometry"].get("z_wall"))
                if frames["frame"].nunique() != log["frames"]:
                    raise ValueError(f"{frames['frame'].nunique()} frames in the file, the log says {log['frames']}")
            except Exception as error:
                errors.append((task_id, protein, f"{type(error).__name__}: {error}"))
                continue
            frames.insert(0, "replicate", replicate)
            added.append(frames)
            if profile is not None:
                profile.insert(0, "replicate", replicate)
                profiles.append(profile)
            analysed.add(task_id)

        if added:
            # the profile first: a replicate counts as analysed once it is in the timeseries
            if profiles:
                add_rows(output_path(run_dir, "profiles", protein, physics_name), profiles, ["replicate", "z"])
            order = ["replicate", "frame"] + (["chain"] if "chain" in added[0] else [])
            table = add_rows(path, added, order)
            new += len(added)
            if sys.stdout.isatty():
                print(f"\r  {new} replicate(s) analysed", end="")
        if table is not None:
            summary.append(summary_row(protein, physics_name, physics[physics_name],
                                       len(group["sequence"].iloc[0]), table))
    if new and sys.stdout.isatty():
        print()

    print(f"{completed} of {len(tasks)} tasks completed; {new} newly analysed, {len(analysed)} analysed in total")
    if errors:
        print(f"  {len(errors)} could not be analysed (their trajectory is kept):")
        for task_id, protein, message in errors[:10]:
            print(f"    task {task_id} ({protein}): {message}")
    if missing:
        print(f"  no trajectory here for {len(missing)} completed task(s) (moved to storage before being "
              f"analysed?): copy them back into runs/{run_dir.name}/trajectories/ to analyse them")
    if summary:
        add_reduced_density(summary, physics)
        short = [r for r in summary if r["independent_samples"] < FEW_SAMPLES]
        if short:
            print(f"  [!] {len(short)} protein x physics row(s) have fewer than {FEW_SAMPLES} independent samples "
                  f"per chain (column independent_samples): their runs were too short to trust, e.g. "
                  + ", ".join(f"{r['sequence_name']} ({r['physics']}, {r['independent_samples']:.0f})" for r in short[:3]))
        write_atomically(pd.DataFrame(summary), run_dir / "analysis" / "summary.csv")
        print(f"  -> analysis/summary.csv ({len(summary)} protein x physics rows), analysis/timeseries/"
              + (", analysis/profiles/" if any(r["chains"] > 1 for r in summary) else ""))
    return analysed


def main():
    analyze(run_from_command_line("./noodle analyze <run name>"))


if __name__ == "__main__":
    main()
