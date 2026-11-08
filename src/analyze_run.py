"""./noodle analyze <run name>  --  the shape of every finished chain, frame by frame and on average.

For each protein and physics set, the frames of all its completed replicates are pooled
into one small compressed table:
    analysis/timeseries/<physics set>/<protein>.csv.gz   one row per frame of every replicate
and from those the summary is rebuilt every time:
    analysis/summary.csv   one row per protein and physics set: averages and error bars

The timeseries keep everything a statistic needs, so the trajectories themselves can
go to long-term storage (./noodle clean) and the analysis can still be redone, at
~60 KB per replicate instead of megabytes. Columns (lengths in nm):
    replicate     which replicate (independent simulation) the frame comes from
    frame         index of the saved frame within that replicate
    rg            radius of gyration, mass-weighted
    ree_x/y/z     end-to-end vector, from the first residue to the last
    asphericity   0 for a sphere-like shape, up to 1 for a rod
    com_height    grafted only: height of the chain's centre of mass above the surface

Replicates already in a table are skipped, and a replicate that finishes later is added
to its table, so this is quick to rerun while a run is going.
"""

import gzip
import os
import sys

import gsd.fl
import gsd.hoomd
import numpy as np
import pandas as pd

from runs import read_log, read_tasks, run_from_command_line

BLOCKS = 10  # a lone replicate's frames are cut into this many blocks to estimate its error bar
OBSERVABLES = ["rg", "rg_sq", "ree_sq", "asphericity", "ree_z", "com_height"]


# ── One trajectory ────────────────────────────────────────────────────────────
def trajectory_metrics(gsd_path, z_wall=None):
    """One row per saved frame of one trajectory (a single chain).

    The file holds positions wrapped into the periodic box. The chain is made whole
    again by walking along it: each bond is shorter than half the box, so its true
    vector is its minimum-image vector, and residue i sits at residue 0 plus the sum of
    the first i bonds.
    """
    with gsd.hoomd.open(gsd_path, "r") as traj:
        first = traj[0]
        masses = first.particles.mass.astype(np.float64)
        box = np.asarray(first.configuration.box[:3], dtype=np.float64)
    with gsd.fl.open(gsd_path, "r") as f:  # positions only: much faster than whole frames
        x = np.stack([f.read_chunk(i, "particles/position") for i in range(f.nframes)]).astype(np.float64)

    bonds = np.diff(x, axis=1)
    bonds -= box * np.round(bonds / box)
    x = np.concatenate([x[:, :1], x[:, :1] + np.cumsum(bonds, axis=1)], axis=1)

    weights = masses / masses.sum()
    com = np.einsum("i,fij->fj", weights, x)
    rel = x - com[:, None, :]
    gyration = np.einsum("i,fij,fik->fjk", weights, rel, rel)  # gyration tensor per frame
    moments = np.linalg.eigvalsh(gyration)
    rg_sq = moments.sum(axis=1)
    asphericity = 1.5 * np.sum((moments - moments.mean(axis=1, keepdims=True)) ** 2, axis=1) / rg_sq**2
    ree = x[:, -1] - x[:, 0]

    table = pd.DataFrame(dict(frame=np.arange(len(x)), rg=np.sqrt(rg_sq),
                              ree_x=ree[:, 0], ree_y=ree[:, 1], ree_z=ree[:, 2], asphericity=asphericity))
    if z_wall is not None:
        table["com_height"] = com[:, 2] - z_wall  # z is not periodic for a grafted chain
    return table


def timeseries_path(run_dir, protein, physics):
    return run_dir / "analysis" / "timeseries" / physics / f"{protein}.csv.gz"


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


def summary_row(protein, physics, n_residues, table):
    """Averages over all pooled frames. The error bar is the scatter between the
    replicates' own averages (independent by construction) when there are several, and
    the block error of the single replicate otherwise."""
    table = table.assign(rg_sq=table["rg"] ** 2,  # theory uses <Rg^2>, not <Rg>^2
                         ree_sq=table["ree_x"] ** 2 + table["ree_y"] ** 2 + table["ree_z"] ** 2)
    replicates = table.groupby("replicate")
    row = dict(sequence_name=protein, physics=physics, n_residues=n_residues,
               replicates=replicates.ngroups, frames=len(table))
    for name in OBSERVABLES:
        if name in ("ree_z", "com_height") and "com_height" not in table:
            row[name] = row[f"{name}_err"] = np.nan  # no meaning for a free chain
            continue
        means = replicates[name].mean()
        row[name] = means.mean()
        row[f"{name}_err"] = (means.std(ddof=1) / np.sqrt(len(means)) if len(means) > 1
                              else block_error(table[name].to_numpy()))
    return row


# ── The command ───────────────────────────────────────────────────────────────
def analyze(run_dir):
    """Add every completed replicate not yet in its protein's table, then rebuild the summary.

    Returns the ids of every task whose frames are in a table."""
    tasks = pd.DataFrame(read_tasks(run_dir))
    logs = {int(i): read_log(run_dir, i) for i in tasks["task_id"]}
    completed = sum(log.get("status") == "completed" for log in logs.values())

    analysed, summary, new, missing, errors = set(), [], 0, [], []
    for (protein, physics), group in tasks.groupby(["sequence_name", "physics"], sort=False):
        path = timeseries_path(run_dir, protein, physics)
        table = pd.read_csv(path) if path.exists() else None
        have = set(table["replicate"]) if table is not None else set()

        added = []
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
                frames = trajectory_metrics(trajectory, log["geometry"].get("z_wall"))
                if len(frames) != log["frames"]:
                    raise ValueError(f"{len(frames)} frames in the file, the log says {log['frames']}")
            except Exception as error:
                errors.append((task_id, protein, f"{type(error).__name__}: {error}"))
                continue
            frames.insert(0, "replicate", replicate)
            added.append(frames)
            analysed.add(task_id)

        if added:
            table = pd.concat(([table] if table is not None else []) + added, ignore_index=True)
            table = table.sort_values(["replicate", "frame"], kind="stable")
            write_atomically(table, path)
            new += len(added)
            if sys.stdout.isatty():
                print(f"\r  {new} replicate(s) analysed", end="")
        if table is not None:
            summary.append(summary_row(protein, physics, len(group["sequence"].iloc[0]), table))
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
        write_atomically(pd.DataFrame(summary), run_dir / "analysis" / "summary.csv")
        print(f"  -> analysis/summary.csv ({len(summary)} protein x physics rows), analysis/timeseries/")
    return analysed


def main():
    analyze(run_from_command_line("./noodle analyze <run name>"))


if __name__ == "__main__":
    main()
