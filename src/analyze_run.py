"""./noodle analyze <run name>  --  per-frame shape of every completed trajectory.

For every protein and physics set, pools all completed replicates into
    runs/<run name>/timeseries/<physics set>/<protein>_pooled_timeseries.csv
with one row per frame (all lengths in nm):
    rg            radius of gyration, mass-weighted
    rg_sq         rg squared
    ree_x/y/z     end-to-end vector, from the first residue to the last
    ree_sq        its length squared
    asphericity   0 for a sphere-like shape, up to 1 for a rod
Tasks that are not completed are skipped, so it can be run on a run in progress.
"""

import MDAnalysis as mda
import numpy as np
import pandas as pd

from runs import read_log, read_tasks, run_from_command_line


def trajectory_metrics(gsd_path):
    """One row per saved frame of one trajectory."""
    universe = mda.Universe(str(gsd_path))
    chain = universe.select_atoms("all")
    rows = []
    for ts in universe.trajectory:
        universe.atoms.unwrap(reference="com")  # undo periodic wrapping, chain kept whole
        positions = chain.positions
        rg = chain.radius_of_gyration()
        ree = positions[-1] - positions[0]
        rows.append(dict(
            frame=ts.frame,
            rg=rg,
            rg_sq=rg**2,
            ree_sq=np.sum(ree**2), ree_x=ree[0], ree_y=ree[1], ree_z=ree[2],
            asphericity=chain.asphericity(),
        ))
    return pd.DataFrame(rows)


def main():
    run_dir = run_from_command_line("./noodle analyze <run name>")
    tasks = pd.DataFrame(read_tasks(run_dir))

    for (protein, physics), group in tasks.groupby(["sequence_name", "physics"], sort=False):
        pooled = []
        for _, task in group.iterrows():
            log = read_log(run_dir, task["task_id"])
            if log.get("status") != "completed":
                continue
            metrics = trajectory_metrics(run_dir / log["trajectory"])
            metrics.insert(0, "seed", task["seed"])
            metrics.insert(0, "task_id", task["task_id"])
            metrics.insert(0, "physics", physics)
            metrics.insert(0, "sequence_name", protein)
            pooled.append(metrics)
        if not pooled:
            print(f"  {protein} ({physics}): no completed tasks yet")
            continue
        out = run_dir / "timeseries" / physics / f"{protein}_pooled_timeseries.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        pd.concat(pooled, ignore_index=True).to_csv(out, index=False)
        print(f"  {protein} ({physics}): {len(pooled)} task(s) -> {out.relative_to(run_dir)}")


if __name__ == "__main__":
    main()
