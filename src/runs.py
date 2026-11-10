"""Small helpers shared by the commands (./noodle new, launch, status, analyze, clean) and simulate.py.

A task's log (logs/task_0001.json) has a "status": "running" while it runs (or until
it is stopped), "completed" when done, or "failed" if it ended with an error (the
error message is then in its "error" field).

Every file in a run folder has exactly one program that writes it, so programs
running at the same time (./noodle launch on the login node, tasks on the cluster,
./noodle status) can never overwrite each other's information:
    run.yaml                      you
    tasks.csv, physics.yaml       launch_run.py, once, on the first launch
    slurm/submissions.csv         launch_run.py, one line per task each time it is submitted
    logs/task_0001.json           simulate.py, while that task runs
    analysis/                     analyze_run.py (also called by ./noodle clean)
    trajectories.sha256           clean_run.py, which also deletes trajectories of analysed tasks
"""

import csv
import getpass
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent  # the noodle-md folder
RUNS = ROOT / "runs"
SEQUENCES = ROOT / "sequences"  # the shared sequence sets


# ── Choosing a run ────────────────────────────────────────────────────────────
def list_runs():
    """Print every run folder with its description."""
    runs = sorted(p for p in RUNS.glob("*") if (p / "run.yaml").exists()) if RUNS.exists() else []
    if not runs:
        print("There are no runs yet. Make one with ./noodle new <run name>")
        return
    print("Existing runs:")
    for run_dir in runs:
        settings = yaml.safe_load((run_dir / "run.yaml").read_text()) or {}
        tasks = len(read_tasks(run_dir)) if (run_dir / "tasks.csv").exists() else 0
        state = f"{tasks} tasks" if tasks else "not launched"
        print(f"  {run_dir.name:32s} {state:>14s}   {settings.get('description') or ''}")


def run_from_command_line(usage):
    """The run folder named on the command line, or print usage and the runs and stop."""
    if len(sys.argv) != 2:
        print(f"Usage: {usage}\n")
        list_runs()
        sys.exit(1)
    run_dir = RUNS / sys.argv[1]
    if not (run_dir / "run.yaml").exists():
        print(f"There is no run called '{sys.argv[1]}'.\n")
        list_runs()
        sys.exit(1)
    return run_dir


# ── Reading and writing a run's files ─────────────────────────────────────────
def read_tasks(run_dir):
    """The rows of tasks.csv, as dicts (all values are strings)."""
    with open(run_dir / "tasks.csv") as f:
        return list(csv.DictReader(f))


def log_path(run_dir, task_id):
    return run_dir / "logs" / f"task_{int(task_id):04d}.json"


def read_log(run_dir, task_id):
    """A task's log as a dict ({} if the task has never started)."""
    path = log_path(run_dir, task_id)
    return json.loads(path.read_text()) if path.exists() else {}


def update_log(run_dir, task_id, **fields):
    """Add or change fields in a task's log. Only simulate.py calls this.

    The log is written to a temporary file which then replaces the old one in a single
    step, so ./noodle status never reads a half-written log.
    """
    path = log_path(run_dir, task_id)
    path.parent.mkdir(exist_ok=True)
    log = read_log(run_dir, task_id)
    log.update(fields)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(log, indent=2))
    os.replace(tmp, path)


def read_submissions(run_dir):
    """{task id: SLURM job id of its latest submission}, e.g. {17: '1234_17'}."""
    path = run_dir / "slurm" / "submissions.csv"
    if not path.exists():
        return {}
    with open(path) as f:
        # a task submitted several times has several lines: the last one wins
        return {int(row["task_id"]): row["slurm_job"] for row in csv.DictReader(f)}


# ── Brush size, run length and time ───────────────────────────────────────────
def brush_size(physics):
    """(chains, box width in nm): how many chains a task simulates, and how wide its box is.

    Bulk, and grafted with grafting_density 0: one chain in a box box_size wide.
    A brush (grafting_density > 0): `chains` chains on a square lattice, grafting_density
    chains per nm^2, so the box is sqrt(chains / grafting_density) wide and repeats
    periodically in x and y, like an infinite brush. At high density that box could be
    narrower than twice the interaction range (plus a margin), so that a residue would
    feel two copies of the same neighbour through the periodic boundary; then the lattice
    grows to the next square number of chains until it is wide enough.
    """
    density = physics.get("grafting_density") or 0  # absent in runs launched before brushes
    if physics["geometry"] != "grafted" or density == 0:
        return 1, physics["box_size"]
    side = math.isqrt(physics["chains"])
    cutoff = max(physics["lj_r_cut"], physics["yukawa_r_cut"], physics["wall_r_cut"])
    while side / math.sqrt(density) < 2 * (cutoff + 0.5):
        side += 1
    return side**2, side / math.sqrt(density)


def run_lengths(physics, n_residues):
    """Equilibration steps, production steps, and steps between saved frames.

    Lengths scale like the chain's relaxation time, ~ N^2.2 (Rouse scaling for a chain
    with excluded volume): steps = multiplier * N^2.2 / dt. A number given for
    equilibration_steps / production_steps overrides that.
    """
    dt = physics["dt"]
    equilibration = physics.get("equilibration_steps")
    if equilibration is None:
        equilibration = physics["rouse_multiplier"] * n_residues**2.2 / dt
    production = physics.get("production_steps")
    if production is None:
        production = physics["production_multiplier"] * n_residues**2.2 / dt
    save_every = int(production) // int(physics["frames"])
    if save_every < 1:
        raise ValueError(f"production_steps ({int(production)}) is smaller than frames ({physics['frames']})")
    return int(equilibration), int(production), save_every


def duration(seconds):
    """3725 -> '1h 02m', 200000 -> '2d 07h'."""
    if seconds < 60:
        return "<1m"
    minutes = int(seconds // 60)
    days, hours, minutes = minutes // 1440, minutes // 60 % 24, minutes % 60
    if days:
        return f"{days}d {hours:02d}h"
    return f"{hours}h {minutes:02d}m" if hours else f"{minutes}m"


def compress(ids):
    """[1, 2, 3, 7] -> '1-3,7', the way sbatch --array wants it."""
    parts, start = [], ids[0]
    for previous, current in zip(ids, ids[1:] + [None]):
        if current != previous + 1:
            parts.append(f"{start}-{previous}" if previous != start else str(start))
            start = current
    return ",".join(parts)


def file_size(n_bytes):
    """1.23e9 -> '1.2 GB'."""
    for size, unit in ((1e12, "TB"), (1e9, "GB"), (1e6, "MB")):
        if n_bytes >= size:
            return f"{n_bytes / size:.1f} {unit}"
    return f"{n_bytes / 1e3:.0f} kB"


def about(seconds):
    """An estimated duration: 'under a minute', '~1h 02m', '~2d 07h'."""
    return "under a minute" if seconds < 60 else f"~{duration(seconds)}"


# ── The SLURM queue ───────────────────────────────────────────────────────────
def slurm_queue(run_name):
    """{job id: state} for this run's tasks in the SLURM queue, e.g.
    {'1234_17': 'RUNNING', '1234_18': 'PENDING'}; None on a computer without SLURM.

    A job array's elements are listed one by one ('<job>_<array index>'), and a job
    leaves the queue as soon as it ends, for whatever reason.
    """
    if shutil.which("squeue") is None:
        return None
    out = subprocess.run(["squeue", "--user", getpass.getuser(), "--noheader", "--array",
                          "--format=%i %T", f"--name={run_name}"], capture_output=True, text=True)
    return dict(line.split() for line in out.stdout.splitlines() if line.strip())
