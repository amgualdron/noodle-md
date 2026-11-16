"""Small helpers shared by the commands (./noodle new, launch, status, analyze, clean, ...) and simulate.py.

Tasks and jobs. A task is one simulation (one protein, one physics set, one replicate).
A job is what SLURM runs: one element of a job array, one CPU core, running a few tasks
one after another (or a single long one). ./noodle launch groups the tasks into jobs.

Every file in a run folder has one program that writes it at any given time, so programs
running at the same time (./noodle commands on the login node, jobs on the cluster) never
overwrite each other's information:
    run.yaml                     you
    tasks.csv, physics.yaml      launch_run.py, once, at the first launch
    manifest.csv                 launch_run.py creates it; analyze_run.py adds the averages
    jobs/submitted.csv           launch_run.py, one line per job it submits
    jobs/job_0007.csv            launch_run.py writes it just before submitting job 7; from
                                 then on only job 7 itself (simulate.py) updates it
    jobs/job_0007_profiles.csv   job 7, the profiles of each task it analysed
    jobs/job_0007.out            everything job 7 printed
    log.csv                      status.py: every task's current state, merged from jobs/
    analysis/                    analyze_run.py
    trajectories.sha256          clean_run.py

A task's current state is its row in the newest job file that lists it: a task that is
resumed or retried gets a row in a new job, and its rows in older jobs are history.

runs/runs.csv, the list of every run with its number, is written only by ./noodle commands
(new, launch, status, analyze), never by jobs.
"""

import csv
import getpass
import math
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent  # the noodle-md folder
RUNS = ROOT / "runs"
SEQUENCES = ROOT / "sequences"  # the shared sequence sets

# Speed of one CPU core, for rough estimates before a run has measured its own: residues x
# steps per second. Measured with this force field on one core of a recent laptop: ~1.0
# million for single chains of 25-100 residues, 0.6-1.2 million for brushes (denser is
# slower). A cluster core may differ by a factor ~2 either way.
RESIDUE_STEPS_PER_SECOND = 1.0e6


# ── The list of runs: runs/runs.csv ───────────────────────────────────────────
REGISTRY = RUNS / "runs.csv"
REGISTRY_COLUMNS = ["number", "name", "description", "created", "launched", "finished", "proteins",
                    "physics_sets", "replicates", "tasks", "jobs", "done", "failed", "analysed",
                    "core_hours", "trajectories_gb", "last_checked"]


def registry():
    """The rows of runs/runs.csv: one per run ever made in this copy of noodle-md, with its
    number. A run folder without a row (made by an older version, or copied in) gets one.
    Numbers are never reused: a run's row stays when its folder is deleted."""
    rows = read_csv(REGISTRY) if REGISTRY.exists() else []
    known = {row["name"] for row in rows}
    folders = [p for p in RUNS.glob("*") if (p / "run.yaml").exists() and p.name not in known]
    for run_dir in sorted(folders, key=lambda p: p.stat().st_mtime):
        add_run(rows, run_dir)
    if folders:
        write_csv(rows, REGISTRY, REGISTRY_COLUMNS)
    return rows


def add_run(rows, run_dir):
    """Give a run folder the next number (a new row in `rows`); returns the number."""
    text = (run_dir / "run.yaml").read_text()
    created = re.search(r"\(created (\d{4}-\d\d-\d\d)\)", text)  # from its first line
    number = max((row["number"] for row in rows), default=0) + 1
    rows.append(dict(number=number, name=run_dir.name, description=(yaml.safe_load(text) or {}).get("description"),
                     created=created.group(1) if created else None))
    return number


def registry_row(rows, name):
    """The row of the run called `name` (the newest, if a deleted run had the same name)."""
    matches = [row for row in rows if row["name"] == name]
    return matches[-1] if matches else None


def update_registry(run_dir, **fields):
    rows = registry()
    registry_row(rows, run_dir.name).update(fields)
    write_csv(rows, REGISTRY, REGISTRY_COLUMNS)


def refresh_registry(run_dir, states=None):
    """Update the run's numbers in runs/runs.csv: tasks done, failed, analysed, core-hours..."""
    if run_dir.parent != RUNS:  # not a run of this folder (./noodle check's temporary one)
        return
    settings = yaml.safe_load((run_dir / "run.yaml").read_text()) or {}
    fields = dict(description=settings.get("description"), last_checked=datetime.now().isoformat(timespec="seconds"))
    if (run_dir / "tasks.csv").exists():
        tasks = read_tasks(run_dir)
        states = task_states(run_dir) if states is None else states
        statuses = [state["status"] for state in states.values()]
        reanalysed = run_dir / "analysis" / "reanalysed.csv"
        analysed = {i for i, state in states.items() if state.get("analysed_at")}
        analysed |= {row["task_id"] for row in read_csv(reanalysed)} if reanalysed.exists() else set()
        finished = [state["completed_at"] for state in states.values() if state.get("completed_at")]
        fields.update(
            proteins=len({t["sequence_name"] for t in tasks}), physics_sets=len({t["physics"] for t in tasks}),
            replicates=max(t["replicate"] for t in tasks), tasks=len(tasks), jobs=len(job_numbers(run_dir)),
            done=statuses.count("completed"), failed=statuses.count("failed"), analysed=len(analysed),
            core_hours=round(sum(state.get("run_time_s") or 0 for state in states.values()) / 3600, 2),
            trajectories_gb=round(sum(p.stat().st_size for p in (run_dir / "trajectories").glob("*.gsd")) / 1e9, 3),
            finished=max(finished) if statuses.count("completed") == len(tasks) else None)
    update_registry(run_dir, **fields)


# ── Choosing a run ────────────────────────────────────────────────────────────
def where_you_are():
    """The folder ./noodle was typed in (the wrapper passes it on, then moves to ROOT)."""
    return Path(os.environ.get("NOODLE_CWD", ".")).resolve()


def run_here():
    """The run folder you are in (or in a subfolder of), or None."""
    here = where_you_are()
    if here.is_relative_to(RUNS) and here != RUNS:
        run_dir = RUNS / here.relative_to(RUNS).parts[0]
        if (run_dir / "run.yaml").exists():
            return run_dir
    return None


def find_run(text):
    """The run folder that `text` stands for: a run name, a run number from runs/runs.csv,
    or a path to the run folder (what tab completes to). None if there is no such run."""
    if (RUNS / text / "run.yaml").exists():
        return RUNS / text
    path = (where_you_are() / text).resolve()
    if path.is_relative_to(RUNS) and path != RUNS:
        run_dir = RUNS / path.relative_to(RUNS).parts[0]
        if (run_dir / "run.yaml").exists():
            return run_dir
    if text.isdigit():
        rows = registry()
        row = next((row for row in rows if row["number"] == int(text)), None)
        if row and row is registry_row(rows, row["name"]) and (RUNS / row["name"] / "run.yaml").exists():
            return RUNS / row["name"]
        if row:  # (a later run may have the same name: it has a number of its own)
            sys.exit(f"Run {text} ({row['name']}) no longer exists: its folder was deleted.")
    return None


def list_runs():
    """Print every run with its number, how far along it is and its description."""
    rows = registry()
    rows = [row for row in rows if row is registry_row(rows, row["name"]) and (RUNS / row["name"] / "run.yaml").exists()]
    if not rows:
        print("There are no runs yet. Make one with: noodle new <run name>")
        return
    width = max(len(row["name"]) for row in rows)
    print(f"  {'#':>3}  {'run':{width}}  {'progress':24}  description")
    for row in rows:
        if not row.get("tasks"):
            state = "not launched yet"
        else:
            state = f"{row['done']} of {row['tasks']} done" + (f", {row['failed']} failed" if row.get("failed") else "")
        print(f"  {row['number']:>3}  {row['name']:{width}}  {state:24}  {row.get('description') or ''}")


def run_from_command_line(usage):
    """The run named on the command line (by name, number or path), or the run folder you
    are in when none is named. Otherwise print the usage and the runs, and stop."""
    args = sys.argv[1:]
    run_dir = find_run(args[0]) if len(args) == 1 else run_here() if not args else None
    if run_dir is None:
        if len(args) == 1:
            print(f"There is no run called '{args[0]}'.\n")
        else:
            print(f"Usage: {usage}   (a run's name or number, or nothing inside its folder)\n")
        list_runs()
        sys.exit(1)
    return run_dir


# ── CSV tables ────────────────────────────────────────────────────────────────
# Columns that always hold text, even when a value looks like a number (a protein named 171)
TEXT_COLUMNS = {"name", "description", "sequence_name", "sequence", "physics", "error", "analysis_error", "run",
                "slurm_job", "profile", "x_unit", "value_unit"}


def read_csv(path):
    """The rows of a CSV table as dicts. Values that look like numbers become numbers
    (except in TEXT_COLUMNS); empty cells become None."""
    def value(column, text):
        if text == "":
            return None
        if column in TEXT_COLUMNS or "_" in text:  # Python reads "1234_7" as the number 12347
            return text
        for kind in (int, float):
            try:
                return kind(text)
            except ValueError:
                pass
        return text

    with open(path, newline="") as f:
        return [{column: value(column, text) for column, text in row.items()} for row in csv.DictReader(f)]


def write_csv(rows, path, columns=None):
    """Write a list of dicts as a CSV table, in a single step (a temporary file that then
    replaces the old one), so a reader never sees half a file. The columns are `columns`,
    then any other key of the rows in the order they first appear. None stays empty;
    decimal numbers keep 6 significant digits."""
    def cell(v):
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return ""
        if isinstance(v, float):
            return f"{v:.6g}"
        return v

    columns = list(columns or [])
    for row in rows:
        columns += [key for key in row if key not in columns]
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows([cell(row.get(c)) for c in columns] for row in rows)
    os.replace(tmp, path)


# ── A run's tasks and jobs ────────────────────────────────────────────────────
def read_tasks(run_dir):
    """The rows of tasks.csv (task_id, sequence_name, physics, replicate, seed, length, sequence)."""
    return read_csv(run_dir / "tasks.csv")


def trajectory_path(run_dir, task):
    return run_dir / "trajectories" / f"task_{task['task_id']:04d}_{task['sequence_name']}.gsd"


def checkpoint_path(run_dir, task_id):
    return run_dir / "checkpoints" / f"task_{task_id:04d}.gsd"


# The columns of a job file, in order: who the task is, how it is going, and when it was
# analysed; the analysis results (observables.py) follow these
JOB_COLUMNS = ["task_id", "sequence_name", "physics", "replicate", "N", "chains", "status", "phase",
               "step", "total_steps", "steps_per_second", "run_time_s", "attempts", "frames", "host",
               "started_at", "last_update", "completed_at", "stop_reason", "error", "analysed_at",
               "analysis_error"]
# status of a task in a job file:
#   submitted  in a job that has not reached it yet
#   running    running now (or the job died without saying so: ./noodle status tells)
#   stopped    stopped itself when asked (time limit, scancel, Ctrl-C), after a checkpoint
#   failed     ended with an error (in the column error; the full message in the job's .out)
#   completed  done; analysed right away by the job (analysed_at, or analysis_error if that failed)


def job_file(run_dir, job):
    return run_dir / "jobs" / f"job_{job:04d}.csv"


def profiles_file(run_dir, job):
    return run_dir / "jobs" / f"job_{job:04d}_profiles.csv"


def job_numbers(run_dir):
    """The numbers of every job of the run, in order."""
    return sorted(int(p.stem.split("_")[1]) for p in (run_dir / "jobs").glob("job_[0-9][0-9][0-9][0-9].csv"))


def read_job(run_dir, job):
    return read_csv(job_file(run_dir, job))


def write_job(run_dir, job, rows):
    write_csv(rows, job_file(run_dir, job), JOB_COLUMNS)


def task_states(run_dir):
    """{task id: its current state}: its row in the newest job file that lists it, plus the
    job's number (key "job"). Tasks never submitted are left out."""
    states = {}
    for job in job_numbers(run_dir):
        for row in read_job(run_dir, job):
            states[row["task_id"]] = dict(row, job=job)
    return states


def read_submissions(run_dir):
    """{job number: SLURM job id}, e.g. {7: '1234_7'}."""
    path = run_dir / "jobs" / "submitted.csv"
    return {row["job"]: str(row["slurm_job"]) for row in read_csv(path)} if path.exists() else {}


# ── Physics conditions, brush size, run length and time ───────────────────────
# The physics settings copied into manifest.csv, so every row says under which conditions
# it was simulated, whatever the physics set is called (grafted-only ones stay empty for bulk).
CONDITION_SETTINGS = ["geometry", "temperature_K", "ionic_concentration_M", "graft_terminus",
                      "wall_charge_density", "wall_attraction", "grafting_density"]
GRAFTED_ONLY = ["graft_terminus", "wall_charge_density", "wall_attraction", "grafting_density"]


def conditions(physics):
    """{setting: value} of CONDITION_SETTINGS for one physics set, plus `chains`: how many
    chains each task simulates (1, or a brush's)."""
    grafted = physics["geometry"] == "grafted"
    row = {key: (physics.get(key) if grafted or key not in GRAFTED_ONLY else None) for key in CONDITION_SETTINGS}
    row["chains"] = brush_size(physics)[0]
    return row


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
    density = physics.get("grafting_density") or 0
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


def measured_speeds(tasks, states):
    """[(physics set, N, steps per second)] of every task of the run that has measured its
    speed on the machine it ran on: finished tasks (all their steps over all their run
    time) and started ones (their latest measurement)."""
    speeds = []
    for task in tasks:
        state = states.get(task["task_id"], {})
        if state.get("status") == "completed" and state.get("run_time_s") and state.get("total_steps"):
            speed = state["total_steps"] / state["run_time_s"]
        elif state.get("steps_per_second"):
            speed = state["steps_per_second"]
        else:
            continue
        speeds.append((task["physics"], task["length"], speed))
    return speeds


def seconds_left(task, state, physics, speeds, rough=False):
    """The run time a task still needs, in seconds, at its own measured speed, or else at
    that of the measured task closest in length (in the same physics set if there is one).
    With no measured speed at all: the rough core speed RESIDUE_STEPS_PER_SECOND if `rough`,
    else None."""
    n = task["length"]
    total = state.get("total_steps") or sum(run_lengths(physics, n)[:2])
    steps = total - (state.get("step") or 0)
    speed = state.get("steps_per_second")
    if not speed and speeds:
        similar = [s for s in speeds if s[0] == task["physics"]] or speeds
        speed = min(similar, key=lambda s: abs(s[1] - n))[2]
    if speed:
        return steps / speed
    return steps * brush_size(physics)[0] * n / RESIDUE_STEPS_PER_SECOND if rough else None


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
    """{job id: state} for this run's jobs in the SLURM queue, e.g.
    {'1234_17': 'RUNNING', '1234_18': 'PENDING'}; None on a computer without SLURM.

    A job array's elements are listed one by one ('<job>_<array index>'), and a job
    leaves the queue as soon as it ends, for whatever reason.
    """
    if shutil.which("squeue") is None:
        return None
    out = subprocess.run(["squeue", "--user", getpass.getuser(), "--noheader", "--array",
                          "--format=%i %T", f"--name={run_name}"], capture_output=True, text=True)
    return dict(line.split() for line in out.stdout.splitlines() if line.strip())
