"""./launch-run <run name>  --  start (or continue) every unfinished task of a run.

The first time, it reads run.yaml + your CSV, checks the settings, shows the plan (how
many tasks, how long, how much disk) and asks before doing anything. If you say yes, it
writes the run's task list and starts the tasks:
    tasks.csv      one row per task: protein x physics set x replicate, with its seed
    physics.yaml   the physics sets this run uses, fully written out and frozen
Later launches show what is done, what is still in the queue and what will be
(re)started, then start only the tasks that are neither completed nor in the queue.
Tasks that were stopped resume from their last checkpoint.

Tasks run on a cluster (when sbatch exists) as one SLURM job array, one array element
per task; on your own computer, one after another, right here.
"""

import csv
import difflib
import math
import random
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from statistics import median

import yaml

from runs import (ROOT, about, read_log, read_submissions, read_tasks, run_from_command_line,
                  run_lengths, slurm_queue)

# The 20 standard amino acids: the same letters as the residue table in forcefield.py
# (repeated here so this script does not have to load HOOMD).
AMINO_ACIDS = set("ARNDCQEGHILKMFPSTWYV")

RUN_SETTINGS = ["description", "sequences_csv", "name_column", "sequence_column", "rows",
                "physics", "replicates", "time_limit_days", "max_running", "extra_sbatch"]
# The sets in physics.yaml that list every setting a bulk / grafted simulation needs
REFERENCE_SETS = {"bulk": "baseline", "grafted": "grafted"}
CHOICES = {"geometry": ["bulk", "grafted"], "graft_terminus": ["N", "C"],
           "wall_potential": ["lj93", "lj126"], "wall_attraction_scale": ["hydropathy", "uniform"]}

SPEED_TEST_SECONDS = 3      # length of the test simulation used to estimate run times
BYTES_PER_RESIDUE = 12      # trajectory size per residue per frame: one position (3 x 4 bytes)
BYTES_PER_FRAME = 100       # plus the step number and the file's index, per frame


# ── Checking the settings ─────────────────────────────────────────────────────
def did_you_mean(word, options):
    close = difflib.get_close_matches(word, options, n=1)
    return f" (did you mean '{close[0]}'?)" if close else ""


def stop_if(problems):
    if problems:
        print("Nothing was launched. Please fix:")
        for problem in problems:
            print(f"  - {problem}")
        sys.exit(1)


def check_run_settings(settings):
    """Clear mistakes in run.yaml: unknown or missing settings, impossible values."""
    problems = [f"run.yaml: unknown setting '{key}'{did_you_mean(key, RUN_SETTINGS)}"
                for key in settings if key not in RUN_SETTINGS]
    problems += [f"run.yaml: the setting '{key}' is missing"
                 for key in RUN_SETTINGS if key not in settings and key != "description"]
    stop_if(problems)
    if not isinstance(settings["replicates"], int) or settings["replicates"] < 1:
        problems.append("run.yaml: replicates must be a whole number, 1 or more")
    if not isinstance(settings["time_limit_days"], (int, float)) or settings["time_limit_days"] <= 0:
        problems.append("run.yaml: time_limit_days must be a number of days, e.g. 3 or 0.5")
    if not isinstance(settings["extra_sbatch"], list):
        problems.append('run.yaml: extra_sbatch must be a list, e.g. [] or ["--partition=normal"]')
    stop_if(problems)


def check_physics_set(name, physics, all_sets):
    """Clear mistakes in a physics set: misspelled or missing settings, impossible choices."""
    geometry = physics.get("geometry")
    if geometry not in CHOICES["geometry"]:
        return [f"physics set '{name}': geometry must be bulk or grafted, not '{geometry}'"]
    every_setting = resolve_physics(all_sets, REFERENCE_SETS["grafted"])
    needed = resolve_physics(all_sets, REFERENCE_SETS[geometry])
    problems = [f"physics set '{name}': unknown setting '{key}'{did_you_mean(key, every_setting)}"
                for key in physics if key not in every_setting]
    problems += [f"physics set '{name}': the setting '{key}' is missing"
                 for key in needed if key not in physics]
    for key, options in CHOICES.items():
        if key in physics and physics[key] not in options:
            problems.append(f"physics set '{name}': {key} must be {' or '.join(options)}, not '{physics[key]}'")
    return problems


# ── Building the task list (first launch only) ────────────────────────────────
def select_rows(rows, setting):
    """rows: all | 7 | 101-200 (counted from 1)."""
    text = str(setting).strip()
    if text == "all":
        return rows
    match = re.fullmatch(r"(\d+)(?:\s*-\s*(\d+))?", text)
    if not match:
        sys.exit(f"rows: '{setting}' is not understood. Use all, a number like 7, or a range like 101-200.")
    first = int(match.group(1))
    last = int(match.group(2) or first)
    if not 1 <= first <= last <= len(rows):
        sys.exit(f"rows: {text} is outside the file, which has {len(rows)} rows.")
    return rows[first - 1:last]


def read_sequences(run_dir, settings):
    """({name: sequence}, [names skipped]) for the selected rows of the run's CSV."""
    csv_path = run_dir / settings["sequences_csv"]
    if not csv_path.exists():
        sys.exit(f"Can't find the sequence file {csv_path}. Copy it into runs/{run_dir.name}/ "
                 "or fix sequences_csv in run.yaml.")
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    name_col, seq_col = settings["name_column"], settings["sequence_column"]
    for column in (name_col, seq_col):
        if column not in rows[0]:
            sys.exit(f"{csv_path.name} has no column '{column}'. Its columns are: {list(rows[0])}")

    sequences, skipped = {}, []
    for row in select_rows(rows, settings["rows"]):
        # names end up in file names, so keep them to safe characters
        name = re.sub(r"[^A-Za-z0-9._-]+", "_", row[name_col].strip())
        sequence = row[seq_col].strip().upper()
        if name in sequences:
            sys.exit(f"The name '{name}' appears twice in {csv_path.name}; names must be unique.")
        if not sequence or set(sequence) - AMINO_ACIDS:
            skipped.append(name)
            continue
        sequences[name] = sequence
    return sequences, skipped


def resolve_physics(all_sets, name):
    """A physics set with everything it `extends` filled in underneath it."""
    if name not in all_sets:
        sys.exit(f"physics.yaml has no set called '{name}'. Its sets are: {list(all_sets)}")
    settings = dict(all_sets[name])
    parent = settings.pop("extends", None)
    return {**resolve_physics(all_sets, parent), **settings} if parent else settings


def plan_tasks(run_dir, settings):
    """Everything the first launch would create, checked, but not written yet.

    Returns (sequences, skipped, physics, tasks): the proteins, the names of the rows
    that were skipped, the physics sets {name: settings}, and the task list."""
    sequences, skipped = read_sequences(run_dir, settings)
    if not sequences:
        sys.exit("No usable sequences were selected.")
    physics_names = settings["physics"] if isinstance(settings["physics"], list) else [settings["physics"]]
    all_sets = yaml.safe_load((ROOT / "physics.yaml").read_text())
    physics = {name: resolve_physics(all_sets, name) for name in physics_names}

    problems = []
    for name, settings_of_set in physics.items():
        problems += check_physics_set(name, settings_of_set, all_sets)
    stop_if(problems)
    for name, settings_of_set in physics.items():
        try:
            run_lengths(settings_of_set, min(map(len, sequences.values())))
        except ValueError as error:
            problems.append(f"physics set '{name}': {error}")
    stop_if(problems)

    tasks = []
    for name, sequence in sequences.items():
        for physics_name in physics_names:
            for replicate in range(1, settings["replicates"] + 1):
                tasks.append(dict(task_id=len(tasks) + 1, sequence_name=name, physics=physics_name,
                                  replicate=replicate,
                                  seed=random.randint(1, 65535),  # HOOMD seeds are 16-bit numbers
                                  length=len(sequence), sequence=sequence))
    return sequences, skipped, physics, tasks


def write_tasks(run_dir, physics, tasks):
    """Freeze the run: write its physics sets and task list."""
    (run_dir / "physics.yaml").write_text(
        "# The physics sets this run uses, copied from physics.yaml when it was launched.\n"
        + yaml.safe_dump(physics, sort_keys=False))
    with open(run_dir / "tasks.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(tasks[0]))
        writer.writeheader()
        writer.writerows(tasks)


# ── The plan shown before the first launch ────────────────────────────────────
def describe_physics(p):
    """One line saying what a physics set simulates."""
    if p["geometry"] == "bulk":
        where = "free chain in solution"
    else:
        attraction = ("repulsive" if p["wall_attraction"] == 0 else
                      f"attractive {p['wall_attraction']:g} ({p['wall_attraction_scale']})")
        charge = "neutral" if p["wall_charge_density"] == 0 else f"charge {p['wall_charge_density']:+g} e/nm^2"
        where = f"grafted by the {p['graft_terminus']}-terminus to a {attraction}, {charge} {p['wall_potential']} surface"
    return f"{where}; {p['temperature_K']:g} K, {1000 * p['ionic_concentration_M']:g} mM salt"


def in_words(number):
    """1_800_000 -> '1.8 million'."""
    for size, word in ((1e9, "billion"), (1e6, "million"), (1e3, "thousand")):
        if number >= size:
            return f"{number / size:.3g} {word}"
    return str(number)


def file_size(n_bytes):
    for size, unit in ((1e12, "TB"), (1e9, "GB"), (1e6, "MB")):
        if n_bytes >= size:
            return f"{n_bytes / size:.1f} {unit}"
    return f"{n_bytes / 1e3:.0f} kB"


def measure_speed(sequence, physics):
    """Steps per second for this sequence on this computer, from a few seconds of simulation.

    HOOMD is imported only here, so the rest of ./launch-run stays quick.
    """
    from simulate import add_forces, create_simulation, starting_snapshot, system_geometry

    sim = create_simulation(seed=1)
    geometry = system_geometry(sequence, physics)
    sim.create_state_from_snapshot(starting_snapshot(sequence, physics, geometry, seed=1))
    add_forces(sim, physics, geometry)
    sim.run(1000)  # let the starting chain relax a little first
    start, start_step = time.monotonic(), sim.timestep
    while time.monotonic() - start < SPEED_TEST_SECONDS:
        sim.run(1000)
    return (sim.timestep - start_step) / (time.monotonic() - start)


def show_plan(run_dir, settings, sequences, skipped, physics, tasks, on_cluster):
    lengths = sorted(map(len, sequences.values()))
    print(f"New run: {run_dir.name}" + (f'   "{settings["description"]}"' if settings.get("description") else ""))
    print()
    sizes = (f"{lengths[0]} residues" if lengths[0] == lengths[-1] else
             f"{lengths[0]} to {lengths[-1]} residues (median {median(lengths):.0f})")
    print(f"Sequences  {settings['sequences_csv']}, rows {settings['rows']}: {len(sequences)} proteins, {sizes}")
    if skipped:
        print(f"           skipped {len(skipped)} with letters other than the 20 standard amino acids: "
              + ", ".join(skipped[:5]) + (", ..." if len(skipped) > 5 else ""))

    print("Physics")
    for name, p in physics.items():
        print(f"  {name:16s} {describe_physics(p)}")
        if name.startswith("test_"):
            print(f"  {'':16s} [!] a quick-test set: tiny run length, results are NOT physically meaningful")
    print(f"Tasks      {len(sequences)} proteins x {len(physics)} physics set(s) x "
          f"{settings['replicates']} replicate(s) = {len(tasks)} tasks")

    print("Length     steps per task (equilibration + production), for the shortest and longest protein")
    for name, p in physics.items():
        shortest, longest = (sum(run_lengths(p, n)[:2]) for n in (lengths[0], lengths[-1]))
        steps = (in_words(shortest) if shortest == longest else
                 f"{in_words(shortest)} (N={lengths[0]}) to {in_words(longest)} (N={lengths[-1]})")
        print(f"  {name:16s} {steps}, saving {p['frames']} frames")

    # Time: one short test simulation of the longest protein. The cost of a step grows
    # about in proportion to the number of residues, which gives every other protein's speed.
    print(f"Time       (measuring the speed with a {SPEED_TEST_SECONDS}-second test simulation...)", flush=True)
    longest_name = max(sequences, key=lambda name: len(sequences[name]))
    first_set = next(iter(physics.values()))
    speed_times_n = measure_speed(sequences[longest_name], first_set) * lengths[-1]
    seconds = [sum(run_lengths(physics[t["physics"]], t["length"])[:2]) * t["length"] / speed_times_n
               for t in tasks]
    print(f"           {about(min(seconds))} to {about(max(seconds))} per task, "
          f"{about(sum(seconds))} of computing in total")
    print("           (a rough estimate, measured on this computer: cluster nodes may be faster or slower)")
    if on_cluster:
        limit = settings["time_limit_days"] * 86400
        too_long = [s for s in seconds if s > limit]
        if too_long:
            print(f"           [!] {len(too_long)} task(s) need longer than the {settings['time_limit_days']:g}-day time "
                  f"limit: they will stop and need ./launch-run again ({math.ceil(max(seconds) / limit) - 1} "
                  "time(s) for the longest), unless you raise time_limit_days now")
    else:
        print("           here, without SLURM, the tasks run one after another, so the total is the waiting time")

    n_bytes = sum(int(physics[t["physics"]]["frames"]) * (BYTES_PER_RESIDUE * t["length"] + BYTES_PER_FRAME)
                  for t in tasks)
    print(f"Disk       ~{file_size(n_bytes)} of trajectories")

    if on_cluster:
        at_once = f"at most {settings['max_running']} at a time" if settings["max_running"] else "no limit on how many at a time"
        extra = f", extra sbatch options: {' '.join(settings['extra_sbatch'])}" if settings["extra_sbatch"] else ""
        print(f"Cluster    {len(tasks)} SLURM jobs, time limit {settings['time_limit_days']:g} day(s) each, {at_once}{extra}")
    else:
        print("Where      no SLURM here: the tasks run in this terminal (Ctrl-C stops, ./launch-run continues)")
    print()
    print("Once launched, the run is frozen: only time_limit_days, max_running and extra_sbatch can still change.")


def confirm(question):
    """Ask a yes/no question; anything but y or yes means no."""
    try:
        return input(f"{question} [y/N] ").strip().lower() in ("y", "yes")
    except EOFError:  # no keyboard (e.g. input from a pipe): treat as no
        print()
        return False


# ── Deciding what to start (later launches) ───────────────────────────────────
def plan_relaunch(run_dir, on_cluster):
    """Sort the tasks into done / still in the queue / failed / to start, say so, and
    return the ids of the tasks to start (failed ones only if you agree)."""
    tasks = read_tasks(run_dir)
    names = {int(t["task_id"]): t["sequence_name"] for t in tasks}
    logs = {task_id: read_log(run_dir, task_id) for task_id in names}
    done = [i for i in names if logs[i].get("status") == "completed"]
    if len(done) == len(tasks):
        print(f"All {len(tasks)} tasks of '{run_dir.name}' are completed.")
        return []

    unfinished = [i for i in names if i not in done]
    in_queue = []
    if on_cluster:
        queue, submitted = slurm_queue(run_dir.name), read_submissions(run_dir)
        in_queue = [i for i in unfinished if submitted.get(i) in queue]
    failed = [i for i in unfinished if i not in in_queue and logs[i].get("status") == "failed"]
    to_start = [i for i in unfinished if i not in in_queue and i not in failed]
    resuming = [i for i in to_start if (run_dir / "checkpoints" / f"task_{i:04d}.gsd").exists()]

    parts = [f"{len(done)} done"]
    if in_queue:
        parts.append(f"{len(in_queue)} in the queue (left alone)")
    if failed:
        parts.append(f"{len(failed)} failed")
    parts.append(f"{len(to_start)} to start ({len(resuming)} resume from a checkpoint, "
                 f"{len(to_start) - len(resuming)} from the beginning)")
    print(f"{len(tasks)} tasks: " + " · ".join(parts))

    if failed:
        print(f"\n{len(failed)} task(s) failed with an error last time:")
        for i in failed[:10]:
            print(f"  task {i:>4}  {names[i]}: {logs[i]['error']}")
        if len(failed) > 10:
            print(f"  ... and {len(failed) - 10} more (./status {run_dir.name} lists them all)")
        print("They will most likely fail again unless the cause has been fixed"
              + (f" (full messages in runs/{run_dir.name}/slurm/)." if on_cluster else "."))
        if confirm("Start the failed tasks again anyway?"):
            to_start = sorted(to_start + failed)
    return to_start


# ── Starting the tasks ────────────────────────────────────────────────────────
def compress(ids):
    """[1, 2, 3, 7] -> '1-3,7', the way sbatch --array wants it."""
    parts, start = [], ids[0]
    for previous, current in zip(ids, ids[1:] + [None]):
        if current != previous + 1:
            parts.append(f"{start}-{previous}" if previous != start else str(start))
            start = current
    return ",".join(parts)


def max_array_size():
    out = subprocess.run(["scontrol", "show", "config"], capture_output=True, text=True).stdout
    match = re.search(r"MaxArraySize\s*=\s*(\d+)", out)
    return int(match.group(1)) if match else 1001


def submit_to_slurm(run_dir, settings, task_ids):
    name = run_dir.name
    (run_dir / "slurm").mkdir(exist_ok=True)
    time_limit_minutes = round(settings["time_limit_days"] * 24 * 60)
    jobs = {}  # task id -> SLURM job id

    # Array indices must stay below the cluster's MaxArraySize, so a long task list is
    # sent as several arrays, each with an offset: task id = offset + array index.
    chunk = max_array_size() - 1
    for offset in range(0, max(task_ids), chunk):
        ids = [i for i in task_ids if offset < i <= offset + chunk]
        if not ids:
            continue
        array = compress([i - offset for i in ids])
        if settings["max_running"]:
            array += f"%{settings['max_running']}"
        output = "task_%a.out" if offset == 0 else f"task_%a_plus_{offset}.out"
        command = ["sbatch", "--parsable", f"--job-name={name}", f"--array={array}",
                   f"--time={time_limit_minutes}", f"--chdir={ROOT}", "--open-mode=append",
                   f"--output={run_dir / 'slurm' / output}", *settings["extra_sbatch"],
                   str(ROOT / "src" / "slurm_job.sh"), f"runs/{name}", str(offset)]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            sys.exit(f"sbatch refused the job:\n{result.stderr}"
                     + (f"\n({len(jobs)} task(s) had already been submitted.)" if jobs else ""))
        job_id = result.stdout.strip().split(";")[0]
        jobs.update({i: f"{job_id}_{i - offset}" for i in ids})

        # Record which SLURM job runs which task, so ./status can find it in the queue.
        submissions = run_dir / "slurm" / "submissions.csv"
        is_new_file = not submissions.exists()
        with open(submissions, "a", newline="") as f:
            writer = csv.writer(f)
            if is_new_file:
                writer.writerow(["task_id", "slurm_job", "submitted_at"])
            submitted_at = datetime.now().isoformat(timespec="seconds")
            writer.writerows([i, jobs[i], submitted_at] for i in ids)
        print(f"Submitted {len(ids)} task(s) as SLURM job {job_id}.")

    # Check that the jobs really are in the queue
    queue = slurm_queue(name) or {}
    found = sum(job in queue for job in jobs.values())
    if found == len(jobs):
        print(f"Confirmed: {found} of {len(jobs)} task(s) are in the SLURM queue.")
    else:
        print(f"Only {found} of {len(jobs)} task(s) are in the SLURM queue; the others may have "
              f"finished or failed already: see ./status {name}")
    print(f"Check on them any time with:  ./status {name}")


def run_here(run_dir, task_ids):
    print(f"No SLURM here: running {len(task_ids)} task(s) one after another. "
          f"Ctrl-C stops; ./launch-run {run_dir.name} again resumes.\n")
    for n, task_id in enumerate(task_ids, 1):
        print(f"── task {task_id} ({n} of {len(task_ids)}) ──")
        result = subprocess.run([sys.executable, str(ROOT / "src" / "simulate.py"), str(run_dir), str(task_id)])
        if result.returncode != 0:
            print(f"Task {task_id} FAILED (error above).")
        print()
    failed = [i for i in task_ids if read_log(run_dir, i).get("status") == "failed"]
    print(f"Finished: {len(task_ids) - len(failed)} of {len(task_ids)} task(s) completed"
          + (f", {len(failed)} failed (tasks {compress(failed)})." if failed else "."))


def main():
    run_dir = run_from_command_line("./launch-run <run name>")
    settings = yaml.safe_load((run_dir / "run.yaml").read_text())
    check_run_settings(settings)
    on_cluster = shutil.which("sbatch") is not None

    if not (run_dir / "tasks.csv").exists():
        sequences, skipped, physics, tasks = plan_tasks(run_dir, settings)
        show_plan(run_dir, settings, sequences, skipped, physics, tasks, on_cluster)
        if not confirm(f"Launch these {len(tasks)} tasks?"):
            print("Nothing was launched.")
            return
        write_tasks(run_dir, physics, tasks)
        to_start = [t["task_id"] for t in tasks]
    else:
        if (run_dir / "run.yaml").stat().st_mtime > (run_dir / "tasks.csv").stat().st_mtime:
            print("Note: run.yaml was edited after this run was launched. Only time_limit_days, "
                  "max_running and extra_sbatch changes take effect; for anything else make a new run.\n")
        to_start = plan_relaunch(run_dir, on_cluster)

    if not to_start:
        return
    if on_cluster:
        submit_to_slurm(run_dir, settings, to_start)
    else:
        run_here(run_dir, to_start)


if __name__ == "__main__":
    main()
