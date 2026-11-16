"""./noodle status [run name]  --  how a run is going: done, running, waiting, stopped, failed.

With no name, it lists every run with how far along it is.

The report is printed and also saved as runs/<run name>/status.txt. On screen each list
shows at most SHOW_ON_SCREEN tasks; status.txt always has all of them.

A task's state comes from two places:
    its log (logs/task_0001.json)  how far it got, written by the task itself
    the SLURM queue                whether it is waiting or running right now
On a computer without SLURM, a task counts as running if its log was updated recently
(a running task updates it at least once an hour).
"""

import os
import sys
from collections import Counter, defaultdict
from datetime import datetime

import yaml

from runs import (RUNS, about, compress, duration, read_log, read_submissions, read_tasks,
                  run_from_command_line, run_lengths, slurm_queue)

SHOW_ON_SCREEN = 10  # tasks per list on screen (status.txt lists them all)
LIST_DONE_UP_TO = 20  # runs with this many tasks or fewer also list every finished task
STALE_HOURS = 1.5     # a running task that has not reported for this long is stuck or gone

# In the order the report shows them: what needs you first
STATES = ("failed", "stopped", "running", "waiting", "not launched", "done")
COLORS = {"failed": "31", "stopped": "33", "running": "36", "waiting": "2", "not launched": "2",
          "done": "32", "bold": "1", "dim": "2"}  # ANSI codes: red, yellow, cyan, dim, green


def hours_since(timestamp):
    return (datetime.now() - datetime.fromisoformat(timestamp)).total_seconds() / 3600


def task_state(log, slurm_job, queue):
    """The state of one task, with what the report needs to describe it:
        state    one of STATES
        percent  how far it got (None if it never started)
        silent   hours since it last reported, if a running task is overdue (else None)

    log        the task's log ({} if it has never started)
    slurm_job  the SLURM job id it was last submitted as (None if never submitted)
    queue      {job id: SLURM state} on a cluster; None on a computer without SLURM
    """
    info = dict(state=None, percent=None, silent=None)
    if log.get("total_steps"):
        info["percent"] = 100 * log["step"] / log["total_steps"]

    if log.get("status") == "completed":
        return dict(info, state="done", percent=100)

    slurm_state = queue.get(slurm_job) if queue else None  # None: not in the queue (any more)
    if slurm_state == "PENDING":
        return dict(info, state="waiting")
    if log.get("status") == "failed" and slurm_state is None:
        return dict(info, state="failed")
    if log.get("status") == "stopped" and slurm_state is None:
        return dict(info, state="stopped")  # it stopped itself (time limit, scancel, Ctrl-C)

    if queue is None:  # own computer: has it reported recently?
        if not log:
            return dict(info, state="not launched")
        alive = hours_since(log["last_update"]) < STALE_HOURS
    else:              # cluster: is its job still in the queue?
        if not log and slurm_job is None:
            return dict(info, state="not launched")
        alive = slurm_state is not None
        if alive and log and hours_since(log["last_update"]) > STALE_HOURS:
            info["silent"] = hours_since(log["last_update"])  # in the queue but not reporting

    return dict(info, state="running" if alive else "stopped")


# ── Time left, from speeds measured in this run ────────────────────────────────
def measured_speeds(rows):
    """[(physics set, N, steps per second)] for every task of this run that has measured
    its speed on the machine it ran on: finished tasks (all their steps over all their run
    time) and started ones (their latest hourly measurement)."""
    speeds = []
    for r in rows:
        log = r["log"]
        if r["state"] == "done" and log.get("run_time_s") and log.get("total_steps"):
            speed = log["total_steps"] / log["run_time_s"]
        elif log.get("steps_per_second"):
            speed = log["steps_per_second"]
        else:
            continue
        speeds.append((r["task"]["physics"], int(r["task"]["length"]), speed))
    return speeds


def add_time_left(rows, physics, speeds):
    """Give every unfinished, not failed task r["left"]: the run time it still needs, in
    seconds. It uses the task's own measured speed, or else that of the measured task
    closest in length (in the same physics set if there is one). None until some task of
    the run has measured its speed."""
    for r in rows:
        r["left"] = None
        if r["state"] in ("done", "failed") or not speeds:
            continue
        task, log = r["task"], r["log"]
        n = int(task["length"])
        total = log.get("total_steps") or sum(run_lengths(physics[task["physics"]], n)[:2])
        speed = log.get("steps_per_second")
        if not speed:
            similar = [s for s in speeds if s[0] == task["physics"]] or speeds
            speed = min(similar, key=lambda s: abs(s[1] - n))[2]
        r["left"] = (total - log.get("step", 0)) / speed


# ── The report ────────────────────────────────────────────────────────────────
class Report:
    """Builds the report's lines, in color for the screen or plain for status.txt."""

    def __init__(self, rows, color, limit):
        self.lines, self.color, self.limit = [], color, limit
        tasks = [r["task"] for r in rows]
        self.id_width = max(len(t["task_id"]) for t in tasks)
        self.name_width = min(24, max(len(t["sequence_name"]) for t in tasks))
        self.physics_width = max(len(t["physics"]) for t in tasks)

    def paint(self, text, color):
        return f"\033[{COLORS[color]}m{text}\033[0m" if self.color else text

    def add(self, *lines):
        self.lines += lines

    def task(self, r, details):
        t = r["task"]
        name = t["sequence_name"]
        if len(name) > self.name_width:
            name = name[:self.name_width - 1] + "…"
        self.add(f"  {t['task_id']:>{self.id_width}}  {name:{self.name_width}}  N={t['length']:<4}  "
                 f"{t['physics']:{self.physics_width}}  {details}")

    def tasks(self, rows, details):
        """One line per task, cut off on screen after self.limit lines (unless only a
        couple more would follow)."""
        shown = rows if self.limit is None or len(rows) <= self.limit + 2 else rows[:self.limit]
        for r in shown:
            self.task(r, details(r))
        if len(rows) > len(shown):
            self.add(self.paint(f"  … and {len(rows) - len(shown)} more (all are listed in status.txt)", "dim"))

    def heading(self, state, count, text=""):
        self.add("", self.paint(f"{state.upper()} {count}", state) + (f"  {text}" if text else ""))


def progress_bar(counts, total, width=30):
    """█ done, ▓ running, ░ everything else."""
    done = round(width * counts["done"] / total)
    running = round(width * (counts["done"] + counts["running"]) / total) - done
    return "█" * done + "▓" * running + "░" * (width - done - running)


def build_report(run_dir, settings, rows, speeds, color, limit):
    name = run_dir.name
    by_state = {s: [r for r in rows if r["state"] == s] for s in STATES}
    counts = {s: len(by_state[s]) for s in STATES}
    physics_sets = list(dict.fromkeys(r["task"]["physics"] for r in rows))
    launched = datetime.fromtimestamp((run_dir / "tasks.csv").stat().st_mtime)
    report = Report(rows, color, limit)

    # ── Summary ──
    if counts["done"] == len(rows):  # all done: from launch to the last task finishing
        finished = max((datetime.fromisoformat(r["log"]["completed_at"]) for r in rows
                        if "completed_at" in r["log"]), default=datetime.now())
        elapsed = (report.paint(f"took {duration((finished - launched).total_seconds())}", "bold")
                   + report.paint(f"   launched {launched:%Y-%m-%d %H:%M} · last task done {finished:%Y-%m-%d %H:%M}", "dim"))
    else:
        elapsed = (report.paint(f"{duration((datetime.now() - launched).total_seconds())} since launch", "bold")
                   + report.paint(f"   launched {launched:%Y-%m-%d %H:%M} · checked {datetime.now():%Y-%m-%d %H:%M}", "dim"))
    report.add(report.paint(name, "bold") + f"  {settings.get('description') or '(no description)'}",
               elapsed,
               "",
               f"{progress_bar(counts, len(rows))}  {counts['done']} of {len(rows)} tasks done "
               f"({100 * counts['done'] / len(rows):.0f}%)",
               "  " + " · ".join(report.paint(f"{counts[s]} {s}", s) for s in STATES if counts[s]))

    to_go = [r for r in rows if r["state"] not in ("done", "failed")]
    if to_go:
        if speeds:
            left = [r["left"] for r in to_go]
            report.add("", f"Time left  {about(sum(left))} of computing for the {len(to_go)} unfinished "
                           f"task{'s' * (len(to_go) > 1)}, the longest needs {about(max(left))}",
                       report.paint(f"           from the speeds measured by {len(speeds)} task"
                                    f"{'s' * (len(speeds) > 1)} of this run; how long you wait also "
                                    "depends on the queue", "dim"))
        else:
            report.add("", "Time left  not known yet: each task measures its speed a minute after it starts")

    if len(physics_sets) > 1:  # the same counts for each physics set, with done as done/total
        shown = [s for s in STATES if counts[s]]
        table = {}
        for physics in physics_sets:
            in_set = [r for r in rows if r["task"]["physics"] == physics]
            per_state = Counter(r["state"] for r in in_set)
            table[physics] = [f"{per_state[s]}/{len(in_set)}" if s == "done" else str(per_state[s] or "·")
                              for s in shown]
        widths = [max(7, len(s), *(len(cells[i]) for cells in table.values())) for i, s in enumerate(shown)]
        report.add("", "  " + " " * report.physics_width + "".join(f"  {s:>{w}}" for s, w in zip(shown, widths)))
        for physics, cells in table.items():
            report.add(f"  {physics:{report.physics_width}}" + "".join(f"  {c:>{w}}" for c, w in zip(cells, widths)))

    # ── Failed: grouped by error, since one cause often fails many tasks ──
    if counts["failed"]:
        report.heading("failed", counts["failed"], f"fix the cause first · full messages in runs/{name}/slurm/")
        by_error = defaultdict(list)
        for r in by_state["failed"]:
            by_error[r["log"].get("error", "(no error recorded)")].append(int(r["task"]["task_id"]))
        groups = sorted(by_error.items(), key=lambda g: -len(g[1]))
        for error, ids in groups if limit is None else groups[:limit]:
            report.add(f"  {len(ids)} × {error}", report.paint(f"      task{'s' * (len(ids) > 1)} {compress(sorted(ids)).replace(',', ', ')}", "dim"))
        if limit is not None and len(groups) > limit:
            report.add(report.paint(f"  … and {len(groups) - limit} other errors (all are listed in status.txt)", "dim"))

    # ── Stopped ──
    if counts["stopped"]:
        report.heading("stopped", counts["stopped"], f"./noodle launch {name} resumes them from their last checkpoint")
        report.tasks(by_state["stopped"], lambda r: "never started (see the slurm/ folder)" if r["percent"] is None else
                     f"stopped at {r['percent']:3.0f}% after {duration(r['log'].get('run_time_s', 0))} of run time"
                     + (f" ({r['log']['stop_reason']})" if r["log"].get("stop_reason") else "")
                     + (f", needs {about(r['left'])} more" if r["left"] is not None else ""))

    # ── Running ──
    if counts["running"]:
        report.heading("running", counts["running"])

        def running_details(r):
            if r["percent"] is None:
                return "starting"
            text = (f"{r['log'].get('phase', 'starting'):13s} {r['percent']:3.0f}%  "
                    f"{duration(r['log'].get('run_time_s', 0))} so far")
            if r["left"] is not None:
                text += f", {about(r['left'])} left"
            if r["silent"]:
                text += report.paint(f"  no report for {duration(r['silent'] * 3600)}: stuck?", "stopped")
            return text
        report.tasks(by_state["running"], running_details)

    # ── Waiting and not launched: just how many ──
    if counts["waiting"]:
        resuming = sum(r["percent"] is not None for r in by_state["waiting"])
        report.heading("waiting", counts["waiting"], "in the SLURM queue" +
                       (f" ({resuming} will resume from a checkpoint)" if resuming else ""))
    if counts["not launched"]:
        report.heading("not launched", counts["not launched"], f"./noodle launch {name} starts them")

    # ── Done ──
    if counts["done"]:
        done = sorted(by_state["done"], key=lambda r: r["log"]["run_time_s"])
        fastest, slowest = done[0], done[-1]
        total = sum(r["log"]["run_time_s"] for r in done)
        report.heading("done", counts["done"],
                       f"run time {duration(fastest['log']['run_time_s'])} (N={fastest['task']['length']}) "
                       f"to {duration(slowest['log']['run_time_s'])} (N={slowest['task']['length']}), "
                       f"{duration(total)} in total")
        if len(rows) <= LIST_DONE_UP_TO:
            report.tasks(done, lambda r: f"{r['log']['frames']} frames in {duration(r['log']['run_time_s'])}")

    # ── What to do next ──
    report.add("")
    if counts["done"] == len(rows):
        report.add(report.paint("All tasks are done.", "done") + f" Next: ./noodle analyze {name}")
    elif counts["stopped"] or counts["not launched"]:
        todo = [f"resumes the {counts['stopped']} stopped" if counts["stopped"] else "",
                f"starts the {counts['not launched']} not launched" if counts["not launched"] else ""]
        report.add(f"Next: ./noodle launch {name} " + " and ".join(t for t in todo if t)
                   + f" task{'s' * (counts['stopped'] + counts['not launched'] > 1)}"
                   + (" (and asks before retrying the failed ones)." if counts["failed"] else "."))
    elif counts["failed"] and not (counts["running"] or counts["waiting"]):
        report.add(f"Next: fix the errors above, then ./noodle launch {name} retries the failed tasks.")
    else:
        report.add(f"Nothing to do: check back later. ./noodle analyze {name} already works on the finished tasks."
                   if counts["done"] else "Nothing to do: check back later.")
    return report.lines


# ── Every run at a glance ─────────────────────────────────────────────────────
def list_all_runs():
    runs = sorted(p for p in RUNS.glob("*") if (p / "run.yaml").exists()) if RUNS.exists() else []
    if not runs:
        print("There are no runs yet. Make one with ./noodle new <run name>")
        return
    width = max(len(p.name) for p in runs)
    for run_dir in runs:
        settings = yaml.safe_load((run_dir / "run.yaml").read_text()) or {}
        if (run_dir / "tasks.csv").exists():
            tasks = read_tasks(run_dir)
            statuses = Counter(read_log(run_dir, t["task_id"]).get("status") for t in tasks)
            state = f"{statuses['completed']} of {len(tasks)} done"
            if statuses["failed"]:
                state += f", {statuses['failed']} failed"
        else:
            state = "not launched yet"
        print(f"  {run_dir.name:{width}}  {state:24s}  {settings.get('description') or ''}")
    print("\nDetails of one run: ./noodle status <run name>")


def main():
    if len(sys.argv) == 1:
        list_all_runs()
        return
    run_dir = run_from_command_line("./noodle status <run name>")
    settings = yaml.safe_load((run_dir / "run.yaml").read_text())
    if not (run_dir / "tasks.csv").exists():
        print(f"'{run_dir.name}' has not been launched yet: ./noodle launch {run_dir.name}")
        return

    queue = slurm_queue(run_dir.name)
    submitted = read_submissions(run_dir)
    rows = []
    for task in read_tasks(run_dir):
        log = read_log(run_dir, task["task_id"])
        rows.append(dict(task=task, log=log, **task_state(log, submitted.get(int(task["task_id"])), queue)))
    speeds = measured_speeds(rows)
    add_time_left(rows, yaml.safe_load((run_dir / "physics.yaml").read_text()), speeds)

    color = sys.stdout.isatty() and "NO_COLOR" not in os.environ
    print("\n".join(build_report(run_dir, settings, rows, speeds, color, SHOW_ON_SCREEN)))
    (run_dir / "status.txt").write_text("\n".join(build_report(run_dir, settings, rows, speeds, False, None)) + "\n")


if __name__ == "__main__":
    main()
