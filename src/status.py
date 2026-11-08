"""./status <run name>  --  everything about one run: done, running, waiting, stopped, failed.

The report is printed and also saved as runs/<run name>/status.txt.

A task's state comes from two places:
    its log (logs/task_0001.json)  how far it got, written by the task itself
    the SLURM queue                whether it is waiting or running right now
On a computer without SLURM, a task counts as running if its log was updated recently
(a running task updates it at least once an hour).
"""

from datetime import datetime

import yaml

from runs import about, duration, read_log, read_submissions, read_tasks, run_from_command_line, slurm_queue

LIST_ALL_UP_TO = 20  # runs with this many tasks or fewer list every task
STALE_HOURS = 1.5    # without SLURM: a running task that has not reported for this long has stopped


def hours_since(timestamp):
    return (datetime.now() - datetime.fromisoformat(timestamp)).total_seconds() / 3600


def task_state(log, slurm_job, queue):
    """(state, details) of one task.

    log        the task's log ({} if it has never started)
    slurm_job  the SLURM job id it was last submitted as (None if never submitted)
    queue      {job id: SLURM state} on a cluster; None on a computer without SLURM
    """
    if log.get("status") == "completed":
        return "done", f"{log['frames']} frames, {duration(log['run_time_s'])} of run time"

    slurm_state = queue.get(slurm_job) if queue else None  # None: not in the queue (any more)
    if slurm_state == "PENDING":
        return "waiting", "in the SLURM queue"
    if log.get("status") == "failed" and slurm_state is None:
        return "failed", log["error"]

    if queue is None:  # own computer: has it reported recently?
        if not log:
            return "not launched", ""
        alive = hours_since(log["last_update"]) < STALE_HOURS
    else:              # cluster: is its job still in the queue?
        if not log and slurm_job is None:
            return "not launched", ""
        alive = slurm_state is not None

    if not log:  # submitted, but the task itself never got going
        if alive:
            return "running", "starting"
        return "stopped", "its job ended before the task started (see the slurm/ folder)"

    percent = 100 * log["step"] / log["total_steps"]
    speed = log.get("steps_per_second")
    left = about((log["total_steps"] - log["step"]) / speed) if speed else "an unknown time"
    if alive:
        return "running", f"{log.get('phase', 'starting')} {percent:.0f}%, {left} left"
    return "stopped", (f"stopped at {percent:.0f}% after {duration(log.get('run_time_s', 0))}, "
                       f"needs {left} more: ./launch-run resumes it")


def main():
    run_dir = run_from_command_line("./status <run name>")
    settings = yaml.safe_load((run_dir / "run.yaml").read_text())
    if not (run_dir / "tasks.csv").exists():
        print(f"'{run_dir.name}' has not been launched yet: ./launch-run {run_dir.name}")
        return

    tasks = read_tasks(run_dir)
    queue = slurm_queue(run_dir.name)
    submitted = read_submissions(run_dir)
    rows = []
    for task in tasks:
        task_id = int(task["task_id"])
        log = read_log(run_dir, task_id)
        state, details = task_state(log, submitted.get(task_id), queue)
        rows.append(dict(task=task, log=log, state=state, details=details))

    counts = {s: sum(r["state"] == s for r in rows)
              for s in ("done", "running", "waiting", "stopped", "failed", "not launched")}
    done_fraction = counts["done"] / len(rows)
    bar = "█" * round(20 * done_fraction) + "░" * (20 - round(20 * done_fraction))
    physics = settings["physics"] if isinstance(settings["physics"], list) else [settings["physics"]]
    launched = datetime.fromtimestamp((run_dir / "tasks.csv").stat().st_mtime)

    lines = [f"{run_dir.name}: {settings.get('description') or '(no description)'}",
             f"physics: {', '.join(physics)}   launched {launched:%Y-%m-%d %H:%M}",
             "",
             f"{len(rows)} tasks  {bar}  " + " · ".join(f"{n} {s}" for s, n in counts.items() if n),
             ""]

    def task_line(r):
        t = r["task"]
        return (f"  task {t['task_id']:>4}  {t['sequence_name'][:24]:24s} N={t['length']:<5} "
                f"{t['physics'][:14]:14s} {r['state']:8s} {r['details']}")

    if len(rows) <= LIST_ALL_UP_TO:
        lines += [task_line(r) for r in rows]
    else:
        unfinished = [r for r in rows if r["state"] != "done"]
        if unfinished:
            lines.append("NOT FINISHED")
            # failed and stopped tasks first: those are the ones that need you
            order = ["failed", "stopped", "running", "waiting", "not launched"]
            lines += [task_line(r) for r in sorted(unfinished, key=lambda r: order.index(r["state"]))]
        done = [r for r in rows if r["state"] == "done"]
        if done:
            fastest = min(done, key=lambda r: r["log"]["run_time_s"])
            slowest = max(done, key=lambda r: r["log"]["run_time_s"])
            lines += ["", f"FINISHED  {len(done)} tasks, run time from "
                          f"{duration(fastest['log']['run_time_s'])} (N={fastest['task']['length']}) to "
                          f"{duration(slowest['log']['run_time_s'])} (N={slowest['task']['length']})"]

    if counts["failed"]:
        lines += ["", f"{counts['failed']} task(s) failed with the error shown (full message in "
                      f"runs/{run_dir.name}/slurm/ on the cluster). Fix the cause before relaunching them."]
    if counts["stopped"]:
        lines += ["", f"{counts['stopped']} task(s) stopped. ./launch-run {run_dir.name} resumes them "
                      "(tasks still in the queue are left alone)."]
    if counts["done"] == len(rows):
        lines += ["", "All tasks are done."]

    report = "\n".join(lines)
    print(report)
    (run_dir / "status.txt").write_text(report + "\n")


if __name__ == "__main__":
    main()
