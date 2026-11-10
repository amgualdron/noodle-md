"""./noodle calibrate setup|report|cleanup [config]  --  the brush calibration experiment.

How long must a 10 x 10 brush run? It needs an equilibration long enough to forget its
starting configuration, and a production long enough to give every chain enough
independent conformations. Both scale with the chains' relaxation time, which depends on
the chain length N, the grafting density and the surface. This experiment measures them,
together with how fast the cluster simulates brushes, and writes a self-contained report
with recommended run lengths and costs. The settings are in
experiments/brush_calibration.yaml (or the config file named on the command line).

    setup    picks the sequences, writes runs/<run_name>/ (run.yaml, the sequences, its
             own physics sets in physics_sets.yaml, a frozen copy of the config). Then
             ./noodle launch <run_name> runs it like any other run.
    report   for every completed task: computes what the report needs from its trajectory
             (small CSV/JSON files in results_dir), then deletes the trajectory and any
             checkpoint of that task. Then (re)writes results_dir/report.md from every
             task analysed so far. Safe to run any time, as often as you like.
    cleanup  runs report, then deletes every .gsd file the run still holds (those of
             unfinished or failed tasks), after you type the run name.

Every task records its brush from step 0 (no separate equilibration), so the report
sees it relax from the starting configuration and then fluctuate around equilibrium.
Time is in HOOMD time units (0.489 ps; a step is dt of them), and is also given divided
by N^2.2: those coefficients are exactly the rouse_multiplier / production_multiplier
of physics.yaml (steps = multiplier x N^2.2 / dt).
"""

import json
import math
import socket
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import gsd.fl
import numpy as np
import pandas as pd
import yaml

from analyze_run import trajectory_metrics
from launch_run import check_physics_set, physics_sets, resolve_physics
from runs import ROOT, RUNS, SEQUENCES, brush_size, file_size, read_log, read_tasks, run_lengths

DEFAULT_CONFIG = ROOT / "experiments" / "brush_calibration.yaml"
FROZEN_CONFIG = "calibration.yaml"          # the config, copied into the run folder by setup
SEQUENCE_FILE = "calibration_sequences.csv"  # the picked sequences, in the run folder
MARKER = "RUN.json"                          # in results_dir: which run its files belong to
EXPONENT = 2.2                               # run lengths scale like N^2.2 (physics.yaml)
PS_PER_TIME_UNIT = 0.489
EQUILIBRATION_OBSERVABLES = ["com_height", "ree_z", "rg"]
CHAIN_OBSERVABLES = ["rg", "ree_z", "com_height", "ree_vector"]  # per-chain relaxation times
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
# what each kind of sequence aims for: quantiles (within the sequences of that length) of
# hydropathy, fraction of charged residues and net charge per residue
KINDS = {"typical": dict(hydropathy=0.5, fcr=0.5, ncpr=0.5),
         "hydrophobic": dict(hydropathy=0.9, fcr=0.5, ncpr=0.5),
         "basic": dict(hydropathy=0.5, ncpr=0.9),
         "acidic": dict(hydropathy=0.5, ncpr=0.1)}
REQUIRED = ["run_name", "description", "results_dir", "sequences_csv", "lengths", "kinds", "base_set",
            "chains", "conditions", "frames", "replicates", "equilibration_safety",
            "target_independent_samples", "time_limit_days", "max_running", "extra_sbatch"]


# ── Config ────────────────────────────────────────────────────────────────────
def load_config(path):
    config = yaml.safe_load(Path(path).read_text())
    missing = [key for key in REQUIRED if key not in config]
    unknown = [kind for kind in config.get("kinds", []) if kind not in KINDS]
    if missing or unknown:
        sys.exit(f"{path}: " + "; ".join(([f"missing {', '.join(missing)}"] if missing else [])
                                         + ([f"unknown kinds {unknown} (use {list(KINDS)})"] if unknown else [])))
    return config


def condition_name(condition):
    """cal_s0.03_neutral, cal_s0.1_q-0.1: the physics set of one condition."""
    charge = condition["wall_charge_density"]
    return f"cal_s{condition['grafting_density']:g}_" + ("neutral" if charge == 0 else f"q{charge:+g}")


def git_version():
    def git(*args):
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    commit = git("rev-parse", "--short", "HEAD") or "unknown"
    return commit + (" (with uncommitted changes)" if git("status", "--porcelain", "--", "src", "physics.yaml") else "")


# ── setup ─────────────────────────────────────────────────────────────────────
def sequence_properties(sequence):
    from forcefield import AA_PARAMS  # imports HOOMD; only setup needs it
    n = len(sequence)
    charges = [AA_PARAMS[aa]["charge"] for aa in sequence]
    return dict(N=n, hydropathy=sum(AA_PARAMS[aa]["hydropathy"] for aa in sequence) / n,
                fcr=sum(abs(q) == 1 for q in charges) / n, ncpr=sum(charges) / n)


def pick_sequences(config):
    """One sequence per length and kind (see the config file for what the kinds mean)."""
    path = SEQUENCES / config["sequences_csv"]
    table = pd.read_csv(path)
    table = table[table["sequence"].str.fullmatch(f"[{AMINO_ACIDS}]+", na=False)].reset_index(drop=True)
    properties = pd.DataFrame([sequence_properties(s) for s in table["sequence"]])
    table = pd.concat([table.drop(columns=properties.columns, errors="ignore"), properties], axis=1)
    picked = []
    for length in config["lengths"]:
        tolerance = max(2, round(0.03 * length))
        pool = table[(table["N"] - length).abs() <= tolerance]
        if "tier" in pool and (pool["tier"] == "confident").sum() >= 10:
            pool = pool[pool["tier"] == "confident"]
        if len(pool) < len(config["kinds"]):
            sys.exit(f"Only {len(pool)} sequences of {length} +- {tolerance} residues in {path.name}.")
        spread = pool[["hydropathy", "fcr", "ncpr"]].std().replace(0, 1)
        for kind in config["kinds"]:
            score = 0.1 * (pool["N"] - length).abs() / tolerance  # prefer the exact length a little
            for prop, quantile in KINDS[kind].items():
                score = score + ((pool[prop] - pool[prop].quantile(quantile)) / spread[prop]) ** 2
            score = score[~pool["id"].isin([p["source_id"] for p in picked])]
            row = pool.loc[score.idxmin()]
            picked.append(dict(id=f"L{length}_{kind}_{row['id']}", sequence=row["sequence"], kind=kind,
                               target_length=length, N=row["N"], hydropathy=round(row["hydropathy"], 4),
                               fcr=round(row["fcr"], 4), ncpr=round(row["ncpr"], 4), source_id=row["id"],
                               source_csv=path.name))
    return pd.DataFrame(picked)


def calibration_sets(config):
    sets = {}
    for condition in config["conditions"]:
        charge = condition["wall_charge_density"]
        sets[condition_name(condition)] = dict(
            extends=config["base_set"],
            description=(f"Brush calibration: {math.isqrt(config['chains'])} x {math.isqrt(config['chains'])} at "
                         f"{condition['grafting_density']:g} chains/nm^2, "
                         + ("neutral surface" if charge == 0 else f"surface {charge:+g} e/nm^2")
                         + ", recorded from step 0"),
            chains=config["chains"], grafting_density=condition["grafting_density"],
            wall_charge_density=charge, equilibration_steps=0,
            production_multiplier=condition["run_multiplier"], frames=config["frames"])
    return sets


def setup(config_path):
    config = load_config(config_path)
    run_dir = RUNS / config["run_name"]
    if (run_dir / "tasks.csv").exists():
        sys.exit(f"runs/{config['run_name']} has already been launched. For a new experiment, change "
                 f"run_name (and results_dir) in {config_path}.")
    sequences = pick_sequences(config)
    sets = calibration_sets(config)

    run_dir.mkdir(parents=True, exist_ok=True)
    sequences.to_csv(run_dir / SEQUENCE_FILE, index=False)
    (run_dir / "physics_sets.yaml").write_text(
        "# The physics sets of the brush calibration, written by ./noodle calibrate setup.\n"
        + yaml.safe_dump(sets, sort_keys=False))
    (run_dir / FROZEN_CONFIG).write_text(
        f"# Copied from {Path(config_path).resolve()} by ./noodle calibrate setup "
        f"on {datetime.now():%Y-%m-%d %H:%M}; ./noodle calibrate report reads this copy.\n"
        + Path(config_path).read_text())
    run = dict(description=config["description"], sequences_csv=SEQUENCE_FILE, name_column="id",
               sequence_column="sequence", rows="all", max_residues=0, physics=list(sets),
               replicates=config["replicates"], time_limit_days=config["time_limit_days"],
               max_running=config["max_running"], extra_sbatch=config["extra_sbatch"])
    (run_dir / "run.yaml").write_text(
        "# Written by ./noodle calibrate setup: edit experiments/brush_calibration.yaml and run setup\n"
        "# again instead (until it is launched). time_limit_days, max_running and extra_sbatch\n"
        "# can still be changed here after the launch.\n" + yaml.safe_dump(run, sort_keys=False))

    all_sets = physics_sets(run_dir)
    problems = []
    for name in sets:
        physics = resolve_physics(all_sets, name)
        problems += check_physics_set(name, physics, all_sets)
        try:
            run_lengths(physics, int(sequences["N"].min()))
        except ValueError as error:
            problems.append(f"{name}: {error}")
    if problems:
        sys.exit("The calibration physics sets have problems:\n  " + "\n  ".join(problems))

    print(f"Set up runs/{config['run_name']}/ from {config_path}\n")
    print("Sequences (one per length and kind):")
    for _, s in sequences.iterrows():
        print(f"  {s['id']:34s} N={s['N']:<4} hydropathy {s['hydropathy']:.3f}  fcr {s['fcr']:.2f}  "
              f"ncpr {s['ncpr']:+.2f}")
    print("\nConditions (physics sets, in runs/{0}/physics_sets.yaml):".format(config["run_name"]))
    for condition, name in zip(config["conditions"], sets):
        physics = resolve_physics(all_sets, name)
        chains, width = brush_size(physics)
        print(f"  {name:22s} {chains} chains, box {width:.1f} nm wide, run {condition['run_multiplier']} x N^2.2 "
              f"time units, {config['frames']} frames from step 0")
    print(f"\n{len(sequences)} sequences x {len(sets)} conditions x {config['replicates']} replicate(s) = "
          f"{len(sequences) * len(sets) * config['replicates']} tasks.")
    print("\nNext, on the cluster:")
    print(f"  ./noodle launch {config['run_name']}      shows the cost and disk, asks, submits")
    print(f"  ./noodle status {config['run_name']}      how it's going; launch again resumes tasks that hit the time limit")
    print("  ./noodle calibrate report        any time: analyses finished tasks, deletes their trajectories,")
    print(f"                                   writes {config['results_dir']}/report.md")


# ── report: one task ──────────────────────────────────────────────────────────
def equilibration(times, y):
    """When a chain-averaged observable has relaxed from the start to its plateau.

    plateau: the mean of the second half of the run. The series is smoothed over 1% of
    the run; the tolerance is 4 standard deviations of the smoothed second half (at least
    0.5% of the plateau). t_eq is the time right after the last smoothing window that is
    outside plateau +- tolerance. ok: t_eq lies in the first half of the run (else the run
    was too short to tell). drift: the third and fourth quarters differ by more than the
    tolerance (the "plateau" still moves)."""
    n = len(y)
    w = max(5, n // 100)
    smooth = np.convolve(y, np.ones(w) / w, mode="valid")  # smooth[i] = mean of y[i:i+w]
    plateau = y[n // 2:].mean()
    tolerance = max(4 * smooth[len(smooth) // 2:].std(), 0.005 * abs(plateau))
    outside = np.flatnonzero(np.abs(smooth - plateau) > tolerance)
    i = 0 if len(outside) == 0 else min(outside[-1] + w, n - 1)
    drift = abs(y[3 * n // 4:].mean() - y[n // 2:3 * n // 4].mean()) > tolerance
    return dict(t_eq=float(times[i]), plateau=float(plateau), start=float(y[0]),
                tolerance=float(tolerance), ok=bool(i < n // 2), drift=bool(drift))


def autocorrelation(x):
    """Normalised autocorrelation function of x (frames first, then chains and, for a
    vector, its components), averaged over the chains; vectors use the dot product."""
    x = x - x.mean(axis=0)
    n = len(x)
    spectrum = np.fft.rfft(x, 2 * n, axis=0)
    covariance = np.fft.irfft(spectrum * np.conj(spectrum), axis=0)[:n]
    covariance /= (n - np.arange(n)).reshape((n,) + (1,) * (x.ndim - 1))  # unbiased at every lag
    covariance = covariance.reshape(n, -1).mean(axis=1)  # over chains (and components)
    return covariance / covariance[0]


def relaxation_times(acf, frame_time):
    """(integrated time, 1/e time) of an autocorrelation function, in time units.
    integrated: frame_time x (1/2 + the sum of the ACF up to its first negative value),
    as ./noodle analyze uses for independent_samples. 1/e: where it first drops below
    1/e (interpolated); NaN if it never does within the data."""
    first_negative = int(np.argmax(acf < 0)) or len(acf)
    integrated = frame_time * (0.5 + acf[1:first_negative].sum())
    below = np.flatnonzero(acf < 1 / math.e)
    if len(below) == 0:
        return integrated, np.nan
    k = below[0]
    fraction = (acf[k - 1] - 1 / math.e) / (acf[k - 1] - acf[k])
    return integrated, frame_time * (k - 1 + fraction)


def analyse_task(run_dir, task, log, physics, sequence_info, results):
    """Everything the report needs from one finished task, into results/tasks/:
    task_XXXX.timeseries.csv.gz, task_XXXX.acf.csv.gz, and task_XXXX.json (written
    last: a task counts as analysed once its json exists)."""
    task_id, n = int(task["task_id"]), int(task["length"])
    chains, width = brush_size(physics)
    trajectory = run_dir / log["trajectory"]
    table, _ = trajectory_metrics(trajectory, chains, log["geometry"]["z_wall"])
    with gsd.fl.open(trajectory, "r") as f:
        steps = np.array([int(f.read_chunk(i, "configuration/step")[0]) for i in range(f.nframes)])
    frames = len(steps)
    times = steps * physics["dt"]
    per_chain = {name: table[name].to_numpy().reshape(frames, chains)
                 for name in ("rg", "ree_x", "ree_y", "ree_z", "com_height", "asphericity")}
    ree = np.stack([per_chain["ree_x"], per_chain["ree_y"], per_chain["ree_z"]], axis=-1)

    averaged = pd.DataFrame(dict(
        time=times, step=steps, rg=per_chain["rg"].mean(1), rg_sq=(per_chain["rg"] ** 2).mean(1),
        ree_z=per_chain["ree_z"].mean(1), ree_sq=(ree ** 2).sum(-1).mean(1),
        com_height=per_chain["com_height"].mean(1), asphericity=per_chain["asphericity"].mean(1),
        com_height_sd=per_chain["com_height"].std(1), ree_z_sd=per_chain["ree_z"].std(1)))

    eq = {name: equilibration(times, averaged[name].to_numpy()) for name in EQUILIBRATION_OBSERVABLES}
    equilibrated = all(e["ok"] for e in eq.values())
    t_eq = max(e["t_eq"] for e in eq.values())
    # relaxation times from the frames after equilibration (the second half if it was not reached)
    first = int(np.searchsorted(times, t_eq)) if equilibrated else frames // 2
    frame_time = float(np.median(np.diff(times)))
    series = dict(rg=per_chain["rg"][first:], ree_z=per_chain["ree_z"][first:],
                  com_height=per_chain["com_height"][first:], ree_vector=ree[first:],
                  collective_height=averaged["com_height"].to_numpy()[first:, None])
    acfs = {name: autocorrelation(x) for name, x in series.items()}
    taus = {name: relaxation_times(acf, frame_time) for name, acf in acfs.items()}
    used = frames - first
    lags = used // 2
    pd.DataFrame(dict(lag_time=np.arange(lags) * frame_time, **{k: v[:lags] for k, v in acfs.items()})).to_csv(
        results / "tasks" / f"task_{task_id:04d}.acf.csv.gz", index=False, float_format="%.5g")
    averaged.to_csv(results / "tasks" / f"task_{task_id:04d}.timeseries.csv.gz", index=False, float_format="%.6g")

    tau_chain = max(taus[name][0] for name in CHAIN_OBSERVABLES)
    tau_fastest = min(taus[name][0] for name in CHAIN_OBSERVABLES)
    run_time_s = log.get("run_time_s") or np.nan
    steps_per_second = log["total_steps"] / run_time_s if run_time_s else np.nan
    scale = n ** EXPONENT
    row = dict(
        task_id=task_id, sequence_name=task["sequence_name"], kind=sequence_info.get("kind", ""),
        N=n, physics=task["physics"], grafting_density=physics["grafting_density"],
        wall_charge_density=physics["wall_charge_density"], chains=chains, box_width_nm=round(width, 3),
        hydropathy=sequence_info.get("hydropathy"), fcr=sequence_info.get("fcr"), ncpr=sequence_info.get("ncpr"),
        seed=int(task["seed"]), dt=physics["dt"], frames=frames, frame_time=frame_time,
        run_time_units=float(times[-1]), run_multiplier=float(times[-1] / scale),
        t_eq=t_eq, c_eq=t_eq / scale, equilibrated=equilibrated,
        drift=any(e["drift"] for e in eq.values()), relaxation_from_time=float(times[first]),
        **{f"t_eq_{k}": e["t_eq"] for k, e in eq.items()},
        **{f"start_{k}": e["start"] for k, e in eq.items()},
        **{f"plateau_{k}": e["plateau"] for k, e in eq.items()},
        **{f"tau_{k}": t[0] for k, t in taus.items()},
        **{f"tau_1e_{k}": t[1] for k, t in taus.items()},
        tau_chain=tau_chain, c_tau=tau_chain / scale, c_tau_collective=taus["collective_height"][0] / scale,
        relaxation_times_covered=(times[-1] - times[first]) / tau_chain,
        frames_per_fastest_tau=tau_fastest / frame_time,
        independent_samples_rg=used / (2 * taus["rg"][0] / frame_time),
        total_steps=log["total_steps"], run_time_s=run_time_s, steps_per_second=steps_per_second,
        residue_steps_per_second=steps_per_second * chains * n,
        last_hourly_steps_per_second=log.get("steps_per_second"), host=log.get("host"),
        attempts=log.get("attempts"), trajectory_bytes=trajectory.stat().st_size,
        hoomd_version=log.get("hoomd_version"), completed_at=log.get("completed_at"))
    row = {k: (None if isinstance(v, float) and not math.isfinite(v) else
               v.item() if isinstance(v, np.generic) else v) for k, v in row.items()}
    path = results / "tasks" / f"task_{task_id:04d}.json"
    path.with_suffix(".tmp").write_text(json.dumps(row, indent=1))
    path.with_suffix(".tmp").replace(path)


def task_files(run_dir, task_id):
    """The task's trajectory, checkpoint and half-written files."""
    stem = f"task_{task_id:04d}"
    return [p for folder in ("trajectories", "checkpoints") for p in (run_dir / folder).glob(f"{stem}[_.]*")]


# ── report: all tasks ─────────────────────────────────────────────────────────
def ceil_half(x):
    return max(0.5, math.ceil(2 * x) / 2) if math.isfinite(x) else float("nan")


def md_table(df, digits=3):
    def cell(v):
        if isinstance(v, (float, np.floating)):
            return "" if not math.isfinite(v) else f"{v:.{digits}g}"
        return str(v)
    lines = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
    lines += ["| " + " | ".join(cell(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines)


def fit_exponent(group, column):
    """(a, b) of column = a N^b, least squares in log-log; NaN with fewer than 2 lengths."""
    data = group[[column, "N"]].dropna()
    data = data[data[column] > 0]
    if data["N"].nunique() < 2:
        return np.nan, np.nan
    b, log_a = np.polyfit(np.log(data["N"]), np.log(data[column]), 1)
    return math.exp(log_a), b


def condition_summary(tasks, config):
    """One row per condition: the slowest equilibration and relaxation, exponents, speed,
    and the multipliers they call for."""
    safety, samples = config["equilibration_safety"], config["target_independent_samples"]
    rows = []
    for physics, group in tasks.groupby("physics", sort=False):
        _, b_tau = fit_exponent(group, "tau_chain")
        _, b_eq = fit_exponent(group, "t_eq")
        rows.append(dict(
            physics=physics, density=group["grafting_density"].iloc[0], charge=group["wall_charge_density"].iloc[0],
            tasks=len(group), lengths=",".join(map(str, sorted(group["N"].unique()))),
            c_eq_max=group["c_eq"].max(), c_eq_median=group["c_eq"].median(),
            c_tau_max=group["c_tau"].max(), c_tau_median=group["c_tau"].median(),
            c_tau_collective_max=group["c_tau_collective"].max(), exponent_tau=b_tau, exponent_eq=b_eq,
            not_equilibrated=int((~group["equilibrated"]).sum()), min_times_covered=group["relaxation_times_covered"].min(),
            residue_steps_per_s=group["residue_steps_per_second"].median(),
            rec_rouse_multiplier=ceil_half(safety * group["c_eq"].max()),
            rec_production_multiplier=ceil_half(2 * samples * group["c_tau"].max())))
    return pd.DataFrame(rows)


def cost_table(summary, tasks, config, lengths=(25, 50, 100, 150, 200, 300)):
    """Core-hours of one task at the recommended multipliers, with the speed (residue-steps
    per second) measured for that condition at the nearest measured length."""
    rows = []
    for _, c in summary.iterrows():
        multiplier = c["rec_rouse_multiplier"] + c["rec_production_multiplier"]
        measured = (tasks[tasks["physics"] == c["physics"]].groupby("N")["residue_steps_per_second"]
                    .median().dropna())
        row = {"physics": c["physics"], "multipliers (eq + prod)": f"{c['rec_rouse_multiplier']:g} + "
                                                                  f"{c['rec_production_multiplier']:g}"}
        for n in lengths:
            if measured.empty or not math.isfinite(multiplier):
                row[f"N={n}"] = ""
                continue
            speed = measured.iloc[int(np.argmin(np.abs(measured.index.to_numpy() - n)))]
            steps = multiplier * n ** EXPONENT / 0.02  # dt of every calibration set (checked in report)
            row[f"N={n}"] = f"{config['chains'] * n * steps / speed / 3600:,.1f} h"
        rows.append(row)
    return pd.DataFrame(rows)


def write_report(run_dir, config, physics, tasks, sequences, status, errors, results, gsd_left):
    all_tasks = read_tasks(run_dir)
    launched = datetime.fromtimestamp((run_dir / "tasks.csv").stat().st_mtime)
    current_brush = resolve_physics(yaml.safe_load((ROOT / "physics.yaml").read_text()), config["base_set"])
    lines = [f"# Brush calibration report: `{config['run_name']}`", "",
             f"Written by `./noodle calibrate report` on {datetime.now():%Y-%m-%d %H:%M} on `{socket.gethostname()}`, "
             f"noodle-md {git_version()}. Run launched {launched:%Y-%m-%d %H:%M}.", "",
             f"**{len(tasks)} of {len(all_tasks)} tasks analysed.** "
             + ("All tasks are in: this is the final report." if len(tasks) == len(all_tasks) else
                "Not all tasks are in yet: rerun `./noodle calibrate report` once more have finished."), ""]
    if status:
        lines += ["Not analysed yet: " + ", ".join(f"{n} {s}" for s, n in status.items()) + ".", ""]
    if errors:
        lines += ["Tasks that could not be analysed (their files are kept):", ""]
        lines += [f"- task {i}: {message}" for i, message in errors] + [""]
    if gsd_left:
        lines += [f"Trajectory/checkpoint files still in `runs/{config['run_name']}/`: {gsd_left[0]} "
                  f"({gsd_left[1]}), of tasks not analysed yet.", ""]

    lines += ["## 1. Context", "",
              "**Why.** Brushes in noodle-md are 10 x 10 patches of grafted chains (`chains: 100` in "
              "`physics.yaml`), simulated on one CPU core per task with HOOMD-blue (HPS model, Langevin "
              "dynamics). Run lengths are `multiplier x N^2.2 / dt` steps (`rouse_multiplier` for "
              "equilibration, `production_multiplier` for production). Before this experiment the brush "
              "multipliers were provisional guesses (from earlier 4 x 4 measurements). This experiment "
              "measures, for 10 x 10 brushes, (1) how long a brush takes to equilibrate from the starting "
              "configuration, (2) the chains' relaxation times in equilibrium, (3) the simulation speed on "
              "the cluster, to set the multipliers, `frames` and the cost model.", "",
              "**Design.** Each task is one brush of one sequence under one condition (grafting density, "
              "surface charge), run for `run_multiplier x N^2.2` time units and recorded from step 0 with "
              f"{config['frames']} frames. Starting configuration (`simulate.py`): every chain grows upward "
              "from its anchor as a random walk with non-negative z steps inside its own lattice column, so "
              "the brush starts more stretched than in equilibrium. Anchors (N-terminus) are fixed at the "
              "contact distance of an lj93 repulsive surface.", "",
              "**Units.** Time in HOOMD units of 0.489 ps (dt = 0.02 per step). `c_*` columns are times "
              "divided by N^2.2, i.e. directly comparable with the multipliers in `physics.yaml`.", "",
              f"**The `{config['base_set']}` set in physics.yaml when this report was written** (the "
              "calibration sets extend it):", "", "```yaml",
              yaml.safe_dump({config["base_set"]: current_brush}, sort_keys=False).rstrip(), "```", "",
              "**Config of this experiment** (frozen copy in the run folder):", "", "```yaml",
              (run_dir / FROZEN_CONFIG).read_text().rstrip(), "```", "",
              "**Sequences**", "",
              md_table(sequences[["id", "kind", "N", "hydropathy", "fcr", "ncpr", "source_id"]]), "",
              "**Methods** (`src/calibrate.py`).", "",
              "- Per frame and chain: Rg, end-to-end vector (`ree_z` its height component), centre-of-mass "
              "height above the surface (`com_height`), as `./noodle analyze` computes them.",
              "- *Equilibration time* `t_eq`: for the chain-averaged com_height, ree_z and Rg, smoothed over "
              "1% of the run, the time after the last excursion outside plateau +- tolerance (plateau = "
              "mean of the 2nd half; tolerance = 4 SD of the smoothed 2nd half, at least 0.5% of the plateau). "
              "`t_eq` = the largest of the three. `equilibrated` = false when that is in the second half "
              "(run too short to tell); `drift` = the 3rd and 4th quarters differ by more than the tolerance.",
              "- *Relaxation times* from the frames after `t_eq` (the 2nd half if not equilibrated): per-chain "
              "autocorrelation functions (ACF) of Rg, ree_z, com_height and the end-to-end vector "
              "(`ree_vector`, dot product), averaged over the 100 chains. `tau_X` = integrated time "
              "(frame time x (1/2 + sum of the ACF up to its first negative value)), the same estimator "
              "`./noodle analyze` uses for `independent_samples`; `tau_1e_X` = where the ACF drops below "
              "1/e. `tau_chain` = the largest per-chain `tau_X`. `collective_height` = ACF of the "
              "chain-averaged com_height (the brush height as a whole; noisier, one series per task).",
              "- *Speed*: total steps / run time, from the task log (includes start-up; cumulative over "
              "resumed attempts). `residue_steps_per_second` = that x chains x N.",
              "- *Recommendations*: `rouse_multiplier` = equilibration_safety x max `c_eq`; "
              "`production_multiplier` = 2 x target_independent_samples x max `c_tau` (per chain, "
              "independent samples = production / (2 tau)); both rounded up to 0.5. A frame spacing of "
              "about tau/5 is enough, i.e. `frames` ~ 10 x target_independent_samples.", ""]

    summary = condition_summary(tasks, config)
    lines += ["## 2. Results per condition", "",
              md_table(summary.drop(columns=["rec_rouse_multiplier", "rec_production_multiplier"])), "",
              "`exponent_tau`/`exponent_eq`: fitted b in tau = a N^b over the lengths measured (2.2 is what "
              "the multipliers assume: a larger b means the c values grow with N, so the recommendation for "
              "N beyond the measured lengths should be scaled up by (N/N_max)^(b-2.2)). "
              "`min_times_covered`: the fewest relaxation times any task of the condition covered after "
              "equilibration (below ~30 the taus are rough).", ""]

    lines += ["## 3. Recommendations (derived automatically: check them against section 2 and 4)", "",
              md_table(summary[["physics", "c_eq_max", "c_tau_max", "exponent_tau", "rec_rouse_multiplier",
                                "rec_production_multiplier"]]), "",
              f"Suggested frames: {10 * config['target_independent_samples']} "
              f"(10 x target_independent_samples = {config['target_independent_samples']}).", ""]
    dilute = summary[summary["density"] <= 0.1]
    if len(dilute):
        lines += [f"One setting for every condition up to 0.1 chains/nm^2 (the `brush` default is "
                  f"{current_brush.get('grafting_density')}): `rouse_multiplier: "
                  f"{dilute['rec_rouse_multiplier'].max():g}`, `production_multiplier: "
                  f"{dilute['rec_production_multiplier'].max():g}`.", ""]
    lines += ["**Cost of one task at the recommended multipliers** (one core; residue-steps per second "
              "measured for the condition at the nearest measured length; beyond the measured lengths it "
              "assumes tau ~ N^2.2, see exponent_tau):", "", md_table(cost_table(summary, tasks, config)), "",
              f"Median speed over all tasks: {tasks['residue_steps_per_second'].median():,.0f} residue-steps "
              f"per second per core (launch_run.py assumes {1e6:,.0f} for its rough cost).", ""]

    by_n = tasks.groupby(["physics", "N"]).agg(
        tasks=("task_id", "count"), residue_steps_per_s=("residue_steps_per_second", "median"),
        run_hours=("run_time_s", lambda s: s.median() / 3600), trajectory_MB=("trajectory_bytes", lambda s: s.median() / 1e6),
        c_eq_max=("c_eq", "max"), c_tau_max=("c_tau", "max")).reset_index()
    lines += ["## 4. Per condition and length", "", md_table(by_n), "",
              "## 5. Every task", "",
              md_table(tasks[["task_id", "sequence_name", "physics", "N", "equilibrated", "drift", "t_eq", "c_eq",
                              "tau_rg", "tau_ree_z", "tau_com_height", "tau_ree_vector", "tau_collective_height",
                              "c_tau", "relaxation_times_covered", "frames_per_fastest_tau",
                              "residue_steps_per_second", "host"]]), ""]

    lines += ["## 6. Files, and how to pick this up again", "",
              f"Everything is in `{config['results_dir']}/` (small enough to commit to git):", "",
              "- `report.md`: this file. `tasks.csv`: one row per analysed task, every number above and more "
              "(start and plateau values, 1/e times, hosts, attempts, trajectory sizes). `conditions.csv`: "
              "section 2 and 3 as a table.",
              "- `tasks/task_XXXX.timeseries.csv.gz`: chain-averaged observables of every saved frame from "
              "step 0 (`time` in time units, `step`, `rg`, `rg_sq`, `ree_z`, `ree_sq`, `com_height`, "
              "`asphericity`, and the spread over chains `com_height_sd`, `ree_z_sd`). Enough to redo the "
              "equilibration analysis with another criterion.",
              "- `tasks/task_XXXX.acf.csv.gz`: the ACFs (`lag_time`, `rg`, `ree_z`, `com_height`, "
              "`ree_vector`, `collective_height`) from the equilibrated part, up to half its length. Enough "
              "to refit relaxation times (e.g. exponential tails).",
              "- `tasks/task_XXXX.json`: the task's row of `tasks.csv`.",
              f"- In the run folder `runs/{config['run_name']}/` (cluster): `tasks.csv`, `physics.yaml` (the "
              "sets as run, fully written out), `logs/` (per-task speed, host), `slurm/`. The trajectories "
              "are deleted once analysed: this experiment keeps no .gsd files.", "",
              "**Next steps this report is meant for:**", "",
              "1. Set `rouse_multiplier`, `production_multiplier` and `frames` of the `brush` set in "
              "`physics.yaml` from section 3 (dense brushes may need their own set), and update the README "
              "(section 3 cost table, section 4 run-length note).",
              "2. Update the brush cost model: `RESIDUE_STEPS_PER_SECOND` in `src/launch_run.py` (one value "
              "for every geometry now) from section 2/4, possibly per density.",
              "3. Check `exponent_tau`: if far from 2.2, the N^2.2 rule under- or over-runs long chains; "
              "consider a different exponent for brushes.",
              "4. If `c_tau_collective` is much larger than `c_tau`, the brush height as a whole relaxes "
              "slower than single chains: error bars on brush averages then need several replicates.",
              "5. Conditions with `not_equilibrated` > 0 or `min_times_covered` < 30: rerun them longer "
              "(a larger run_multiplier in the config, a new run_name) before trusting their numbers.", ""]
    (results / "report.md").write_text("\n".join(lines))
    summary.to_csv(results / "conditions.csv", index=False)
    return summary


def report(config_path):
    config = load_config(config_path)
    run_dir = RUNS / config["run_name"]
    if not (run_dir / "tasks.csv").exists():
        sys.exit(f"runs/{config['run_name']} has not been launched yet: ./noodle calibrate setup, then "
                 f"./noodle launch {config['run_name']}")
    config = yaml.safe_load((run_dir / FROZEN_CONFIG).read_text())  # the config as it was run
    results = ROOT / config["results_dir"]
    (results / "tasks").mkdir(parents=True, exist_ok=True)
    launched = datetime.fromtimestamp((run_dir / "tasks.csv").stat().st_mtime).isoformat(timespec="seconds")
    marker = dict(run_name=config["run_name"], launched=launched)
    if (results / MARKER).exists() and json.loads((results / MARKER).read_text()) != marker:
        sys.exit(f"{config['results_dir']} holds the results of another run ({(results / MARKER).read_text()}). "
                 "Give this experiment its own results_dir.")
    (results / MARKER).write_text(json.dumps(marker))

    physics = yaml.safe_load((run_dir / "physics.yaml").read_text())
    if any(p["dt"] != 0.02 for p in physics.values()):
        print("Note: the cost table assumes dt = 0.02; this run used another dt.")
    sequences = pd.read_csv(run_dir / SEQUENCE_FILE)
    info = {row["id"]: row for row in sequences.to_dict("records")}
    status, errors, analysed, freed = {}, [], 0, 0
    tasks = read_tasks(run_dir)
    for task in tasks:
        task_id = int(task["task_id"])
        log = read_log(run_dir, task_id)
        if log.get("status") != "completed":
            state = log.get("status") or "not started"
            state = "unfinished" if state == "running" else state
            status[state] = status.get(state, 0) + 1
            continue
        if not (results / "tasks" / f"task_{task_id:04d}.json").exists():
            if not (run_dir / log["trajectory"]).exists():
                errors.append((task_id, "completed, but its trajectory is not here"))
                continue
            print(f"  analysing task {task_id} ({task['sequence_name']}, {task['physics']})", flush=True)
            try:
                analyse_task(run_dir, task, log, physics[task["physics"]], info.get(task["sequence_name"], {}), results)
            except Exception as error:
                errors.append((task_id, f"{type(error).__name__}: {error}"))
                continue
            analysed += 1
        for path in task_files(run_dir, task_id):  # analysed: its trajectory is no longer needed
            freed += path.stat().st_size
            path.unlink()

    rows = [json.loads(p.read_text()) for p in sorted((results / "tasks").glob("task_*.json"))]
    left = [p for folder in ("trajectories", "checkpoints") for p in (run_dir / folder).glob("*") if p.is_file()]
    gsd_left = (len(left), file_size(sum(p.stat().st_size for p in left))) if left else None
    print(f"{len(rows)} of {len(tasks)} tasks analysed ({analysed} new); deleted {file_size(freed)} of trajectories.")
    if errors:
        print(f"  {len(errors)} could not be analysed: " + "; ".join(f"task {i}: {m}" for i, m in errors[:5]))
    if not rows:
        print("No finished task yet: nothing to report.")
        return
    table = pd.DataFrame(rows).sort_values("task_id")
    table.to_csv(results / "tasks.csv", index=False)
    summary = write_report(run_dir, config, physics, table, sequences, status, errors, results, gsd_left)
    print(f"  -> {config['results_dir']}/report.md (and tasks.csv, conditions.csv, tasks/)")
    print(summary[["physics", "tasks", "c_eq_max", "c_tau_max", "rec_rouse_multiplier",
                   "rec_production_multiplier"]].to_string(index=False))
    if gsd_left:
        print(f"  {gsd_left[0]} trajectory/checkpoint file(s) of unfinished tasks remain ({gsd_left[1]}); they go "
              "once those tasks finish and are reported, or with ./noodle calibrate cleanup")
    else:
        print("  no trajectory or checkpoint files remain in the run folder")
    if len(rows) == len(tasks):
        print(f"All tasks are in. Share the results with:  git add {config['results_dir']} && git commit "
              f"-m 'Brush calibration results' && git push   (or copy the folder)")


def cleanup(config_path):
    config = load_config(config_path)
    report(config_path)
    run_dir = RUNS / config["run_name"]
    left = [p for folder in ("trajectories", "checkpoints") for p in (run_dir / folder).glob("*") if p.is_file()]
    if not left:
        print("\nNothing to clean up.")
        return
    print(f"\nDelete {len(left)} file(s) ({file_size(sum(p.stat().st_size for p in left))}) of tasks that are "
          "not finished (or could not be analysed): they can then no longer be resumed.")
    try:
        answer = input(f"Type the run name ({config['run_name']}) to delete them, anything else to cancel: ")
    except EOFError:
        answer = ""
    if answer.strip() != config["run_name"]:
        print("Cancelled: nothing was deleted.")
        return
    for path in left:
        path.unlink()
    print(f"Deleted {len(left)} file(s).")


def main():
    stages = dict(setup=setup, report=report, cleanup=cleanup)
    if len(sys.argv) not in (2, 3) or sys.argv[1] not in stages:
        sys.exit("Usage: ./noodle calibrate setup|report|cleanup [config file]\n"
                 f"       (default config: {DEFAULT_CONFIG.relative_to(ROOT)}; see the top of that file)")
    stages[sys.argv[1]](Path(sys.argv[2]) if len(sys.argv) == 3 else DEFAULT_CONFIG)


if __name__ == "__main__":
    main()
