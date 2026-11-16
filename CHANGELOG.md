# Changelog

What changed between versions of noodle-md, for anyone who pulls an update.
Runs launched with one version can't be continued or analysed by another version
with a different run layout: finish them with the version that launched them
(`git checkout v1`, for example).

## v2 (in progress, branch `release`)

### Added
- `runs/<name>/manifest.csv`, written at the first launch: one row per protein and
  physics set with its name, sequence, where it came from (`source_csv`, `source_row`),
  the physics conditions (`geometry`, `temperature_K`, `ionic_concentration_M`,
  `graft_terminus`, `wall_charge_density`, `wall_attraction`, `grafting_density`,
  `chains`) and the sequence descriptors (`src/descriptors.py`): `N`, `net_charge`,
  `ncpr`, `fcr`, `f_plus`, `f_minus`, `mean_hydro`, `kappa`, `scd`, `scd_pos`,
  `scd_neg`, `shd`. HPS1 charges (histidine +0.5) and idp-rg's formulas and names, so
  the values agree with idp-rg (checked: to float32 precision on both shared sets).
  `kappa` uses Das & Pappu's block-arrangement delta_max (the E/K test sequences give
  0.0009 and 1, as published); it can exceed 1 for ~2% of sequences, mostly ones with
  few charges.

- **Every task is analysed by its own job, right after it finishes** (`src/observables.py`;
  seconds for a chain, about half a minute for a brush). Averages, each with a
  block-average error bar (`<name>_err`), go into the task's row of its job file; profiles
  go into `jobs/job_0007_profiles.csv`. Lengths in nm, times in model units (tu).
  - every chain: `rg_nm`, `rg_sq_nm2`, `rg_sd_nm`, `rh_nm` (Kirkwood), `ree_nm`,
    `ree_sq_nm2`, `asphericity`, `shape_anisotropy` (relative shape anisotropy kappa^2,
    not the sequence kappa), `prolateness`, `transverse_rms_nm`, `nu` and `r0_nm`
    (R_ij = r0 |i-j|^nu for |i-j| >= 5), `bond_cos`, `persistence_nm` (Flory's
    projection), `contacts_per_residue` (|i-j| >= 3, closer than 1.5 x the mean diameter)
  - grafted chains and brushes: `com_height_nm`, `end_height_nm`, `max_height_nm`,
    `height_95_nm`, `rg_perp_nm`, `rg_par_nm`, `end_lateral_rms_nm`, `tilt_cos`,
    `bond_order` (P2 of the bonds with the normal), `adsorbed_fraction`
  - brushes: `brush_height_nm` (2 x mean residue height), `concentration_mg_ml` (below
    height_95), `interchain_contacts_per_residue`, `lateral_heterogeneity`; and, in the
    manifest, `reduced_density` (grafting density x pi <Rg^2> of the same protein as a
    single chain in the same run)
  - diagnostics: `tau_rg_tu`, `tau_ree_tu`, `independent_samples`, `drift`
  - profiles: `internal_distance`, `bond_correlation`, `residue_contacts`; grafted:
    `density`, `mass_density`, `charge_density`, `end_density` (vs height),
    `residue_height`, `residue_adsorption`; brushes: `residue_interchain_contacts`
  Checked on shapes with known answers (rod, ring, random walk, upright grafted rod,
  synthetic brush), interchain contacts against brute force, and Rg, kappa^2 and Ree
  against MDAnalysis on simulated trajectories (identical).
- **Run numbers**: `runs/runs.csv` lists every run with a number (given by `./noodle new`,
  never reused) and run-wide figures: description, created / launched / finished dates,
  proteins, physics sets, replicates, tasks, jobs, done, failed, analysed, core-hours,
  trajectory size (GB), last checked. `new`, `launch`, `status` and `analyze` keep it up
  to date (never the jobs). Runs made before it get a number the first time it is read.
  `./noodle status` without a run prints it.
- **Choosing a run without typing its name**: wherever a run is asked for, give its name,
  its number, a path to its folder (what Tab completes to, e.g. `runs/my_run/`), or
  nothing when you are inside its folder (or a subfolder).
- `./noodle cancel <run>`: counts the run's jobs in the queue, asks, then `scancel`s them.
  Running tasks save a checkpoint first; `./noodle launch` continues the run later.
- `./noodle help` lists the commands; `./noodle help <command>` explains one in detail
  (the description at the top of its program). Unknown commands show the list too.
- Tab completion (`noodle-completion.bash`): commands, then run names. `setup.sh` ends by
  printing two lines to add to `~/.bashrc`: an alias, so `noodle` works from any folder,
  and the line that loads the completion.
- `./noodle export <run> [<run> ...]`: brings each run's analysis up to date, then writes
  `exports/<name>/` with the runs' `manifest.csv`, `replicates.csv` and `profiles.csv`
  merged (first column `run`), their rows of `runs.csv`, and `physics/<run>.yaml`. One
  run: named after it; several: it asks, suggesting the common part of their names. Warns
  about unfinished runs and about a physics set name with different settings in different
  runs. `exports/` is ignored by git.
- scipy is listed in `src/environment.yml` (it was already installed through MDAnalysis).

### Changed
- Run names can't be only digits (those are run numbers now).
- README rewritten for all of the above: the `noodle` alias and Tab completion, run
  numbers, jobs, `log.csv`, cancel, export, and every column of the results tables with
  its definition and caveats (model times, Kirkwood Rh, provisional brush run length,
  kappa > 1). `src/run_template.yaml` explains the per-job time limit.
- `./noodle check` runs its test simulations as one job and checks the new tables
  (in-job analysis, manifest columns per geometry, brush density profiles, reduced density).
- **`./noodle analyze` pools instead of recomputing**: it gathers the job files into
  `analysis/replicates.csv` (one row per task), `analysis/profiles.csv` (each profile
  averaged over replicates, with `err`) and `analysis/report.json` (counts, the range of
  each observable per physics set, rows that drift or have fewer than 20 independent
  samples, tasks that could not be analysed), and fills the results columns of
  `manifest.csv`: `replicates`, then every observable averaged over replicates, with the
  standard error between replicates (or a lone replicate's block error). It takes
  seconds. A task its job could not analyse is analysed by `analyze` instead
  (`analysis/reanalysed.csv`). `analysis/summary.csv`, `analysis/timeseries/` and
  `analysis/profiles/` are gone: the per-task results in `jobs/` replace them, and keep
  working after `./noodle clean` deletes the trajectories.
- **Tasks are grouped into jobs.** A job is one SLURM array element (one core) that runs
  its tasks one after another. `./noodle launch` packs short tasks into jobs of up to
  ~4 h (estimated; at most half of `time_limit_days`), longest first ("first fit
  decreasing"); a longer task gets a job of its own. The longest jobs get the lowest
  numbers, so they start first. A failed task doesn't stop its job. Thousands of
  minute-long tasks now make a few hundred jobs instead of thousands, and HOOMD is loaded
  once per job. `max_running` now counts jobs (one core each). Launching looks the same;
  the plan's Cluster line says how many jobs.
- **New bookkeeping layout** (fewer files): `logs/` (one JSON per task) and `slurm/` (one
  output file per task) are replaced by `jobs/`:
  - `jobs/job_0007.csv`: the tasks of job 7 and how each is going (status, phase, step,
    speed, run time, attempts, host, dates, error), written by launch just before the job
    is submitted, then only by the job itself. A task's current state is its row in the
    newest job file that lists it (a resumed or retried task gets a row in a new job).
  - `jobs/job_0007.out`: everything job 7 printed (all of its tasks).
  - `jobs/submitted.csv`: which SLURM job runs which job.
  - `jobs/slurm.out`: SLURM's own messages (normally empty).
- `./noodle status` writes `log.csv` (one row per task: state, progress, job, SLURM id,
  speed, run time, time left, attempts, host, dates, error) instead of `status.txt`.
  Tasks waiting in a job that is still running other tasks count as waiting.
- A task that is about to hit the cluster time limit stops cleanly: SLURM warns it 10
  minutes ahead (`#SBATCH --signal=B:USR1@600` in `src/slurm_job.sh`), and the task saves
  a checkpoint and records `stopped (time limit)`. Nothing is lost (before: up to an
  hour). `scancel` and Ctrl-C (own computer) stop it the same way; a second Ctrl-C quits
  at once.
- The simulation runs in chunks of about a minute (sized from the measured speed)
  instead of 10,000 steps, so a stop request is always answered within about a minute.

### Fixed
- Replicates of the same protein under the same physics set always get different
  random seeds (before, two could draw the same seed and be the same simulation).

### Removed
- `./noodle calibrate` and `experiments/` (the brush calibration experiment). It was a
  private test harness and is not part of the released tool. It lives on in the tag
  `v1` and the local branch `calibration`.
- A run's own `physics_sets.yaml` is still supported (sets that only one run needs).

## v1 (tag `v1`)

The first version used on the cluster: bulk, grafted and charged-wall single chains,
10 x 10 brushes, `./noodle new | launch | status | analyze | clean | check | calibrate`.
