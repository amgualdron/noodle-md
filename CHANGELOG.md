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

### Changed
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
