# Changelog

What changed between versions of noodle-md, for anyone who pulls an update.
Runs launched with one version can't be continued or analysed by another version
with a different run layout: finish them with the version that launched them
(`git checkout v1`, for example).

## v2 (in progress, branch `release`)

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
