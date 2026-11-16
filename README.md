# noodle-md

Simulations of intrinsically disordered proteins (IDPs). Each amino acid is one bead,
and the beads form a chain (the HPS model). A protein is simulated **free in solution**
(`baseline`), **attached by one end to a neutral surface** (`grafted`) or **to a charged
surface** (`grafted_charged`), or as a **brush**: many copies of it grafted side by side
at a chosen grafting density (`brush`, `brush_charged`). Runs can hold thousands of
proteins on the cluster, or a few on your own computer.

You use it with one command and one settings file. No programming needed.

```bash
noodle new <name>        # make a run folder (it gets a number too)
noodle launch <run>      # start it (and later, continue anything that stopped)
noodle status [run]      # see how it's going (no run: every run, with its number)
noodle analyze <run>     # gather the results into tables
noodle export <run>...   # copy the tables of one or more runs into exports/<name>/
noodle cancel <run>      # take its jobs out of the SLURM queue
noodle clean <run>       # delete trajectories once they are on the NAS, keeping the results
noodle help [command]    # the list of commands, or one explained in detail
```

**`<run>`** is the run's name, its number, or nothing at all when you are inside its
folder (`runs/<name>/`). `noodle` works from any folder, and Tab completes commands
and run names, once you have added two lines to your `~/.bashrc` (section 1). Without
them, type `./noodle` inside the `noodle-md` folder, with run names in full.

```
noodle-md/
  noodle                   the command
  noodle-completion.bash   Tab completion for it (section 1)
  physics.yaml             the named physics sets (section 4)
  setup.sh                 one-time setup (section 1)
  sequences/               the protein sets runs draw from (section 2, step 2)
  src/                     the code
  runs/                    your runs, one folder each, and runs.csv, the list of them
                           (made by noodle new; not part of git)
  exports/                 what noodle export makes (not part of git)
  env/                     the tool's own Python (made by setup.sh; not part of git)
  CHANGELOG.md             what changed between versions
```

---

## 1. Setup (once per computer)

**Get the code.** On the cluster, log in first (`ssh your_username@cluster.address`).
On your own computer (Linux, macOS, or Windows with WSL) you also need conda: if
`conda --version` says "command not found", install
[Miniforge](https://github.com/conda-forge/miniforge).

```bash
git clone <URL of this repository> noodle-md
cd noodle-md
```

**Build the environment** (the tool's own Python, with HOOMD-blue, in the folder `env/`):

| on the cluster | on your own computer |
|---|---|
| `sbatch setup.sh` (15-30 min on a compute node; follow it with `tail setup.log`) | `bash setup.sh` |

**It worked if** the end of the output (or of `setup.log`) says
`All physics checks passed.` and `Setup finished`.

If it didn't:

| message | fix |
|---|---|
| `Could not load the module 'anaconda'` | find your cluster's name for it with `module avail conda`, put it in the line `CONDA_MODULE=...` at the top of `setup.sh` |
| connection or download errors in `setup.log` | the compute nodes have no internet: run `bash setup.sh` on the login node instead |
| `./env already exists` | a previous attempt left a half-built environment: `rm -rf env`, then try again |
| sbatch asks for a partition or account | add it to the command: `sbatch --partition=<name> setup.sh` |

**Then add two lines to your `~/.bashrc`.** The end of the setup output prints them
with your folder filled in; they look like this:

```bash
alias noodle='/home/you/noodle-md/noodle'
source '/home/you/noodle-md/noodle-completion.bash'
```

Open `~/.bashrc` (`nano ~/.bashrc`), paste them at the end, save, and open a new
terminal. Now `noodle` works from any folder, and Tab completes: `noodle st<Tab>`
gives `noodle status`, and `noodle status gr<Tab>` the run names starting with `gr`.

That's all. You never need to load or activate conda again.

---

## 2. Running simulations

### Step 1: make a run

```bash
noodle new grafted_charged_top500
```

```
Created runs/grafted_charged_top500/, run number 7
```

Pick a name that says what it is (letters, digits, `_ - .`; not only digits). This
creates `runs/grafted_charged_top500/` with a settings file, `run.yaml`, and gives the
run a **number**. From now on `noodle launch 7` is the same as
`noodle launch grafted_charged_top500`, and so is `noodle launch` typed inside
`runs/grafted_charged_top500/`. Numbers are never reused, even after a run is deleted.

### Step 2: choose your sequences

Two protein sets come with the code, in `sequences/`. A run uses one of them by naming
it in `run.yaml` (step 3); there's nothing to copy.

| file | what | proteins | lengths |
|---|---|---|---|
| `human_idrs_12k.csv` | IDRs and whole IDPs of the human proteome | 12,000 | 20-599 (median 54) |
| `mobidb_idps_15k.csv` | IDPs from MobiDB (a subsample of its 77k disordered proteins) | 15,000 | 20-500 (median 190) |

Both have the protein's name in the column `id` and its sequence in the column
`sequence`, which is what a new `run.yaml` already says (`name_column: id`,
`sequence_column: sequence`). Their other columns are ignored by the simulations (the
run's own `manifest.csv` recomputes the sequence descriptors, step 6).

- **Both are ordered for diversity.** Each row was picked to be as different as
  possible from all the rows above it (the `radius_at_pick` column shrinks down the
  file). So `rows: 1-500` gives the 500 most mutually different sequences, and
  `rows: 501-1000` the next 500.
- **Except the end of `mobidb_idps_15k.csv`:** rows 11274-15000 (`source: random`) were
  picked at random rather than for diversity.

**Your own sequences:** put a CSV file in the run folder. It needs a column with a
**unique name** per protein and a column with the **one-letter sequence**; other columns
are ignored. Naming them `id` and `sequence` saves editing `run.yaml`. A file in the run
folder takes precedence over one of the same name in `sequences/`.

```bash
cp my_sequences.csv runs/grafted_charged_top500/
# from your own computer to the cluster:
scp my_sequences.csv your_username@cluster.address:noodle-md/runs/grafted_charged_top500/
```

### Step 3: fill in `run.yaml`

```bash
nano runs/grafted_charged_top500/run.yaml     # arrows to move, Ctrl-O Enter to save, Ctrl-X to quit
```

| setting | meaning | example |
|---|---|---|
| `description` | one line about the run, shown by `noodle status` | `"Top-500 human IDRs, -0.1 e/nm^2 wall"` |
| `sequences_csv` | a set from `sequences/`, or your CSV in the run folder | `human_idrs_12k.csv` |
| `name_column`, `sequence_column` | the CSV columns with the names and the sequences; already right for the shared sets | `id`, `sequence` |
| `rows` | `all`, one row (`7`) or a range (`101-200`); row 1 is the first protein | `1-500` |
| `max_residues` | skip proteins longer than this (`0` = no limit); keeps brush runs affordable (section 3) | `200` |
| `physics` | one physics set from `physics.yaml`, or a list of several (section 3) | `grafted_charged` or `[baseline, brush]` |
| `replicates` | independent copies (different random seeds) of each protein | `3` |
| `time_limit_days` | cluster time limit per job (step 4); a task that hits it saves a checkpoint and continues at the next launch | `3` |
| `max_running` | most jobs running at once, each on one core (0 = no limit) | `0` |
| `extra_sbatch` | anything else your cluster needs | `["--partition=normal"]` |

The run will have *proteins x physics sets x replicates* **tasks**. Each task is one
simulation.

### Step 4: launch

```bash
noodle launch grafted_charged_top500
```

First it checks your settings (a misspelled setting is reported with a suggestion,
and nothing is launched). Then it shows the plan and asks before doing anything:

```
New run: grafted_charged_top500   "Top-500 human IDRs, -0.1 e/nm^2 wall"

Sequences  sequences/human_idrs_12k.csv, rows 1-500: 500 proteins, 20 to 596 residues (median 46)
Physics
  grafted              grafted by the N-terminus to a repulsive, neutral lj93 surface; 298 K, 150 mM salt
  grafted_charged      grafted by the N-terminus to a repulsive, charge -0.1 e/nm^2 lj93 surface; 298 K, 150 mM salt
                       surface potential -0.71 kT/e
Tasks      500 proteins x 2 physics set(s) x 3 replicate(s) = 3000 tasks
Length     steps per task (equilibration + production), for the shortest and longest protein
  grafted              874 thousand (N=20) to 1.53 billion (N=596), saving 3000 frames
  grafted_charged      874 thousand (N=20) to 1.53 billion (N=596), saving 3000 frames
Cost       ~34,607 core-hours, the longest task ~10d 13h (rough: 1 million residue-steps per second per core)
           [!] 186 task(s) need more than time_limit_days: each stops at the limit and continues at the
           next ./noodle launch (4 launches for the longest). Consider max_residues, or a longer time limit
Disk       ~11.8 GB of trajectories
Cluster    3000 tasks in 1210 SLURM job(s) of up to ~4 h (short tasks share a job), time limit 3 day(s) each, no limit on how many at a time

Launch these 3000 tasks? [y/N]
```

The cost is rough: it assumes one core does a million residue-steps per second (what a
recent laptop core does with this model; a cluster core may differ by a factor ~2 either
way). It is there to catch a plan that would take years before you launch it. Once tasks
run, `noodle status` estimates the time left from their real, measured speeds. Answer
`y` to launch; anything else changes nothing. Launching also writes `manifest.csv`
(step 6), with the sequence descriptors of every protein. Then:

- **Cluster:** the tasks are grouped into **jobs**, and the jobs are submitted to SLURM
  as job arrays, named after the run. The command checks that they arrived in the queue.
  You can log out.
- **Own computer:** the jobs run one after another in the terminal. Ctrl-C stops them
  (the running task saves a checkpoint first); `noodle launch` again continues.

**Jobs.** A job is one SLURM job on one core. It runs its tasks one after another, and
analyses each one as soon as it finishes (step 6). Short tasks are packed together into
jobs of up to about 4 hours (estimated; at most half of `time_limit_days`), and a longer
task gets a job of its own. So a run of thousands of minute-long tasks becomes a few
hundred jobs, not thousands. The longest jobs start first. A task that fails doesn't
stop its job: the job goes on with its next task.

Launching **freezes** the run: the task list (with every sequence written into it) and
a copy of the physics settings are saved in the run folder. After that only
`time_limit_days`, `max_running` and `extra_sbatch` can still be changed; for anything
else, make a new run. Changing a file in `sequences/` never affects a launched run.

> **Try it first.** For a 1-minute test of every geometry: `noodle new test`,
> then in `run.yaml` keep `sequences_csv: human_idrs_12k.csv` and set `rows: 1-2`,
> `physics: [test_bulk, test_grafted, test_grafted_charged, test_brush]`.

### Step 5: check on it, and continue what stopped

```bash
noodle status grafted_charged_top500
```

```
grafted_charged_top500  Top-500 human IDRs, -0.1 e/nm^2 wall
3d 07h since launch   launched 2026-11-02 09:14 · checked 2026-11-05 16:40

████████████████████████▓░░░░░  2571 of 3000 tasks done (86%)
  2 failed · 14 stopped · 9 running · 404 waiting · 2571 done

Time left  ~92d 14h of computing for the 427 unfinished tasks, the longest needs ~2d 06h
           from the speeds measured by 2589 tasks of this run; how long you wait also depends on the queue

                      failed  stopped  running  waiting       done
  baseline                 1        6        4      189  1300/1500
  grafted_charged          1        8        5      215  1271/1500

FAILED 2  fix the cause first · full messages in runs/grafted_charged_top500/jobs/job_<number>.out
  2 × RuntimeError: Particle out of box: particle 37 at (12.1, -3.0, 410.2)
      tasks 1207, 2214  (jobs 31, 57)

STOPPED 14  ./noodle launch grafted_charged_top500 resumes them from their last checkpoint
   412  Q8XYZ1_idr3   N=412   grafted_charged  stopped at  82% after 2d 23h of run time (time limit), needs ~14h 20m more
   ...

RUNNING 9
   409  P1234_idr1    N=288   grafted_charged  production     33%  21h 40m so far, ~1d 20h left
   ...

WAITING 404  in the SLURM queue, or next in a job that is running

DONE 2571  run time 9m (N=20) to 2d 19h (N=573), 812d 05h in total

Next: ./noodle launch grafted_charged_top500 resumes the 14 stopped tasks (and asks before retrying the failed ones).
```

Tasks that need you come first. Failed tasks are grouped by error, since one cause
often fails many tasks. On screen each list shows at most 10 tasks. With several
physics sets, the counts are also shown per set (done as done/total). A running task
that has not reported for over 1.5 hours is flagged `stuck?`.

**`runs/<name>/log.csv`** has every task, one row each, written by each `noodle status`:
open it in a spreadsheet to sort and filter. Its columns: `task_id`, `sequence_name`,
`physics`, `replicate`, `N`, `state` (as below), `percent`, `phase` (equilibration,
production, done), `job` and `slurm_job` (which job runs it), `steps_per_second`,
`run_time_h`, `time_left_h`, `attempts`, `host`, `started_at`, `last_update`,
`completed_at`, `stop_reason`, `error`.

Time estimates come only from speeds measured in this run, on the cluster's own nodes:
each task's speed (measured a minute after it starts, then every hour) and, for finished
tasks, all their steps over their run time. A task that has not measured its own speed
yet uses that of the measured task closest in length, in the same physics set. Until any
task has measured its speed, the time left shows as not known. "Time left" is computing
time; how long you actually wait also depends on how many jobs the queue runs at once.

| state | meaning |
|---|---|
| `waiting` | its job is in the SLURM queue, or is running the tasks before it |
| `running` | running now; progress is updated every hour |
| `stopped` | not finished and not running: it reached the time limit, was cancelled or Ctrl-C'd (`stop_reason` in `log.csv` says which), or its job ended before reaching it |
| `failed` | ended with an error, shown in the report (full message: `runs/<name>/jobs/job_<number>.out`) |
| `done` | finished (and analysed by its job) |

**Stopped tasks:** run `noodle launch <run>` again. It first tells you what it will do:

```
3000 tasks: 2571 done · 413 in the queue (left alone) · 16 to start (14 resume from a checkpoint, 2 from the beginning)
```

Finished tasks and tasks still in the queue are left alone, so relaunching is always
safe. The tasks to start are grouped into new jobs. Nothing is lost when a task stops at
the **time limit**: SLURM warns it 10 minutes ahead, and it saves a checkpoint and stops
cleanly. The same happens when you cancel it (below) or press Ctrl-C. Only a crash (a
node failing) loses work, at most an hour, back to its last hourly checkpoint. (A
resumed task is just as valid as an uninterrupted one, but not frame-for-frame
identical.)

**Failed tasks** would most likely fail again, so `noodle launch` lists their errors
and asks before starting them again. Fix the cause first.

**Cancelling:** `noodle cancel <run>` counts the run's jobs in the queue, asks, and
cancels them. Running tasks save a checkpoint before they stop, so `noodle launch`
continues the run later, as above. On your own computer, press Ctrl-C in the terminal
the run is going in.

**Every run at once:** `noodle status` with no run (outside a run folder) lists every
run with its number:

```
    #  run                      progress                  description
    5  bulk_human_top2000       2000 of 2000 done         Top-2000 human IDRs, 3 replicates
    7  grafted_charged_top500   2571 of 3000 done, 2 failed   Top-500 human IDRs, -0.1 e/nm^2 wall
    8  brush_scan_test          not launched yet          Density scan, 4 short IDRs
```

The same list, with more columns, is in **`runs/runs.csv`**: `number`, `name`,
`description`, `created`, `launched`, `finished`, `proteins`, `physics_sets`,
`replicates`, `tasks`, `jobs`, `done`, `failed`, `analysed`, `core_hours` (computing
used so far), `trajectories_gb` (disk used now), `last_checked`. It is updated by
`noodle new`, `launch`, `status` and `analyze`. A deleted run keeps its row, so that
its number is never given to another run.

### Step 6: results

Every job **analyses each of its tasks right after it finishes**, on the cluster node,
from the task's trajectory (seconds for a chain, about half a minute for a brush). Then

```bash
noodle analyze grafted_charged_top500
```

gathers those results into tables. It takes seconds, so run it as often as you like,
also while the run is still going:

```
2571 of 3000 tasks analysed (2571 completed, 2 failed)
  grafted              1290 of 1500 proteins   Rg 0.9123 to 6.402 nm (median 1.731) · centre-of-mass height 1.012 to 9.88 nm (median 2.243)
  ...
  [!] 3 protein x physics row(s) were still drifting (column drift): equilibrate longer, e.g. ...
  -> manifest.csv, analysis/replicates.csv, analysis/profiles.csv, analysis/report.json
```

| file | what |
|---|---|
| `runs/<name>/manifest.csv` | **the result**: one row per protein and physics set: its sequence, the conditions, the sequence descriptors, and every observable averaged over its replicates, with error bars |
| `analysis/replicates.csv` | the same observables, one row per task (replicate), to see how much replicates differ |
| `analysis/profiles.csv` | quantities along a coordinate (height, residue, separation along the chain), averaged over replicates |
| `analysis/report.json` | an overview: task counts, the smallest, median and largest value of every observable per physics set, and the rows worth a look (still drifting, too few independent samples, could not be analysed) |

All are plain CSV and JSON: open them in a spreadsheet, or in Python with
`pandas.read_csv`. To work on them on your own computer, use `noodle export` (step 7).

#### `manifest.csv`, column by column

Lengths are in nm. Columns that don't apply to a row stay empty (heights for a free
chain, brush columns for a single chain).

**Who and where**, written at launch:

| column | meaning |
|---|---|
| `name`, `physics`, `sequence` | the protein, its physics set, its sequence |
| `source_csv`, `source_row` | where the sequence came from |
| `geometry`, `temperature_K`, `ionic_concentration_M` | conditions (from the physics set) |
| `graft_terminus`, `wall_charge_density`, `wall_attraction`, `grafting_density` | grafted only: which end is attached, surface charge (e/nm²), stickiness, chains per nm² (0 for a single chain) |
| `chains` | chains per simulation: 1, or a brush's 100 |

**Sequence descriptors**, computed from the sequence at launch with the model's own
charges and hydropathies (HPS: histidine counts +0.5), with the same formulas and names
as idp-rg, so the tables join directly (`src/descriptors.py`):

| column | meaning |
|---|---|
| `N` | residues |
| `net_charge`, `ncpr` | sum of the charges (e), and per residue |
| `fcr` | fraction of charged residues (q ≠ 0, histidine included) |
| `f_plus`, `f_minus` | positive and negative charge per residue |
| `mean_hydro` | mean hydropathy λ |
| `kappa` | how mixed (0) or separated into blocks (1) the + and − charges are (Das & Pappu 2013); empty without charges |
| `scd`, `scd_pos`, `scd_neg` | sequence charge decoration (Sawle & Ghosh 2015), negative when opposite charges sit in blocks; `scd_pos`/`scd_neg` are asinh(5 scd / N) split at 0, as in idp-rg |
| `shd` | sequence hydropathy decoration (Zheng et al. 2020) |

`kappa` compares the sequence with the most segregated block arrangement of the same
charges. For ~2% of sequences, mostly with few charges (`fcr` < 0.2), spreading the
charges out beats every block arrangement, and `kappa` comes out above 1; there it says
little.

**Results**, written by `noodle analyze`. First `replicates` (how many were analysed),
then every observable, each followed by its error bar `<name>_err` (one standard error:
the scatter between the replicates' own averages divided by √replicates; with a single
replicate, a block-average error from its frames). For a brush, each value is the
average over all its chains.

Every chain:

| column | meaning |
|---|---|
| `rg_nm`, `rg_sq_nm2`, `rg_sd_nm` | radius of gyration (mass-weighted), ⟨Rg²⟩, and the width of the Rg distribution |
| `rh_nm` | hydrodynamic radius, Kirkwood double sum, 1/Rh = ⟨Σ 1/r_ij⟩ / N² over pairs i ≠ j (see the caveats below) |
| `ree_nm`, `ree_sq_nm2` | end-to-end distance and ⟨Ree²⟩ |
| `asphericity` | (largest gyration moment − mean of the other two) / Rg²: 0 sphere, 0.25 flat disc, 1 rod |
| `shape_anisotropy` | relative shape anisotropy κ² (not the sequence `kappa`): 0 sphere, 0.25 flat disc, 1 rod |
| `prolateness` | below 0 flattened (oblate), above 0 elongated (prolate) |
| `transverse_rms_nm` | rms distance of the residues from the line through the two ends |
| `nu`, `r0_nm` | scaling exponent of the distances within the chain, R_ij = r0 \|i − j\|^ν, fitted for \|i − j\| ≥ 5: 0.5 ideal chain, 0.588 good solvent, 1/3 collapsed globule |
| `bond_cos` | ⟨cos θ⟩ between consecutive bonds |
| `persistence_nm` | persistence length, Flory's definition: the end-to-end vector projected on the end bond, averaged over both ends |
| `contacts_per_residue` | residues of the same chain, at least 3 apart along it, closer than 1.5 × their mean diameter |

Grafted chains and brushes (heights above the surface):

| column | meaning |
|---|---|
| `com_height_nm` | height of the chain's centre of mass |
| `end_height_nm` | height of the free end |
| `max_height_nm` | height of its highest residue (averaged over frames) |
| `height_95_nm` | height below which 95% of its residues are |
| `rg_perp_nm`, `rg_par_nm` | Rg perpendicular (z) and parallel (x, y) to the surface |
| `end_lateral_rms_nm` | rms sideways distance of the free end from the anchor |
| `tilt_cos` | ⟨cos⟩ of the angle between anchor → free end and the surface normal: 1 straight up, 0 lying flat |
| `bond_order` | ⟨P2⟩ of the bonds with the surface normal: 1 all vertical, 0 random, −0.5 all flat |
| `adsorbed_fraction` | fraction of its residues (the anchor left out) touching the surface: within their contact distance + 0.5 nm |

Brushes:

| column | meaning |
|---|---|
| `brush_height_nm` | twice the mean residue height (the first moment of the density profile) |
| `concentration_mg_ml` | protein concentration below the brush's `height_95` |
| `interchain_contacts_per_residue` | contacts (as above) with residues of other chains |
| `lateral_heterogeneity` | coefficient of variation of the residue counts in a grid of cells one grafting spacing wide: 0 even, larger when chains bundle together |
| `reduced_density` | grafting density × π ⟨Rg²⟩, with the ⟨Rg²⟩ of the same protein as a single chain at the same temperature and salt in this run (a bulk set, else a single grafted chain; empty if the run has neither). Below ~1 the chains barely touch ("mushrooms"); well above 1 it is a true brush |

Diagnostics (is the run long enough?):

| column | meaning |
|---|---|
| `tau_rg_tu`, `tau_ree_tu` | correlation time of Rg and of the end-to-end vector, in model time units (average over replicates) |
| `independent_samples` | independent Rg values per chain, all replicates together: frames / (2 × correlation time). **Below 20 the run was too short to trust** (`noodle analyze` warns); hundreds is comfortable |
| `drift` | how many replicates were still drifting: Rg (or, grafted, the centre-of-mass height) differs between the two halves of production by more than 3 error bars. Should be 0 |

#### `replicates.csv` and `profiles.csv`

`replicates.csv` has one row per analysed task: `task_id`, `name`, `physics`,
`replicate`, `seed`, `N`, `chains`, `frames`, `run_time_h`, then the observables above
with block-average error bars, and the diagnostics of that replicate alone.

`profiles.csv` is a long table: `name`, `physics`, `profile`, `x`, `x_unit`, `value`,
`value_unit`, `err` (standard error between replicates; empty with one), `replicates`.
To plot one: `df[(df.name == "P12345") & (df.physics == "brush") & (df.profile == "density")]`.

| profile | x | value |
|---|---|---|
| `internal_distance` | separation \|i − j\| along the chain (about 25 values, log-spaced) | √⟨R_ij²⟩, nm (the curve `nu` is fitted to) |
| `bond_correlation` | bonds apart (log-spaced) | ⟨cos⟩ between the two bonds |
| `residue_contacts` | residue (1 = N-terminus) | its contacts within the chain |
| `density`, `mass_density`, `charge_density`, `end_density` | grafted: height above the surface (bins of 0.25 nm) | residues, mass, charge and free ends per volume (a brush: per nm³ of box, mass in mg/mL) or per nm of height (a single chain) |
| `residue_height` | grafted: residue | its mean height, nm |
| `residue_adsorption` | grafted: residue | probability it touches the surface |
| `residue_interchain_contacts` | brushes: residue | its contacts with other chains |

#### Caveats

- **Times are model times.** The Langevin thermostat's friction sets how fast the chain
  moves, and it is not the friction of water: `tau_rg_tu` and `tau_ree_tu` are right for
  comparing sequences and conditions with each other, not as real nanoseconds.
- **Rh** from the Kirkwood formula is a known approximation for flexible chains; compare
  trends between proteins rather than absolute values with DLS or NMR.
- **Brush run lengths are provisional** (section 3): check `drift` and
  `independent_samples` of brush rows before trusting them.

**If a job could not analyse a task** (the error is in `analysis_error` of its job
file), `noodle analyze` analyses it itself, from the trajectory, into
`analysis/reanalysed.csv`. Rows it still can't analyse are listed in `report.json`.

The trajectories (`runs/<name>/trajectories/*.gsd`) open in [OVITO](https://www.ovito.org),
in Python with the `gsd` package, or with MDAnalysis (installed in `env/`) for anything
the tables don't have.

### Step 7: take the results home

```bash
noodle export grafted_charged_top500
noodle export 5 7            # several runs together
```

It first brings each run's analysis up to date, then writes `exports/<name>/`:

| file | what |
|---|---|
| `manifest.csv`, `replicates.csv`, `profiles.csv` | the runs' tables, one after another, with a first column `run` |
| `runs.csv` | the runs' rows of `runs/runs.csv` |
| `physics/<run>.yaml` | the physics sets each run used, fully written out |

For one run the folder is named after it; for several, it asks, suggesting what their
names have in common. It warns about runs that aren't finished, and about a physics set
name with different settings in different runs. Copy the folder to your own computer:

```bash
rsync -av your_username@cluster.address:noodle-md/exports/grafted_charged_top500/ ./grafted_charged_top500/
```

### Step 8: free the disk space

Trajectories take most of a run's space (megabytes per task). Once a run is finished
and analysed, move the trajectories to the group NAS and delete them from the cluster:

1. Copy `runs/<name>/trajectories/` to the NAS with Globus. (Globus checks every file
   it transfers.)
2. Run

   ```bash
   noodle clean grafted_charged_top500
   ```

   It first gathers the results (as `noodle analyze`), then shows what it would delete
   and asks you to type the run name:

   ```
   Delete     the trajectories of 2571 completed, analysed task(s): 9.8 GB
   Keep       analysis/ (180 MB), trajectories.sha256, jobs/, manifest.csv, tasks.csv, physics.yaml, run.yaml
              the files of unfinished tasks (1.3 GB): needed to resume them
   ```

Only trajectories of tasks that are **finished and analysed** are deleted (plus any
checkpoint a finished task left behind). Unfinished tasks keep their trajectory and
checkpoint, so `noodle launch` can still resume them; run `noodle clean` again
once they finish. Everything small stays, and the results live in `jobs/`, so
`noodle status`, `analyze` and `export` keep working.

Before deleting, it writes `runs/<name>/trajectories.sha256`, the checksum of every
trajectory it deletes. It is the record of exactly which files went to the NAS, and
lets you check a copy at any time: put it in the copy's folder and run
`sha256sum -c trajectories.sha256` (each file should say `OK`).

To compute something new from the trajectories later, copy them back from the NAS into
`runs/<name>/trajectories/`.

---

## 3. Bulk, grafted, charged wall, brushes, and batches

The geometry is just the physics set you name in `run.yaml`:

```yaml
physics: baseline                              # free chains
physics: grafted                               # grafted by the N-terminus to a neutral surface
physics: grafted_reversed                      # grafted by the C-terminus to a neutral surface
physics: grafted_charged                       # grafted to a surface with -0.1 e/nm^2
physics: brush                                 # a 10 x 10 patch of chains grafted at 0.05 chains/nm^2
physics: brush_charged                         # the same brush on a -0.1 e/nm^2 surface
physics: [baseline, grafted, grafted_charged]  # all three, for every protein, in one run
```

### Brushes

A brush set grafts a **10 x 10 patch** of copies of the protein (`chains: 100`) on a square
lattice, `grafting_density` chains per nm² apart, in a box that repeats sideways: each
chain sees neighbours all around it, as in an infinite brush. The grafting density is the
setting you choose; the box follows from it: it is sqrt(100 / grafting_density) wide
(0.01 chains/nm²: 100 nm, 0.05: 44.7 nm, 0.3: 18.3 nm, 1: 10 nm), and as tall as for a
single grafted chain (the fully stretched chain plus `wall_headroom`, then one interaction
range of empty space). Every other setting (surface, charge, attraction, `graft_terminus`)
works as for a single grafted chain. A box narrower than twice the interaction range (8 nm)
would let a residue feel the same neighbour twice through the periodic boundary; a
10 x 10 patch never gets there (grafting_density is at most 1), but a set with fewer
`chains` does at high density, and then gets more chains automatically (the launch plan
says so).

**Which densities matter** depends on the protein's size: the chains start to touch
around grafting density × π Rg² ≈ 1. For a 100-residue IDP (Rg ≈ 3 nm) that is ~0.03
chains/nm², for a 50-residue one ~0.1. A scan from ~0.01 to ~0.3 chains/nm² covers the
mushroom regime, the crossover and the brush. Put a single-chain set in the same run
(`physics: [baseline, brush, ...]`), so `noodle analyze` can give each brush its
`reduced_density`. A density scan is a batch, like the charges below:

```yaml
# in physics.yaml
brush_s001:
  extends: brush
  description: "Brush, 0.01 chains/nm^2"
  grafting_density: 0.01

brush_s02:
  extends: brush
  description: "Brush, 0.2 chains/nm^2"
  grafting_density: 0.2
```

```yaml
# in runs/<name>/run.yaml
physics: [baseline, brush_s001, brush, brush_s02]
max_residues: 200
```

**How long a brush must run is not settled yet.** The `brush` set's run length
(`rouse_multiplier: 2`, `production_multiplier: 3`, `frames: 500`) is provisional: it comes
from relaxation times measured earlier on 4 x 4 brushes, not yet on 10 x 10 ones. Until
it is, check the `drift` and `independent_samples` columns of the manifest for every
brush, and expect dense brushes (~0.3 chains/nm²) to need a longer production.

**What it costs.** A brush task has 100x the residues of a single chain, and run length
grows with the chain's relaxation time (~N^2.2), so its cost grows like **N^3.2**. Every
frame averages over 100 chains, which is why a brush gets by with a much shorter
production than a single chain (`production_multiplier` 3 instead of 22). With the
provisional run length, on one core (rough: 1 million residue-steps per second):

| residues N | brush (10 x 10) | single chain | brush trajectory |
|---|---|---|---|
| 50 | 1.9 h | 5 min | 30 MB |
| 100 | 17 h | 50 min | 60 MB |
| 150 | 2.6 days | 3 h | 90 MB |
| 200 | 6.7 days | 8 h | 120 MB |
| 300 | 24 days | 1.2 days | 180 MB |
| 400 | 61 days | 3 days | 240 MB |

So **set `max_residues` for brush runs**. One density, one replicate, of the shared sets:

| set, max_residues | proteins | core-hours | longest task |
|---|---|---|---|
| human_idrs_12k, 150 | 9,801 (82%) | ~66,000 | 2.6 days |
| human_idrs_12k, 200 | 10,501 (88%) | ~139,000 | 6.7 days |
| mobidb_idps_15k, 150 | 6,010 (40%) | ~99,000 | 2.6 days |
| mobidb_idps_15k, 200 | 7,881 (53%) | ~299,000 | 6.7 days |

With 7,500 cores, 66,000 core-hours is ~9 hours of the whole cluster. Tasks longer than
`time_limit_days` simply continue at the next `noodle launch` (the plan warns about
them). These are estimates for one recent CPU core; `noodle status` measures the real
speed of your run.

### Charged wall

`grafted_charged` is `grafted` plus a uniform surface charge of **-0.1 e/nm²** (one
elementary charge per 10 nm²). A negative surface pulls lysine and arginine (and
histidine, +0.5) toward it and pushes aspartate and glutamate away. The charge is
screened by salt with the same Debye length as the residue-residue electrostatics (0.79
nm at 150 mM).

**Choosing the charge.** The surface is modelled with linearised Debye-Hückel theory,
which is accurate while the surface potential stays below about 1 kT/e. Launching shows
that potential for every charged set, and warns above 1:

| surface charge, e/nm² | 50 mM salt | 150 mM | 500 mM |
|---|---|---|---|
| ±0.05 | 0.61 kT/e | 0.35 | 0.19 |
| ±0.1 | 1.22 | 0.71 | 0.39 |
| ±0.2 | 2.44 | 1.41 | 0.77 |

**A batch over several charges** (or signs, or salts): add one set per value at the
end of `physics.yaml` (section 4), then list them all in one run:

```yaml
# in physics.yaml
grafted_charged_m005:
  extends: grafted_charged
  description: "Grafted, -0.05 e/nm^2 surface"
  wall_charge_density: -0.05

grafted_charged_p01:
  extends: grafted_charged
  description: "Grafted, +0.1 e/nm^2 surface"
  wall_charge_density: 0.1
```

```yaml
# in runs/<name>/run.yaml
physics: [grafted, grafted_charged_m005, grafted_charged, grafted_charged_p01]
```

Every protein then runs against every surface. Each set gets its own rows in
`manifest.csv` (column `physics`, and the condition columns say what it was).

### Several runs at once

You can have as many runs as you like going at the same time, of any geometry. Each run
lives in its own folder, its jobs are named after the run, and every job writes only
into its own run folder, so runs never interfere. A large set can be split into several
runs by `rows` (e.g. `1-1000`, `1001-2000`, ...) to launch and check it in pieces, and
put back together with `noodle export`.

---

## 4. Changing the physics

All physics settings are named sets in **`physics.yaml`**. To try something new,
**don't edit the existing sets** (other people's new runs use them). Add a set at the
end that `extends` an existing one and lists only what changes:

```yaml
grafted_sticky_negative:
  extends: grafted_charged
  description: "Weakly attractive surface, -0.05 e/nm^2"
  wall_attraction: 0.3
  wall_charge_density: -0.05
```

Then write `physics: grafted_sticky_negative` in a run's `run.yaml`. Use a new name with
no spaces, and indent with exactly two spaces. A run keeps its own frozen copy of the
sets it uses, so editing `physics.yaml` never changes a run that has been launched.
Sets only one run needs can instead go in `runs/<name>/physics_sets.yaml` (same format,
and they may `extend` the sets of `physics.yaml`, but not reuse their names).

**Any geometry**

| setting | meaning | default |
|---|---|---|
| `geometry` | `bulk` or `grafted` | `bulk` |
| `temperature_K` | temperature, K | `298` |
| `ionic_concentration_M` | salt, mol/L (screens charges) | `0.150` |
| `epsilon_lj` | strength of the attraction between residues, kcal/mol | `0.1` |
| `lj_r_cut`, `yukawa_r_cut` | range of the residue attraction and of the electrostatics, nm | `2.0`, `3.5` |
| `harmonic_k`, `harmonic_r0` | bond stiffness (kcal/mol/nm²) and length (nm) | `1920`, `0.38` |
| `box_size` | box size, nm (a single grafted chain: its width) | `500` |
| `dt` | time step, in units of 0.489 ps | `0.02` |
| `rouse_multiplier`, `production_multiplier` | equilibration and production length, see below | `2`, `22` |
| `equilibration_steps`, `production_steps` | a number here sets the length directly instead | empty |
| `frames` | frames saved, evenly spaced over production | `3000` |

**Grafted only**

| setting | meaning | default |
|---|---|---|
| `graft_terminus` | attached end: `N` (first residue) or `C` (last) | `N` |
| `wall_potential` | surface repulsion: `lj93` (smooth solid) or `lj126` (like two residues) | `lj93` |
| `wall_epsilon`, `wall_sigma` | strength (kcal/mol) and size (nm) of the residue-surface interaction | `0.1`, `0.5` |
| `wall_attraction` | 0 = surface only repels; up to 1 = increasingly sticky | `0.0` |
| `wall_attraction_scale` | `hydropathy`: hydrophobic residues stick more; `uniform`: all the same | `hydropathy` |
| `wall_r_cut` | range of the surface attraction, nm | `3.5` |
| `wall_charge_density` | surface charge, e/nm². Negative attracts K and R, repels D and E | `0.0` (`grafted_charged`: `-0.1`) |
| `wall_headroom` | space between the fully stretched chain and the ceiling of the box, nm | `2.0` |
| `grafting_density` | chains per nm². `0` = a single chain; above 0, a brush | `0` (`brush`: `0.05`) |
| `chains` | a brush's number of chains, a square number (`100` = a 10 x 10 patch); more are added if the box would be too narrow | `100` |

**Run length.** Steps = multiplier x N^2.2 / dt (N = number of residues), because a
chain's relaxation time grows like N^2.2. The defaults are a rule of thumb from bulk
chains; the `drift` and `independent_samples` columns of the manifest show whether a run
was long enough. (Measured relaxation times are ~0.01-0.04 x N^2.2 time units for single
chains and brushes up to 0.1 chains/nm², so the single-chain default covers several
hundred relaxation times. The brush default is provisional, section 3.) The `test_*`
sets use tiny fixed lengths, for trying things out only.

**Speed.** `lj_r_cut: 2.0` and `dt: 0.02` (~10 fs, the time step of the original HPS
model) were chosen for speed: the LJ tail beyond 2 nm is worth less than 0.005 kT per
residue, and dt 0.02 gives the same temperature and energies as 0.01. Together with
skipping the electrostatics of uncharged pairs (exactly zero anyway), runs are ~3x
(single chains) to ~4x (brushes) faster than with the earlier `lj_r_cut: 3.5`, `dt: 0.01`.
Runs launched before keep their own frozen settings.

---

## 5. Updating, and where things are

**Updating the code:** `git pull`. Your runs (`runs/`) and the environment (`env/`) are
not part of the code, so `git pull` never touches them. If `src/environment.yml`
changed, rebuild: `rm -rf env`, then section 1 again. `CHANGELOG.md` says what changed.
A run launched with one version can't always be continued by the next (when the run
folder's layout changed, the CHANGELOG says so): finish it with the version that
launched it, e.g. `git checkout v1`, then `git checkout main` again.

**A run folder**

```
runs/<name>/
  run.yaml, <your>.csv   your settings (and your own sequences, if any)
  physics_sets.yaml      optional: physics sets only this run uses
  tasks.csv              one row per task: protein, physics set, replicate, seed, sequence
  physics.yaml           the physics sets used, frozen at the first launch
  manifest.csv           one row per protein and physics set: descriptors, conditions, results (step 6)
  log.csv                every task's state, written by noodle status (step 5)
  jobs/
    job_0007.csv           job 7's tasks: how each is going, and its results once analysed
    job_0007_profiles.csv  the profiles of job 7's tasks
    job_0007.out           everything job 7 printed (errors in full)
    submitted.csv          which SLURM job runs which job
    slurm.out              SLURM's own messages (normally empty)
  analysis/              the noodle analyze tables (step 6)
  trajectories/          the saved frames, one .gsd file per task (until noodle clean)
  checkpoints/           save points of unfinished tasks
  trajectories.sha256    checksums of the trajectories noodle clean deleted
```

Each file has a single writer at any time (launch, a job, status, analyze), so commands
on the login node and jobs on the cluster never overwrite each other. A task's current
state is its row in the newest job file that lists it; a task that is continued or
retried gets a row in a new job.

**The code** (start reading at `src/forcefield.py`, whose header explains how HOOMD
represents forces)

```
noodle                  the command: starts the matching program below
noodle-completion.bash  Tab completion for it
src/forcefield.py       the model: residue table, electrostatics, every force
src/simulate.py         one job: for each of its tasks, build the chain or brush, equilibrate,
                        record, checkpoint, then analyse it
src/observables.py      what a task measured: every column of the results
src/descriptors.py      the sequence descriptors of the manifest
src/launch_run.py       run.yaml + CSV -> tasks -> jobs; submits them to SLURM or runs them here
src/new_run.py          noodle new          src/run_template.yaml   the run.yaml it copies
src/status.py           noodle status       src/analyze_run.py      noodle analyze
src/export_runs.py      noodle export       src/cancel_run.py       noodle cancel
src/clean_run.py        noodle clean
src/runs.py             shared helpers      src/slurm_job.sh        the SLURM job script
src/check_physics.py    noodle check: the physics self-test (setup.sh runs it)
src/environment.yml     the packages setup.sh installs into env/
```

**GPU workstation:** change `cpu*` to `gpu*` in `src/environment.yml` before `setup.sh`.

---

## 6. Troubleshooting

| problem | fix |
|---|---|
| `The environment is missing` | setup hasn't been done on this computer, or failed (see `setup.log`) |
| `noodle: command not found` | the alias isn't in your `~/.bashrc` yet, or the terminal was open before you added it (section 1). Meanwhile: `./noodle` inside the `noodle-md` folder |
| Tab doesn't complete | the `source .../noodle-completion.bash` line isn't in your `~/.bashrc`, or the terminal is older than it |
| `Permission denied` when running `./noodle` | `chmod +x noodle` |
| `There is no run called ...` | check the spelling or the number; `noodle status` lists the runs |
| `Can't find the sequence file` / `has no column` | the CSV must be in the run folder or in `sequences/`, and its file and column names must match `run.yaml` exactly (the column names of the shared sets are in section 2, step 2) |
| `physics.yaml has no set called ...` | typo in `physics:`; the message lists the sets |
| `[!] above ~1 kT/e` when launching a charged set | the surface charge is too high for the linearised model at this salt: see section 3 |
| changed `run.yaml` but nothing happened | after the first launch, only `time_limit_days`, `max_running`, `extra_sbatch` count; make a new run |
| delete a run | `noodle cancel <run>` first, then `rm -rf runs/<name>` (its number stays taken in `runs/runs.csv`) |

---

## 7. The model

- **Residues:** one bead each, with its own mass, charge, size and hydropathy λ (0 to 1),
  from the HPS model (Dignon et al., *PLoS Comput. Biol.* 2018).
- **Residue pairs:** Lennard-Jones attraction scaled by the pair's average λ
  (Ashbaugh-Hatch form), plus screened electrostatics (Debye-Hückel). Both the Bjerrum
  length and the Debye length follow from the temperature (through water's dielectric
  constant) and the salt concentration.
- **Bonds:** harmonic springs. **Temperature:** Langevin thermostat.
- **Grafted:** a smooth surface at the bottom of the box, with the anchor residue fixed
  at its contact distance from it. A brush has one anchor per chain, on a square
  lattice, in a box periodic in x and y; its chains start grown upward each in its own
  lattice column, so they never overlap at the start. The surface has a repulsive core, an optional
  attraction scaled like hydropathy, and an optional charge, screened like the chain
  (linearised Poisson-Boltzmann, U(h) = 4π l_B kT q σ / κ · e^(−κh) for a residue of
  charge q at height h). A repulsive ceiling, out of the chain's reach, stops beads
  wrapping through the box.
- **Seeds:** every task has its own random seed (in `tasks.csv`), and the replicates of a
  protein under a physics set always have different ones.

`noodle check` (also run by `setup.sh`) checks these against their formulas, runs short
bulk, grafted, charged-wall and brush simulations as one job, and checks their analysis.
Every formula is explained in `src/forcefield.py`.
