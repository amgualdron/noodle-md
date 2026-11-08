# noodle-md

Simulations of intrinsically disordered proteins (IDPs). Each amino acid is one bead,
and the beads form a chain (the HPS model). A chain is simulated in one of three ways:
**free in solution** (`baseline`), **attached by one end to a neutral surface**
(`grafted`), or **attached to a charged surface** (`grafted_charged`). Runs can hold
thousands of proteins on the cluster, or a few on your own computer.

You use it with one command and one settings file. No programming needed.

```bash
./noodle new <name>        # make a run folder
./noodle launch <name>     # start it (and resume anything that stopped)
./noodle status <name>     # see how it's going
./noodle analyze <name>    # turn finished simulations into CSV tables
./noodle clean <name>      # delete trajectories once they are on the NAS, keeping the results
```

Every command in this README is typed in a terminal **inside the `noodle-md` folder**.

```
noodle-md/
  noodle         the command
  physics.yaml   the named physics sets (section 4)
  setup.sh       one-time setup (section 1)
  sequences/     the protein sets runs draw from (section 2, step 2)
  src/           the code
  runs/          your runs, one folder each (made by ./noodle new; not part of git)
  env/           the tool's own Python (made by setup.sh; not part of git)
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

That's all. You never need to load or activate conda again.

---

## 2. Running simulations

### Step 1: make a run

```bash
./noodle new grafted_charged_top500
```

Pick a name that says what it is (letters, digits, `_ - .`). This creates
`runs/grafted_charged_top500/` with a settings file, `run.yaml`.

### Step 2: choose your sequences

Two protein sets come with the code, in `sequences/`. A run uses one of them by naming
it in `run.yaml` (step 3); there's nothing to copy.

| file | what | proteins | lengths |
|---|---|---|---|
| `human_idrs_12k.csv` | IDRs and whole IDPs of the human proteome | 12,000 | 20-599 (median 54) |
| `mobidb_idps_15k.csv` | IDPs from MobiDB (a subsample of its 77k disordered proteins) | 15,000 | 20-500 (median 190) |

Both have the protein's name in the column `id` and its sequence in the column
`sequence`, which is what a new `run.yaml` already says (`name_column: id`,
`sequence_column: sequence`). Their other columns (`fcr`, `scd_pos`/`scd_neg`,
`mean_hydro`, ...) describe each sequence and are ignored by the simulations.

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
| `description` | one line about the run, shown by `./noodle status` | `"Top-500 human IDRs, -0.1 e/nm^2 wall"` |
| `sequences_csv` | a set from `sequences/`, or your CSV in the run folder | `human_idrs_12k.csv` |
| `name_column`, `sequence_column` | the CSV columns with the names and the sequences; already right for the shared sets | `id`, `sequence` |
| `rows` | `all`, one row (`7`) or a range (`101-200`); row 1 is the first protein | `1-500` |
| `physics` | one physics set from `physics.yaml`, or a list of several (section 3) | `grafted_charged` or `[grafted, grafted_charged]` |
| `replicates` | independent copies (different random seeds) of each protein | `3` |
| `time_limit_days` | cluster time limit per task; a task that hits it can be resumed (step 5) | `3` |
| `max_running` | most tasks running at once (0 = no limit) | `0` |
| `extra_sbatch` | anything else your cluster needs | `["--partition=normal"]` |

The run will have *proteins x physics sets x replicates* tasks. Each task is one
simulation.

### Step 4: launch

```bash
./noodle launch grafted_charged_top500
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
  grafted              1.75 million (N=20) to 3.06 billion (N=596), saving 3000 frames
  grafted_charged      1.75 million (N=20) to 3.06 billion (N=596), saving 3000 frames
Time       estimated by ./noodle status once the first tasks have measured their speed
Disk       ~11.8 GB of trajectories
Cluster    3000 SLURM jobs, time limit 3 day(s) each, no limit on how many at a time

Launch these 3000 tasks? [y/N]
```

There is no time estimate yet: a test on your laptop or the login node says little about
the cluster's nodes, so `./noodle status` estimates it from the run itself. Answer `y` to launch;
anything else changes nothing. Then:

- **Cluster:** every task is submitted to SLURM as its own job, and the command checks
  that they arrived in the queue. You can log out.
- **Own computer:** the tasks run one after another in the terminal. Ctrl-C stops them;
  `./noodle launch` again continues.

Launching **freezes** the run: the task list (with every sequence written into it) and
a copy of the physics settings are saved in the run folder. After that only
`time_limit_days`, `max_running` and `extra_sbatch` can still be changed; for anything
else, make a new run. Changing a file in `sequences/` never affects a launched run.

> **Try it first.** For a 1-minute test of all three geometries: `./noodle new test`,
> then in `run.yaml` keep `sequences_csv: human_idrs_12k.csv` and set `rows: 1-2`,
> `physics: [test_bulk, test_grafted, test_grafted_charged]`.

### Step 5: check on it, and resume what stopped

```bash
./noodle status grafted_charged_top500
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

FAILED 2  fix the cause first · full messages in runs/grafted_charged_top500/slurm/
  2 × RuntimeError: Particle out of box: particle 37 at (12.1, -3.0, 410.2)
      tasks 1207, 2214

STOPPED 14  ./noodle launch grafted_charged_top500 resumes them from their last checkpoint
   412  Q8XYZ1_idr3   N=412   grafted_charged  stopped at  82% after 2d 23h of run time, needs ~14h 20m more
   ...

RUNNING 9
   409  P1234_idr1    N=288   grafted_charged  production     33%  21h 40m so far, ~1d 20h left
   ...

WAITING 404  in the SLURM queue

DONE 2571  run time 9m (N=20) to 2d 19h (N=573), 812d 05h in total

Next: ./noodle launch grafted_charged_top500 resumes the 14 stopped tasks (and asks before retrying the failed ones).
```

Tasks that need you come first. Failed tasks are grouped by error, since one cause
often fails many tasks. On screen each list shows at most 10 tasks; `runs/<name>/status.txt`
has the full report. With several physics sets, the counts are also shown per set (done as
done/total). A running task that has not reported for over 1.5 hours is flagged `stuck?`.
Time estimates come only from speeds measured in this run, on the cluster's own nodes:
each task's speed (measured a minute after it starts, then every hour) and, for finished
tasks, all their steps over their run time. A task that has not measured its own speed
yet uses that of the measured task closest in length, in the same physics set. Until any
task has measured its speed, the time left shows as not known. "Time left" is computing
time; how long you actually wait also depends on how many tasks the queue runs at once.

| state | meaning |
|---|---|
| `waiting` | in the SLURM queue |
| `running` | running now; progress is updated every hour |
| `stopped` | not finished and not running (time limit, cluster problem, Ctrl-C) |
| `failed` | ended with an error, shown in the report (full message: its file in `runs/<name>/slurm/`) |
| `done` | finished |

**Stopped tasks:** run `./noodle launch <name>` again. Each task saves a checkpoint
every hour, so it continues from there, losing at most an hour. It first tells you what
it will do:

```
3000 tasks: 2571 done · 413 in the queue (left alone) · 16 to start (14 resume from a checkpoint, 2 from the beginning)
```

Finished tasks and tasks still in the queue are left alone, so relaunching is always
safe. (A resumed task is just as valid as an uninterrupted one, but not frame-for-frame
identical.)

**Failed tasks** would most likely fail again, so `./noodle launch` lists their errors
and asks before starting them again. Fix the cause first.

`./noodle status` without a name lists all runs with how many tasks are done.

### Step 6: results

```bash
./noodle analyze grafted_charged_top500
```

Everything is per protein and physics set, with the replicates pooled. It writes, in
`runs/<name>/analysis/`:

| file | what |
|---|---|
| `summary.csv` | **the result**: one row per protein and physics set, averages with error bars |
| `timeseries/<physics set>/<protein>.csv.gz` | every saved frame of all its replicates (~60 KB per replicate) |

You can run it while the run is still going. Each replicate is analysed once, when it
has finished, and added to its protein's timeseries; rerunning skips what is already
there, so it is quick.

The timeseries (gzip-compressed CSV; pandas reads them directly with `pd.read_csv`) hold
the per-frame values every statistic is made from, so a new statistic (a distribution,
a correlation time, a different average) needs only them, not the trajectories:

| column | meaning (lengths in nm) |
|---|---|
| `replicate` | which replicate (independent simulation) the frame comes from |
| `frame` | the saved frame within that replicate, 0, 1, 2, ... (evenly spaced over production) |
| `rg` | radius of gyration (chain size), mass-weighted |
| `ree_x`, `ree_y`, `ree_z` | end-to-end vector, from the first residue to the last. For grafted chains z points away from the surface (so `ree_z` is negative when the C-terminus is grafted) |
| `asphericity` | shape: 0 = ball, 1 = rod |
| `com_height` | grafted only: height of the chain's centre of mass above the surface |

`summary.csv` gives, besides `n_residues`, `replicates` and `frames` (pooled), for each
of `rg`, `rg_sq` (Rg²), `ree_sq` (end-to-end distance²), `asphericity`, `ree_z` and
`com_height` (the last two for grafted chains only) the average over all pooled frames
and its error bar `<name>_err` (one standard error). The error bar comes from the
scatter between the replicates' own averages: replicates are independent, so this is
the honest one. With a single replicate it is a block error instead (its frames cut
into 10 consecutive blocks, from the scatter of the block averages).

The trajectories (`runs/<name>/trajectories/*.gsd`) open in [OVITO](https://www.ovito.org)
or in Python with the `gsd` package.

### Step 7: free the disk space

Trajectories take most of a run's space (megabytes per task), so once a run's results
are in, move the trajectories to the group NAS and delete them from the cluster:

1. Copy `runs/<name>/trajectories/` to the NAS with Globus. (Globus checks every file
   it transfers.)
2. Run

   ```bash
   ./noodle clean grafted_charged_top500
   ```

   It first analyses every finished task (step 6), then shows what it would delete and
   asks you to type the run name:

   ```
   Delete     the trajectories of 2571 completed, analysed task(s): 9.8 GB
   Keep       analysis/ (180 MB), trajectories.sha256, logs/, slurm/, tasks.csv, physics.yaml, run.yaml
              the files of unfinished tasks (1.3 GB): needed to resume them
   ```

Only trajectories of tasks that are **finished and analysed** are deleted (plus any
checkpoint a finished task left behind). Unfinished tasks keep their trajectory and
checkpoint, so `./noodle launch` can still resume them; run `./noodle clean` again
once they finish. Everything small stays, so `./noodle status` and `./noodle analyze`
keep working.

Before deleting, it writes `runs/<name>/trajectories.sha256`, the checksum of every
trajectory it deletes. It is the record of exactly which files went to the NAS, and
lets you check a copy at any time: put it in the copy's folder and run
`sha256sum -c trajectories.sha256` (each file should say `OK`).

To compute something new from the trajectories later, copy them back from the NAS into
`runs/<name>/trajectories/`.

---

## 3. Bulk, grafted, charged wall, and batches

The geometry is just the physics set you name in `run.yaml`:

```yaml
physics: baseline                              # free chains
physics: grafted                               # grafted by the N-terminus to a neutral surface
physics: grafted_reversed                      # grafted by the C-terminus to a neutral surface
physics: grafted_charged                       # grafted to a surface with -0.1 e/nm^2
physics: [baseline, grafted, grafted_charged]  # all three, for every protein, in one run
```

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

Every protein then runs against every surface. `./noodle analyze` gives each set its own
rows in `analysis/summary.csv` (column `physics`) and its own folder under
`analysis/timeseries/`.

### Several runs at once

You can have as many runs as you like going at the same time, of any geometry. Each run
lives in its own folder, its tasks are separate SLURM jobs named after the run, and
every task writes only into its own run folder, so runs never interfere. A large set can
be split into several runs by `rows` (e.g. `1-1000`, `1001-2000`, ...) to launch and
check it in pieces.

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

**Any geometry**

| setting | meaning | default |
|---|---|---|
| `geometry` | `bulk` or `grafted` | `bulk` |
| `temperature_K` | temperature, K | `298` |
| `ionic_concentration_M` | salt, mol/L (screens charges) | `0.150` |
| `epsilon_lj` | strength of the attraction between residues, kcal/mol | `0.1` |
| `lj_r_cut`, `yukawa_r_cut` | range of the residue attraction and of the electrostatics, nm | `3.5`, `3.5` |
| `harmonic_k`, `harmonic_r0` | bond stiffness (kcal/mol/nm²) and length (nm) | `1920`, `0.38` |
| `box_size` | box size, nm (grafted: its width) | `500` |
| `dt` | time step, in units of 0.489 ps | `0.01` |
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

**Run length.** Steps = multiplier x N^2.2 / dt (N = number of residues), because a
chain's relaxation time grows like N^2.2. The defaults are a rule of thumb from bulk
chains; grafted chains may need others. The `test_bulk` / `test_grafted` /
`test_grafted_charged` sets use tiny fixed lengths, for trying things out only.

---

## 5. Updating, and where things are

**Updating the code:** `git pull`. Your runs (`runs/`) and the environment (`env/`) are
not part of the code, so `git pull` never touches them. If `src/environment.yml`
changed, rebuild: `rm -rf env`, then section 1 again.

**A run folder**

```
runs/<name>/
  run.yaml, <your>.csv   your settings (and your own sequences, if any)
  tasks.csv              one row per task: protein, sequence, physics set, replicate, seed
  physics.yaml           the physics sets used, frozen at the first launch
  logs/                  each task's progress (what ./noodle status reads)
  trajectories/          the saved frames, one .gsd file per task (until ./noodle clean)
  checkpoints/           hourly save points of unfinished tasks
  slurm/                 each cluster task's printed output; submissions.csv
  analysis/              the ./noodle analyze tables
  trajectories.sha256    checksums of the trajectories ./noodle clean deleted
  status.txt             the last ./noodle status report
```

**The code** (start reading at `src/forcefield.py`, whose header explains how HOOMD
represents forces)

```
noodle                  the command: starts the matching program below
src/forcefield.py       the model: residue table, electrostatics, every force
src/simulate.py         one task: build the chain, equilibrate, record, checkpoint
src/launch_run.py       run.yaml + CSV -> tasks; submits them to SLURM or runs them here
src/new_run.py          ./noodle new        src/run_template.yaml   the run.yaml it copies
src/status.py           ./noodle status     src/analyze_run.py      ./noodle analyze
src/clean_run.py        ./noodle clean
src/runs.py             shared helpers      src/slurm_job.sh        the SLURM job script
src/check_physics.py    ./noodle check: the physics self-test (setup.sh runs it)
src/environment.yml     the packages setup.sh installs into env/
```

**GPU workstation:** change `cpu*` to `gpu*` in `src/environment.yml` before `setup.sh`.

---

## 6. Troubleshooting

| problem | fix |
|---|---|
| `The environment is missing` | setup hasn't been done on this computer, or failed (see `setup.log`) |
| `Permission denied` when running `./noodle` | `chmod +x noodle` |
| `There is no run called ...` | check the spelling; `./noodle status` lists the runs |
| `Can't find the sequence file` / `has no column` | the CSV must be in the run folder or in `sequences/`, and its file and column names must match `run.yaml` exactly (the column names of the shared sets are in section 2, step 2) |
| `physics.yaml has no set called ...` | typo in `physics:`; the message lists the sets |
| `[!] above ~1 kT/e` when launching a charged set | the surface charge is too high for the linearised model at this salt: see section 3 |
| changed `run.yaml` but nothing happened | after the first launch, only `time_limit_days`, `max_running`, `extra_sbatch` count; make a new run |
| delete a run | cancel its jobs first (`squeue -u $USER`, `scancel <job id>`), then `rm -rf runs/<name>` |

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
  at its contact distance from it. The surface has a repulsive core, an optional
  attraction scaled like hydropathy, and an optional charge, screened like the chain
  (linearised Poisson-Boltzmann, U(h) = 4π l_B kT q σ / κ · e^(−κh) for a residue of
  charge q at height h). A repulsive ceiling, out of the chain's reach, stops beads
  wrapping through the box.

`./noodle check` (also run by `setup.sh`) checks these against their formulas and runs
short bulk, grafted and charged-wall simulations. Every formula is explained in
`src/forcefield.py`.
