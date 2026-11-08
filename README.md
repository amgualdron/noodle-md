# noodle-md

Simulations of intrinsically disordered proteins (IDPs). Each amino acid is one bead,
and the beads form a chain (the HPS model). A chain is simulated either **free in
solution** (`bulk`) or **attached by one end to a surface** (`grafted`), for hundreds of
proteins at once on the cluster or a few on your own computer.

You use it with four commands and one settings file. No programming needed.

```bash
./new-run <name>        # make a run folder
./launch-run <name>     # start it (and resume anything that stopped)
./status <name>         # see how it's going
./analyze-run <name>    # turn finished simulations into CSV tables
```

Every command in this README is typed in a terminal **inside the `noodle-md` folder**.

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
./new-run grafted_tierA_neutral
```

Pick a name that says what it is (letters, digits, `_ - .`). This creates
`runs/grafted_tierA_neutral/` with a settings file, `run.yaml`.

### Step 2: add your sequences

Put a CSV file in the run folder. It needs a column with a **unique name** per protein
and a column with the **one-letter sequence**; other columns are ignored. Rows with
anything other than the 20 standard amino acids are skipped, with a message.

```bash
cp my_sequences.csv runs/grafted_tierA_neutral/
# from your own computer to the cluster:
scp my_sequences.csv your_username@cluster.address:noodle-md/runs/grafted_tierA_neutral/
```

`examples/` has two ready-made CSV files in this format.

### Step 3: fill in `run.yaml`

```bash
nano runs/grafted_tierA_neutral/run.yaml     # arrows to move, Ctrl-O Enter to save, Ctrl-X to quit
```

| setting | meaning | example |
|---|---|---|
| `description` | one line about the run, shown by `./status` | `"Tier A, neutral wall"` |
| `sequences_csv` | your CSV's file name (it must be in the run folder) | `my_sequences.csv` |
| `name_column`, `sequence_column` | the CSV columns with the names and the sequences | `ID`, `Sequence` |
| `rows` | `all`, one row (`7`) or a range (`101-200`); row 1 is the first protein | `1-50` |
| `physics` | one physics set from `physics.yaml`, or a list of several | `grafted` or `[baseline, grafted]` |
| `replicates` | independent copies (different random seeds) of each protein | `3` |
| `time_limit_days` | cluster time limit per task; a task that hits it can be resumed (step 5) | `3` |
| `max_running` | most tasks running at once (0 = no limit) | `0` |
| `extra_sbatch` | anything else your cluster needs | `["--partition=normal"]` |

The run will have *proteins x physics sets x replicates* tasks. Each task is one
simulation.

### Step 4: launch

```bash
./launch-run grafted_tierA_neutral
```

First it checks your settings (a misspelled setting is reported with a suggestion,
and nothing is launched). Then it shows the plan and asks before doing anything:

```
New run: grafted_tierA_neutral   "Tier A, neutral wall"

Sequences  my_sequences.csv, rows 1-40: 40 proteins, 23 to 497 residues (median 118)
Physics
  grafted          grafted by the N-terminus to a repulsive, neutral lj93 surface; 298 K, 150 mM salt
Tasks      40 proteins x 1 physics set(s) x 3 replicate(s) = 120 tasks
Length     steps per task (equilibration + production), for the shortest and longest protein
  grafted          2.38 million (N=23) to 2.05 billion (N=497), saving 3000 frames
Time       under a minute to ~6d 04h per task, ~81d 03h of computing in total
           [!] 33 task(s) need longer than the 3-day time limit: they will stop and need ./launch-run again ...
Disk       ~750 MB of trajectories
Cluster    120 SLURM jobs, time limit 3 day(s) each, no limit on how many at a time

Launch these 120 tasks? [y/N]
```

The time is a rough estimate from a 3-second test simulation. Answer `y` to launch;
anything else changes nothing. Then:

- **Cluster:** every task is submitted to SLURM as its own job, and the command checks
  that they arrived in the queue. You can log out.
- **Own computer:** the tasks run one after another in the terminal. Ctrl-C stops them;
  `./launch-run` again continues.

Launching **freezes** the run: the task list and a copy of the physics settings are
saved in the run folder. After that only `time_limit_days`, `max_running` and
`extra_sbatch` can still be changed; for anything else, make a new run.

> **Try it first.** For a 1-minute test, make a run with
> `cp examples/benchmark_sequences.csv runs/<name>/`, and in `run.yaml`
> `sequences_csv: benchmark_sequences.csv`, `rows: 1-2`,
> `physics: [test_bulk, test_grafted]`.

### Step 5: check on it, and resume what stopped

```bash
./status grafted_tierA_neutral
```

```
112 tasks  █████████████████░░░  96 done · 9 running · 4 waiting · 3 stopped

NOT FINISHED
  task   10  Q8XYZ1   N=412  grafted  stopped  stopped at 82% after 2d 23h, needs ~14h 20m more: ./launch-run resumes it
  task    9  P1234    N=288  grafted  running  production 33%, ~11h 06m left
```

| state | meaning |
|---|---|
| `waiting` | in the SLURM queue |
| `running` | running now; progress is updated every hour |
| `stopped` | not finished and not running (time limit, cluster problem, Ctrl-C) |
| `failed` | ended with an error, shown in the report (full message: its file in `runs/<name>/slurm/`) |
| `done` | finished |

**Stopped tasks:** run `./launch-run <name>` again. Each task saves a checkpoint every
hour, so it continues from there, losing at most an hour. It first tells you what it
will do:

```
120 tasks: 96 done · 10 in the queue (left alone) · 14 to start (12 resume from a checkpoint, 2 from the beginning)
```

Finished tasks and tasks still in the queue are left alone, so relaunching is always
safe. (A resumed task is just as valid as an uninterrupted one, but not frame-for-frame
identical.)

**Failed tasks** would most likely fail again, so `./launch-run` lists their errors and
asks before starting them again. Fix the cause first.

`./status` without a name lists all runs.

### Step 6: results

```bash
./analyze-run grafted_tierA_neutral
```

This writes `runs/<name>/timeseries/<physics set>/<protein>_pooled_timeseries.csv`, one
row per saved frame of every finished replicate. You can run it while the run is still
going.

| column | meaning (lengths in nm) |
|---|---|
| `sequence_name`, `physics`, `task_id`, `seed`, `frame` | where the row comes from |
| `rg`, `rg_sq` | radius of gyration (chain size), and its square |
| `ree_x`, `ree_y`, `ree_z`, `ree_sq` | end-to-end vector (first to last residue), and its length squared. For grafted chains z points away from the surface. |
| `asphericity` | shape: 0 = ball, 1 = rod |

The trajectories (`runs/<name>/trajectories/*.gsd`) open in [OVITO](https://www.ovito.org)
or in Python with the `gsd` package.

---

## 3. Bulk, grafted, and several runs at once

Bulk or grafted is just the physics set you name in `run.yaml`:

```yaml
physics: baseline              # free chains
physics: grafted               # chains grafted by the N-terminus to a neutral surface
physics: [baseline, grafted]   # both, for every protein, in the same run
```

You can have as many runs as you like going at the same time, bulk and grafted mixed.
Each run lives in its own folder, its tasks are separate SLURM jobs named after the
run, and every task writes only into its own run folder, so runs never interfere.

---

## 4. Changing the physics

All physics settings are named sets in **`physics.yaml`**. To try something new,
**don't edit the existing sets** (other people's new runs use them). Add a set at the
end that `extends` an existing one and lists only what changes:

```yaml
grafted_sticky_negative:
  extends: grafted
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
| `wall_charge_density` | surface charge, e/nm². Negative attracts K and R, repels D and E | `0.0` |
| `wall_headroom` | space between the fully stretched chain and the ceiling of the box, nm | `2.0` |

**Run length.** Steps = multiplier x N^2.2 / dt (N = number of residues), because a
chain's relaxation time grows like N^2.2. The defaults are a rule of thumb from bulk
chains; grafted chains may need others. The `test_bulk` / `test_grafted` sets use tiny
fixed lengths, for trying things out only.

---

## 5. Updating, and where things are

**Updating the code:** `git pull`. Your runs (`runs/`) and the environment (`env/`) are
not part of the code, so `git pull` never touches them. If `environment.yml` changed,
rebuild: `rm -rf env`, then section 1 again.

**A run folder**

```
runs/<name>/
  run.yaml, <your>.csv   your settings and sequences
  tasks.csv              one row per task: protein, physics set, replicate, seed
  physics.yaml           the physics sets used, frozen at the first launch
  logs/                  each task's progress (what ./status reads)
  trajectories/          the saved frames, one .gsd file per task
  checkpoints/           hourly save points of unfinished tasks
  slurm/                 each cluster task's printed output; submissions.csv
  timeseries/            the ./analyze-run tables
  status.txt             the last ./status report
```

**The code** (start reading at `src/forcefield.py`, whose header explains how HOOMD
represents forces)

```
src/forcefield.py      the model: residue table, electrostatics, every force
src/simulate.py        one task: build the chain, equilibrate, record, checkpoint
src/launch_run.py      run.yaml + CSV -> tasks; submits them to SLURM or runs them here
src/status.py          ./status          src/analyze_run.py   ./analyze-run
src/runs.py            shared helpers    src/slurm_job.sh     the SLURM job script
src/check_physics.py   physics self-test (env/bin/python src/check_physics.py)
```

**GPU workstation:** change `cpu*` to `gpu*` in `environment.yml` before `setup.sh`.

---

## 6. Troubleshooting

| problem | fix |
|---|---|
| `The environment is missing` | setup hasn't been done on this computer, or failed (see `setup.log`) |
| `Permission denied` when running a command | `chmod +x new-run launch-run status analyze-run` |
| `There is no run called ...` | check the spelling; `./status` lists the runs |
| `Can't find the sequence file` / `has no column` | the CSV must be in the run folder, and its file and column names must match `run.yaml` exactly |
| `physics.yaml has no set called ...` | typo in `physics:`; the message lists the sets |
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
  (linearised Poisson-Boltzmann). A repulsive ceiling, out of the chain's reach, stops
  beads wrapping through the box.

`src/check_physics.py`, which `setup.sh` runs, checks these against their formulas and
runs short bulk and grafted simulations. Every formula is explained in
`src/forcefield.py`.
