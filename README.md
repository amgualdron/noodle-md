# noodle-md

Computer simulations of disordered proteins. You give it a list of protein sequences;
it simulates each one, on the university cluster or on your own computer, and gives you
a table of how each protein behaves: how big it is, what shape it takes, how high it
stands when attached to a surface, and more.

Each protein can be simulated:

- **free** in salt water,
- **attached by one end to a surface** (neutral or charged),
- as a **brush**: 100 copies attached side by side.

You only type short commands and fill in one settings file. No programming needed.

| section | |
|---|---|
| [1. Set up](#1-set-up-once) | once per computer |
| [2. Your first run](#2-your-first-run) | from sequences to results, in 6 steps |
| [3. All commands](#3-all-commands) | the list |
| [4. Choosing the physics](#4-choosing-the-physics) | free, grafted, charged, brush |
| [5. Your results](#5-your-results) | which file has what |
| [6. Good to know](#6-good-to-know) | time limits, stopping, disk space |
| [7. If something goes wrong](#7-if-something-goes-wrong) | common messages and fixes |

**How it works in detail**, every formula and every file: [REFERENCE.md](REFERENCE.md).

---

## 1. Set up (once)

You type everything below in a **terminal** (on the cluster: after logging in with
`ssh your_username@cluster.address`). On your own computer you also need conda: if
`conda --version` says "command not found", install
[Miniforge](https://github.com/conda-forge/miniforge) first.

**1. Download the code**

```bash
git clone <URL of this repository> noodle-md
cd noodle-md
```

**2. Install it** (15 to 30 minutes)

| on the cluster | on your own computer |
|---|---|
| `sbatch setup.sh`, then check on it with `tail setup.log` | `bash setup.sh` |

It worked when the end says `All physics checks passed.` and `Setup finished`.

**3. Make `noodle` work from any folder.** The end of the setup prints two lines like
these. Copy them into the file `~/.bashrc` (open it with `nano ~/.bashrc`, paste at the
end, save with Ctrl-O Enter, leave with Ctrl-X), then open a new terminal:

```bash
alias noodle='/home/you/noodle-md/noodle'
source '/home/you/noodle-md/noodle-completion.bash'
```

Now `noodle` works anywhere, and the **Tab** key completes commands and run names.
(Without these lines, go into the `noodle-md` folder and type `./noodle` instead.)

---

## 2. Your first run

A **run** is one study: some proteins, some physics, in its own folder. Every run also
gets a **number**, which you can type instead of its name.

**Step 1. Make the run**

```bash
noodle new my_first_run
```

```
Created runs/my_first_run/, run number 1
```

**Step 2. Fill in its settings.** Open `runs/my_first_run/run.yaml` (for example with
`nano runs/my_first_run/run.yaml`). The main settings:

| setting | what to write | example |
|---|---|---|
| `description` | one line saying what the run is for | `"Top 500 human IDRs, charged wall"` |
| `sequences_csv` | the file with your sequences (see below) | `human_idrs_12k.csv` |
| `rows` | which proteins of that file: `all`, one (`7`), a range with a dash (`1-500`), or several in brackets (`[1-50, 75]`) | `1-500` |
| `physics` | how to simulate them (section 4); a list for several | `baseline` or `[baseline, grafted]` |
| `replicates` | how many independent copies of each simulation | `3` |
| `max_residues` | skip proteins longer than this (0 = no limit) | `200` |

Square brackets always mean "several": several rows, several physics sets.

The other settings are for the cluster and are explained in the file itself. One is
worth knowing: **`max_jobs`** (default 1000), the most cluster jobs a run is split
into. A small run gets one job per simulation, all running side by side. A large run
is split by protein length: the longest proteins get a job each, and the short ones
share jobs, so the cluster queue isn't flooded.

**Sequences.** Two sets come with the code, in the `sequences/` folder:
`human_idrs_12k.csv` (12,000 human disordered regions, 20-599 residues) and
`mobidb_idps_15k.csv` (15,000 disordered proteins, 20-500 residues). Both are ordered
so that the first rows are the most different from each other: `rows: 1-500` gives
500 very varied proteins. To use **your own**, put a CSV file in the run folder with a
column `id` (a unique name per protein) and a column `sequence` (one-letter code).

> **Try a test first:** `physics: [test_bulk, test_grafted]` and `rows: 1-2` runs in
> about a minute and shows you every step, with results that mean nothing.

**Step 3. Launch it**

```bash
noodle launch my_first_run
```

It checks your settings, shows what it will do (how many simulations, roughly how
long, how much disk space) and asks before starting. On the cluster you can then log
out. On your own computer the simulations run in that terminal.

**Step 4. Check on it**, as often as you like:

```bash
noodle status my_first_run
```

```
████████████████████████▓░░░░░  2571 of 3000 tasks done (86%)
  2 failed · 14 stopped · 9 running · 404 waiting · 2571 done
```

followed by what needs your attention, the time left, and what to do next.

**Step 5. Get the results**, also while the run is still going:

```bash
noodle analyze my_first_run
```

The results are in `runs/my_first_run/manifest.csv`: one row per protein (section 5).

**Step 6. Take them home.** `noodle export my_first_run` puts everything into one
folder, `exports/my_first_run/`. Copy it to your computer from there.

---

## 3. All commands

`<run>` is a run's **name**, its **number**, or nothing at all when you are inside the
run's folder.

| command | what it does |
|---|---|
| `noodle new <name>` | makes a new run folder with a settings file |
| `noodle launch <run>` | starts the run; later, continues whatever stopped |
| `noodle status <run>` | how the run is going |
| `noodle status` | every run, with its number |
| `noodle analyze <run>` | gathers the results into tables |
| `noodle export <run> <run> ...` | copies the results of one or more runs into `exports/` |
| `noodle cancel <run>` | stops the run's simulations on the cluster (`launch` continues them later) |
| `noodle clean <run>` | deletes the large simulation files, keeping the results (section 6) |
| `noodle check` | tests that everything works |
| `noodle help` | this list; `noodle help <command>` explains one command in detail |

---

## 4. Choosing the physics

Write the name of one of these in `physics:` in `run.yaml`, or a list of several in
brackets: every protein is then simulated under each.

| name | what is simulated |
|---|---|
| `baseline` | the protein free in salt water (room temperature, 150 mM salt) |
| `grafted` | attached by its first residue (N-terminus) to a neutral surface |
| `grafted_reversed` | attached by its last residue (C-terminus) |
| `grafted_charged` | attached to a negatively charged surface (−0.1 charges per nm²) |
| `brush` | 100 copies attached side by side, 0.05 per nm² |
| `brush_charged` | the same brush on the charged surface |
| `test_bulk`, `test_grafted`, `test_grafted_charged`, `test_brush` | tiny versions, to try things out |

**To change something** (temperature, salt, surface charge, brush density...), don't
edit these. Add your own set at the end of `physics.yaml`, starting from one of them and
listing only what changes:

```yaml
grafted_hot:
  extends: grafted
  description: "Grafted, at 320 K"
  temperature_K: 320
```

Then use `physics: grafted_hot` in your run. Every setting is listed, with its meaning,
in `physics.yaml` itself and in [REFERENCE.md](REFERENCE.md#3-the-model).

**Brushes are expensive**: a brush is 100 proteins at once, and the time grows quickly
with length (about 17 hours for 100 residues, a week for 200). Set `max_residues` for
brush runs. The launch shows the cost before you agree.

---

## 5. Your results

`noodle analyze` writes these files in the run folder:

| file | what is in it |
|---|---|
| **`manifest.csv`** | **the main result**: one row per protein and physics set |
| `analysis/replicates.csv` | the same, one row per single simulation, to see how much copies differ |
| `analysis/profiles.csv` | curves, e.g. how density changes with height above the surface |
| `analysis/report.json` | a summary, and which rows to be careful with |

All are ordinary tables: open them in Excel, LibreOffice or Python. Each row of
`manifest.csv` has the protein's name and sequence, the conditions, numbers describing
the sequence (charge, hydrophobicity, charge patterning), then the measured averages,
each with its uncertainty in the column ending in `_err`. The ones most people want
(lengths in nm):

| column | meaning |
|---|---|
| `rg_nm` | radius of gyration: the overall size of the protein |
| `ree_nm` | distance between its two ends |
| `rh_nm` | hydrodynamic radius (the size seen by light scattering, roughly) |
| `asphericity` | shape: 0 = round, 1 = rod-like |
| `nu` | how open the chain is: about 0.33 compact, 0.5 relaxed, 0.6 swollen |
| `com_height_nm` | attached proteins: average height above the surface |
| `adsorbed_fraction` | attached proteins: share of the protein touching the surface |
| `brush_height_nm` | brushes: thickness of the brush |
| `independent_samples` | quality check: **below 20, the simulation was too short to trust** |
| `drift` | quality check: should be 0 (otherwise still changing at the end) |

Every column is defined in [REFERENCE.md](REFERENCE.md#7-observables-and-statistics).
`noodle analyze` also warns you about rows with poor quality checks.

---

## 6. Good to know

- **Nothing is lost when a simulation stops.** On the cluster each job has a time limit.
  Shortly before it, a simulation saves its progress and stops. The same happens with
  `noodle cancel` or Ctrl-C. Run `noodle launch <run>` again, and it continues
  from where it stopped. Running `launch` again is always safe: finished and still
  queued simulations are left alone.
- **Failed simulations** are shown by `noodle status` with their error. Fix the cause
  first; `noodle launch` asks before trying them again.
- **A launched run is frozen.** Changing `run.yaml` or `physics.yaml` afterwards
  doesn't change it (except the cluster settings). For something different, make a
  new run.
- **Many runs at once** are fine; they never interfere.
- **Disk space.** The simulation files (`trajectories/`) are large. When a run is done
  and analysed, copy `runs/<name>/trajectories/` to the group storage (NAS), then
  `noodle clean <run>` deletes them from the cluster and keeps all results. It asks
  you to type the run's name first.
- **Updating the code:** `git pull` (your runs are never touched). See `CHANGELOG.md`
  for what changed.

---

## 7. If something goes wrong

| problem | what to do |
|---|---|
| `noodle: command not found` | the two lines aren't in `~/.bashrc` yet, or you need a new terminal (section 1). Meanwhile, use `./noodle` inside `noodle-md` |
| `The environment is missing` | the setup didn't finish: look at the end of `setup.log`, then run it again |
| setup: `Could not load the module 'anaconda'` | find your cluster's name for it with `module avail conda` and write it in `CONDA_MODULE=` at the top of `setup.sh` |
| setup: download errors | the cluster's compute nodes have no internet: run `bash setup.sh` instead of `sbatch` |
| setup: `./env already exists` | an earlier attempt was left half done: `rm -rf env`, then again |
| `There is no run called ...` | check the name or number: `noodle status` lists them |
| `Can't find the sequence file` / `has no column` | the CSV must be in the run folder or in `sequences/`, with the names written exactly as in `run.yaml` |
| `physics.yaml has no set called ...` | a typo in `physics:`; the message lists the right names |
| changed `run.yaml` but nothing happened | a launched run is frozen: make a new run |
| a task failed | `noodle status <run>` shows the error; the full message is in `runs/<name>/jobs/job_<number>.out` |
| delete a run | `noodle cancel <run>`, then `rm -rf runs/<name>` |
