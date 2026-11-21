# ArX DFT pipeline (same workflow as the ArBr library)

Takes one SDF per molecule to finished Gaussian logs with the settings used for the ArBr library: MacroModel conformer search → Kelley clustering → Gaussian 16 opt+freq and single points → log checks and restarts. Used for the ArCl and ArI analogues, and for ArF, ArOTs and ArOTf.

An updated version of Haas and Hardy's *automated_workflow_for_SMILES_to_logs* (v1.1.0). Descriptors are then collected with SigmanGroup Get_Properties v1.0.3 (Zenodo, DOI 10.5281/zenodo.10651727) plus `lp_align.py`. This version contains automated clustering when clustering fails.

## Contents

| File | Step | What it does |
|---|---|---|
| `config.sh` | all | Schrödinger and Python modules, `MAX_CONF`, `N_CLUSTERS`, `MAX_PARALLEL` |
| `mm_reference.com` | 1 | MacroModel search settings (OPLS4, gas phase, 21 kJ/mol window) |
| `01_conformer_search.sbatch` | 1 | MacroModel conformer search for every SDF in the folder |
| `02_cluster_conformers.sh` + `kelley_schrodinger.py` | 2 | Keep ≤ 20 conformers, or RMSD-cluster at the minimum Kelley penalty |
| `method_config.py` | 3–5 | Gaussian routes, cores/memory, iodine basis, expected terminations |
| `03_make_gaussian_inputs.py` | 3 | Writes one linked-job `.com` per conformer into `gaussian/` |
| `04_submit_gaussian.sbatch` | 4 | SLURM job array, one conformer per task; skips finished logs |
| `05_check_logs.py` | 5 | Sorts logs into `logs/`, `term_error/`, `imag/`, `resubmit/` |
| `06_restart_from_log.py` | 6 | Rebuilds inputs for failed or imaginary-frequency conformers |
| `check_progress.sh` | any | Read-only progress count; safe while jobs are running |

## Methods (as in the ArBr library SI)

| ArBr SI | Setting here | File |
|---|---|---|
| MacroModel conformational search, OPLS4, gas phase | `FFLD 16`, no solvent line | `mm_reference.com` |
| Convergence 0.001, max 10,000 iterations | `CONV … 0.0010`, `MINI 1 0 10000` | `mm_reference.com` |
| Conformers within 5.02 kcal/mol | `MCSS`/`DEMX` 21 kJ/mol | `mm_reference.com` |
| > 20 conformers → RMSD clustering at the minimum Kelley penalty, centroids kept | `MAX_CONF=20`, `N_CLUSTERS=0`, `conformer_cluster.py -n 0 -r -rep` | `config.sh`, step 2 |
| Gaussian 16 C.01 | `module load gaussian16` | step 4 |
| Opt + freq B3LYP-D3(BJ)/6-31G(d,p), ultrafine grid | link 1 | `method_config.py` |
| SP M06-2X/def2-TZVP, NBO7, polarizability, EFG, volume | link 2 | `method_config.py` |
| SP M06-2X/def2-TZVP, NMR (GIAO) | link 3 | `method_config.py` |

| Link | Route |
|---|---|
| 1 | `#N empiricaldispersion=GD3BJ OPT FREQ=NORAMAN Integral(Ultrafine) B3LYP 6-31G(d,p)` |
| 2 | `#N int=(grid=ultrafine) pop=nbo7 polar prop=efg volume guess=read geom=check M062X def2TZVP` |
| 3 | `#N int=(grid=ultrafine) nmr=giao guess=read geom=check M062X def2TZVP` |

`%NProcShared=16`, `%Mem=16GB`, as in the ArBr inputs. Schrödinger 2024-4 is used (2021-4 can no longer check out a licence on CHPC); state the release in the SI.

**Iodine.** 6-31G(d,p) is defined only for H–Kr, so ArI inputs use GenECP: 6-31G(d,p) on light atoms with def2-SVP + def2 ECP on I for the optimization, and def2-TZVP + def2 ECP for the single points (`HEAVY` in `method_config.py`). Report this in the SI. Cl, F and S need no change.

**OTs / OTf.** The leaving group is attached through O, so the descriptor "site" for these is the ester O (see `lp_align.py`). The extra C–O–S–C torsions give more conformers, so more molecules reach clustering and DFT takes longer.

## One-time setup on CHPC

```bash
mkdir -p ~/bin && cd ~/bin
git clone git@github.com:<your-github-username>/<this-repo>.git arx_pipeline   # or: unzip ~/arx_pipeline.zip
chmod 744 ~/bin/arx_pipeline/*
echo 'export PIPE_DIR=~/bin/arx_pipeline' >> ~/.bashrc
source ~/.bashrc
```

Step 3 needs Python with RDKit (e.g. a conda environment). If a CHPC Python module is loaded at login, it hides the environment's Python:

```bash
conda activate <env with rdkit>
which python                     # must end in .../envs/<env>/bin/python
module unload python; hash -r    # if it shows /uufs/.../installdir/python/...
python -c "import rdkit; print(rdkit.__version__)"
```

If step 3 fails with a NumPy `_ARRAY_API` error, run `export PYTHONNOUSERSITE=1` first. Use `$CONDA_PREFIX/bin/python` if in doubt.

## Running it

One working folder per leaving group, each holding only that group's SDFs (from the selection notebooks):

```bash
mkdir -p ~/arx_r2/ArCl && cp sdf/Cl/*.sdf ~/arx_r2/ArCl/
cd ~/arx_r2/ArCl
```

### Step 1: conformer search

```bash
sbatch $PIPE_DIR/01_conformer_search.sbatch
squeue -u $USER                                # gone from the list = finished
ls */*_out.maegz 2>/dev/null | wc -l           # molecules done
cat mm_failed.txt                              # should be empty
```

The job runs on `notchpeak-shared-short` (8 h limit). If `sacct` shows `TIMEOUT`, submit it again: finished molecules are skipped. Several conformer searches at once share MacroModel licence tokens, so one set at a time is often faster overall. If the first lines of `slurm-*.out-*` show structconvert didn't finish within 10 min, it's a licence problem (`SCHRODINGER_MODULE` in `config.sh`).

### Step 2: keep or cluster conformers

```bash
sbatch $PIPE_DIR/02_cluster_conformers.sh
cat cluster_failed.txt                          # should be empty
column -s, -t conformer_summary.csv | less      # k_default vs k_used per molecule
```

Kelley penalty plots: `kelley/<name>_kelley.png`.

### Step 3: Gaussian inputs

```bash
python $PIPE_DIR/03_make_gaussian_inputs.py
cat gaussian/problems.csv                       # header only = no problems
ls gaussian/*.com | wc -l
```

### Step 4: submit

```bash
cd gaussian
ls *.com > jobs.txt
sbatch --array=1-$(wc -l < jobs.txt)%50 $PIPE_DIR/04_submit_gaussian.sbatch
```

- `%50` limits how many run **at the same time**; the rest wait and start as others finish. Change it on a queued array with `scontrol update JobId=<id> ArrayTaskThrottle=200`.
- `(JobArrayTaskLimit)` in `squeue` = waiting behind that limit; `(Resources)` = no free nodes.
- A preempted task is requeued; finished logs are never rerun.
- A different list file can be passed as the first argument (used below for chunks and restarts).

**Other clusters.** Command-line options override the `#SBATCH` lines, and `-M` picks the cluster (check names with `myallocation`):

```bash
sbatch -M granite --partition=granite-guest --qos=granite-guest --account=sigman \
       --array=1-$(wc -l < jobs.txt)%50 $PIPE_DIR/04_submit_gaussian.sbatch
squeue -M granite -u $USER                      # jobs on other clusters need -M
```

Never queue the same folder on two clusters at once: both arrays would run the unfinished inputs and overwrite each other's logs. Test one task (`--array=1`) on a new cluster first, in case its Gaussian module has another name.

**Large folders.** Array indices can't exceed `MaxArraySize`, and queued + running tasks count against a per-user limit (`QOSMaxSubmitJobPerUserLimit`):

```bash
scontrol show config | grep -i MaxArraySize
sacctmgr show qos notchpeak-guest format=Name,MaxSubmitPU
squeue -u $USER -r -h | wc -l                   # what you already have queued + running
split -l 400 -d jobs.txt todo_
sbatch --array=1-$(wc -l < todo_00)%200 $PIPE_DIR/04_submit_gaussian.sbatch todo_00
```

Submit `todo_01`, `todo_02` … as the queue drains. To list inputs that never started (only once the folder has nothing left in the queue, since queued tasks have no log yet either):

```bash
for f in *.com; do [[ -e "${f%.com}.log" ]] || echo "$f"; done > not_started.txt
```

### Progress (any time)

```bash
bash $PIPE_DIR/check_progress.sh                       # the gaussian/ folder you're in
bash $PIPE_DIR/check_progress.sh -v ~/arx_r2/*/gaussian  # several folders; -v lists failed jobs and error links
```

Counts finished (4 normal terminations), failed, in progress and not started. Read-only. SLURM may label a task `FAILED` even when its log finished normally; the log is what counts.

### Step 5: sort the logs (after the folder's jobs have left the queue)

```bash
python $PIPE_DIR/05_check_logs.py
ls logs | wc -l
ls term_error imag resubmit
grep ",False$" molecule_status.csv               # molecules not yet complete
```

Don't run step 5 while that folder still has running tasks: a log still being written looks like a failure.

### Step 6: fix failures

| Where | Cause | Fix (in `gaussian/`) |
|---|---|---|
| `resubmit/` | input never ran or was killed | `(cd resubmit && ls *.com) > jobs.txt; mv resubmit/*.com .` then submit |
| `term_error/`, l9999 | optimization ran out of steps | `--all-l9999` (last geometry, `Opt=(MaxCycles=200,CalcFC)`) |
| `term_error/`, "oscillating" force | flat torsion | `--opt "(MaxCycles=300,CalcAll)" <log>`; if still stuck `--opt "(MaxCycles=300,CalcFC,MaxStep=5)"` |
| `term_error/`, l103 / FormBX | internal coordinates broke | `--opt "(MaxCycles=300,CalcFC,Cartesian)" <log>` |
| `term_error/`, no error line | killed (memory, time, preemption) | check `sacct ... State,MaxRSS`; resubmit, with `--mem=40G` if `OUT_OF_MEMORY` |
| `imag/` | saddle point | `--displace` (below) |

Each restart writes new inputs into the current folder and lists them in `jobs_restart.txt`:

```bash
python $PIPE_DIR/06_restart_from_log.py --all-l9999
cp jobs_restart.txt jobs.txt
sbatch --array=1-$(wc -l < jobs.txt) $PIPE_DIR/04_submit_gaussian.sbatch
```

Submit each restart batch before running the next restart command (each one rewrites `jobs_restart.txt`), and check that `jobs_restart.txt` exists before copying it. Then rerun step 5.

**Imaginary frequencies.** `--displace` moves the optimized geometry along the imaginary mode (largest atom displacement `--amp`, default 0.15 Å) and re-optimizes; `--minus` pushes the other way. `--amp` is a script option, not a Gaussian keyword; Gaussian receives ordinary Cartesian coordinates.

| Imaginary frequency | Likely cause | Try |
|---|---|---|
| 0 to about −15 cm⁻¹ | flat rotor (e.g. CH₃, CF₃) | `--displace --opt "(MaxCycles=300,CalcFC,Tight)"`; if it persists, keep or drop by a stated rule and report it |
| about −20 to −40 | shallow torsional saddle | `--displace --amp 0.3` |
| below about −70 | true saddle point | `--displace --amp 0.3`, then `--minus`; drop the conformer if the molecule has others |

`Tight` tightens only the optimizer's stopping criteria (same level of theory); report the conformers that needed it. Don't use `Opt=Loose`: geometries must stay comparable with the ArBr library.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError: rdkit` although `conda list rdkit` shows it | a CHPC Python module is ahead of the environment | `module unload python; hash -r`, or `$CONDA_PREFIX/bin/python` |
| `Invalid job array specification` | `jobs.txt` empty, or more lines than `MaxArraySize` | check `wc -l < jobs.txt`; split into chunks |
| `QOSMaxSubmitJobPerUserLimit` | too many tasks queued | submit chunks as the queue drains, or use another cluster |
| step 1 `TIMEOUT` after 8 h | partition time limit | resubmit (finished molecules skipped) |
| kernel crash in a selection notebook | out of memory (13k × 13k similarity matrix) | close other kernels; more memory for the session |

Useful SLURM commands: `squeue -u $USER`, `squeue -u $USER -r` (every array task), `sacct -u $USER -S today --format=JobID,JobName%25,State,Elapsed,MaxRSS`, `scancel <jobid>`, `scancel -u $USER -t PENDING`. Cancel only one project's pending jobs by name: `squeue -u $USER -h -o "%i %j" | awk '$2 ~ /^arots/ {print $1}' | xargs scancel`.

## Single-cluster results (automated in step 2)

With `-n 0`, `conformer_cluster.py` picks the number of clusters with the lowest Kelley penalty, which can be 1. When the default run returns one representative, step 2 reads the penalty table (`<name>/<name>_cluster_ligand1.txt`, 6th column), takes the lowest penalty at **2 or more** clusters (the "second lowest" rule of the old manual workflow), moves the first run's files to `<name>/default_n0/`, and reruns with `-n <that number>`. Molecules whose default is ≥ 2 clusters are unchanged, as in the ArBr library. Both numbers are in `conformer_summary.csv` (`k_default`, `k_used`).

## Changes from the original scripts

| Original | Problem | Now |
|---|---|---|
| `SMILES_to_sdf.py` | bonds rebuilt from xyz by Open Babel (charges can be wrong), random embedding | SDFs written by RDKit in the selection notebooks, fixed seed, stereo shared across leaving groups |
| `prep_MM_npsh.bash` | throttle never matched (curly quotes); hard-coded 2021-4 | step 1: `bmin -WAIT` with `MAX_PARALLEL`, skips finished molecules, version in `config.sh` |
| `MM_clustering_and_converting.bash` | conformer count from digits in log lines | step 2 counts structures; same clustering command |
| manual single-cluster fix (several scripts) | by hand, deletes files | automatic in step 2, penalty curve plotted, first run kept |
| `check_manual_clustering.bash` + `g16_inpgen_arbrlib.py` | calls `babel` (gone in Open Babel 3); 6-31G(d,p) on iodine | step 3: complete linked inputs, GenECP for I, atom/charge checks |
| post-processing scripts | expected terminations edited by hand; fixed-column parsing | step 5: count from `method_config.py`, summary CSVs, safe to rerun |
| — | no submission, restart or progress scripts | steps 4 and 6, `check_progress.sh` |

## Notes

- File names follow the ArBr workflow: `arcl3949_conf-1.log`, or `arcl3949_clust-1.log` for clustered molecules.
- Followed ArBr library workflow (see SI) and Get_Properties (Zenodo, DOI 10.5281/zenodo.10651727) for the descriptor collection.
