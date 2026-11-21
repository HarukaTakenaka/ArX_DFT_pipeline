#!/bin/bash
#SBATCH --partition=notchpeak-guest      # owner-guest (preemptable); check names with: myallocation
#SBATCH --account=owner-guest
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16                       # = %NProcShared in method_config.py
#SBATCH --mem=20G                         # a few GB above %Mem
#SBATCH --requeue                         # if preempted, go back in the queue and start again
#SBATCH -o slurm-%A_%a.out
#SBATCH -e slurm-%A_%a.err
#
# Step 4: run every .com in this folder as a SLURM job array, one conformer per task.
#
# Run from the gaussian/ folder made by step 3:
#     ls *.com > jobs.txt
#     sbatch --array=1-$(wc -l < jobs.txt)%50 $PIPE_DIR/04_submit_gaussian.sbatch
#     sbatch --array=1-$(wc -l < todo_00)%50 $PIPE_DIR/04_submit_gaussian.sbatch todo_00   # another list file
# (%50 = at most 50 jobs running at once.) To resubmit only failed jobs, write their names to
# jobs.txt and submit again with the new count.
#
# Owner-guest jobs can be stopped at any time when the node's owner needs it. With --requeue, SLURM
# restarts them from the beginning; a conformer whose log already has every expected "Normal termination"
# is skipped, so a requeued task never redoes finished work.

set -uo pipefail
module load gaussian16                    # set to the module from "module avail gaussian" (ArBr: G16 C.01)

EXPECTED=4                                # "Normal termination" lines in a finished log (method_config.py)
LIST=${1:-jobs.txt}                       # list of inputs; default jobs.txt (pass another file as the first argument)
com=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
name="${com%.com}"

if [[ -f "$name.log" ]] && (( $(grep -c "Normal termination" "$name.log") >= EXPECTED )); then
    echo "$name already finished; skipping"; exit 0
fi

export GAUSS_SCRDIR="/scratch/local/$USER/$SLURM_JOB_ID"   # node-local scratch
mkdir -p "$GAUSS_SCRDIR"
g16 < "$com" > "$name.log"
rm -rf "$GAUSS_SCRDIR"
