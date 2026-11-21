#!/bin/bash
#SBATCH --partition=notchpeak-shared-short
#SBATCH --account=notchpeak-shared-short
#SBATCH --time=8:00:00
#SBATCH --ntasks=4
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
# Step 1: MacroModel conformational search for every .sdf in this folder.
# Replaces prep_MM_npsh.bash. Settings come from mm_reference.com (unchanged from the ArBr library).
#
# Run from the folder that holds the .sdf files (e.g. ArCl/):
#     sbatch $PIPE_DIR/01_conformer_search.sbatch
# Safe to resubmit: molecules that already have an output file are skipped.
#
# Changes from prep_MM_npsh.bash:
#   - runs at most MAX_PARALLEL bmin jobs at once using "bmin -WAIT" (the old throttle used
#     curly quotes in grep, so it never matched and never limited anything)
#   - no hard-coded Schrödinger path inside the script (see config.sh), no "sleep 10000" at the end
#   - failed molecules are listed in mm_failed.txt instead of being silently missing

set -uo pipefail
PIPE_DIR="${PIPE_DIR:-$HOME/bin/arx_pipeline}"
WORK="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$WORK"
source "$PIPE_DIR/config.sh"
load_schrodinger

REF="$PIPE_DIR/mm_reference.com"
[[ -f "$REF" ]] || { echo "ERROR: $REF not found"; exit 1; }
shopt -s nullglob
sdfs=( *.sdf )
(( ${#sdfs[@]} )) || { echo "ERROR: no .sdf files in $WORK"; exit 1; }

mkdir -p output/clustered
: > mm_failed.txt

# Quick check that Schrödinger starts (licence and install) before preparing every molecule.
# Without this, a licence problem looks like a silent hang (SCHRODINGER_LICENSE_RETRY keeps retrying).
echo "Checking Schrödinger with ${sdfs[0]} ..."
if ! timeout 600 "$SCHRODINGER/utilities/structconvert" "${sdfs[0]}" schrodinger_check.mae > schrodinger_check.txt 2>&1; then
    echo "ERROR: structconvert did not finish within 10 min or failed. Messages:"; cat schrodinger_check.txt
    echo "Likely a licence or install problem with this version ($SCHRODINGER). Try another version in config.sh."
    exit 1
fi
rm -f schrodinger_check.mae schrodinger_check.txt
echo "Schrödinger OK"

todo=()
for sdf in "${sdfs[@]}"; do
    name="${sdf%.sdf}"
    if [[ -s "$name/${name}_out.maegz" ]]; then continue; fi       # already done
    mkdir -p "$name"
    "$SCHRODINGER/utilities/structconvert" "$sdf" "$name/$name.mae" > "$name/structconvert.txt" 2>&1 \
        || { echo "$name structconvert_failed (see $name/structconvert.txt)" >> mm_failed.txt; continue; }
    { echo "$name.mae"; echo "${name}_out.maegz"; cat "$REF"; } > "$name/$name.com"
    todo+=( "$name" )
done
echo "${#sdfs[@]} molecules, ${#todo[@]} to run, up to $MAX_PARALLEL at a time"

run_one() {
    local name="$1"
    cd "$WORK/$name" || return 0
    "$SCHRODINGER/bmin" -WAIT "$name" > bmin_submit.txt 2>&1
    if [[ -s "${name}_out.maegz" ]]; then echo "done   $name"
    else echo "FAILED $name"; echo "$name bmin_no_output" >> "$WORK/mm_failed.txt"; fi
}
export -f run_one
export WORK SCHRODINGER SCHRODINGER_LICENSE_RETRY

printf '%s\n' "${todo[@]}" | xargs -r -P "$MAX_PARALLEL" -I{} bash -c 'run_one "$1"' _ {}

n_fail=$(wc -l < mm_failed.txt)
echo "Finished. Failed: $n_fail (see mm_failed.txt). Next: bash $PIPE_DIR/02_cluster_conformers.sh"
