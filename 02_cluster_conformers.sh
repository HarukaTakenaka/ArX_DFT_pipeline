#!/bin/bash
#SBATCH --partition=notchpeak-shared-short
#SBATCH --account=notchpeak-shared-short
#SBATCH --time=8:00:00
#SBATCH --ntasks=4
#SBATCH -o slurm-%j.out-%N
#SBATCH -e slurm-%j.err-%N
#
# Step 2: keep or cluster each conformer ensemble, as for the ArBr library.
#   <= MAX_CONF conformers -> all kept:            output/<name>_conf.sdf
#   >  MAX_CONF conformers -> RMSD clustering at the minimum Kelley penalty,
#                             cluster representatives: output/<name>_clust.sdf
#                             full ensemble kept in:   output/clustered/<name>_conf.sdf
# Same conformer_cluster.py options as MM_clustering_and_converting.bash.
#
# Single-cluster fix (automates the old manual "workflow for clustering if failed"):
#   conformer_cluster.py writes the Kelley penalty for every number of clusters, starting at 1, to
#   <name>_cluster_ligand1.txt. When its default choice (-n 0) is ONE cluster, this script takes the
#   lowest penalty at 2 or more clusters from that same file, moves the first run's files to
#   <name>/default_n0/ (the old workflow deleted them), and reruns conformer_cluster.py with -n <that number>.
#   Every clustered molecule also gets kelley/<name>_kelley.png (Schrödinger's penalty curve).
#
# Run from the same folder as step 1 (bash or sbatch both work):
#     bash $PIPE_DIR/02_cluster_conformers.sh
#
# Other changes from MM_clustering_and_converting.bash:
#   - the conformer count is taken by counting structures in <name>_out.maegz, not by stripping
#     every digit from two lines of the log (digits in names such as arcl3949 could end up in the count)
#   - a summary table (conformer_summary.csv) and a list of problems (cluster_failed.txt)

set -uo pipefail
PIPE_DIR="${PIPE_DIR:-$HOME/bin/arx_pipeline}"
WORK="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$WORK"
source "$PIPE_DIR/config.sh"
load_schrodinger
load_python

mkdir -p output/clustered kelley
: > cluster_failed.txt
echo "name,n_conformers,clustered,k_default,k_used,n_kept" > conformer_summary.csv

run_cluster() {   # $1 = name, $2 = -n value; run inside the molecule folder
    "$SCHRODINGER/run" conformer_cluster.py -j "${1}_cluster" -n "$2" -r -rep -WAIT "${1}_out.maegz" > cluster.txt 2>&1
    local rep="${1}_cluster_ligand1_representatives.maegz"
    rm -f "${1}_clust.sdf"
    [[ -s "$rep" ]] && "$SCHRODINGER/utilities/structconvert" "$rep" "${1}_clust.sdf" >/dev/null 2>&1
}
count() { grep -c '^\$\$\$\$' "$1" 2>/dev/null || echo 0; }

for dir in */; do
    name="${dir%/}"
    [[ -f "$name.sdf" ]] || continue                     # only molecule folders made by step 1
    out="$name/${name}_out.maegz"
    [[ -s "$out" ]] || { echo "$name no_conformer_search_output" >> cluster_failed.txt; continue; }

    ( cd "$name"
      "$SCHRODINGER/utilities/structconvert" "${name}_out.maegz" "${name}_conf.sdf" >/dev/null 2>&1 )
    n=$(count "$name/${name}_conf.sdf")
    if (( n == 0 )); then echo "$name empty_ensemble" >> cluster_failed.txt; continue; fi

    if (( n <= MAX_CONF )); then
        cp "$name/${name}_conf.sdf" "output/${name}_conf.sdf"
        echo "$name,$n,no,,,$n" >> conformer_summary.csv
        echo "$name: $n conformers, kept all"
        continue
    fi

    ( cd "$name"; rm -rf default_n0 "${name}"_cluster* "${name}_clust.sdf"; run_cluster "$name" "$N_CLUSTERS" )   # clear any earlier run
    k_default=$(count "$name/${name}_clust.sdf")
    k_used=$k_default
    table="$name/${name}_cluster_ligand1.txt"

    if [[ -s "$table" ]]; then
        k_best=$("$PYTHON" "$PIPE_DIR/kelley_schrodinger.py" "$table" \
                 --plot "kelley/${name}_kelley.png" --title "$name: $n conformers") || k_best=""
        cp "$table" "kelley/${name}_cluster_ligand1.txt"
    else
        k_best=""; echo "$name no_kelley_table" >> cluster_failed.txt
    fi

    if (( k_default == 1 )) && [[ "$N_CLUSTERS" == "0" ]]; then
        if [[ -n "$k_best" ]]; then
            ( cd "$name"
              mkdir -p default_n0
              mv "${name}"_cluster* "${name}_clust.sdf" cluster.txt default_n0/ 2>/dev/null
              run_cluster "$name" "$k_best" )
            k_used=$(count "$name/${name}_clust.sdf")
            echo "$name: default gave 1 cluster; reran with -n $k_best (lowest Kelley penalty at >= 2 clusters)"
            (( k_used == k_best )) || echo "$name rerun_gave_${k_used}_not_${k_best}" >> cluster_failed.txt
        else
            echo "$name clustered_to_one_no_table_cluster_manually" >> cluster_failed.txt
        fi
    fi

    if (( k_used > 0 )); then
        cp "$name/${name}_clust.sdf" "output/${name}_clust.sdf"
        cp "$name/${name}_conf.sdf"  "output/clustered/${name}_conf.sdf"
        echo "$name,$n,yes,$k_default,$k_used,$k_used" >> conformer_summary.csv
        echo "$name: $n conformers, clustered to $k_used"
    else
        echo "$name clustering_failed" >> cluster_failed.txt
    fi
done

echo "Summary: conformer_summary.csv | Kelley plots: kelley/ | problems: $(wc -l < cluster_failed.txt) (cluster_failed.txt)"
echo "Next: python $PIPE_DIR/03_make_gaussian_inputs.py"
