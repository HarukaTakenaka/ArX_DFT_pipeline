# Settings shared by the pipeline scripts. Edit here, not in the scripts.

# --- Schrödinger -------------------------------------------------------------
# 2024-4: 2021-4 (the path in the original prep_MM_npsh.bash) is still installed on CHPC but can no
# longer check out a licence ("FATAL -15: Could not check out a license for mmlibs version 5.6", Nov 2026).
# Pin one exact version so every ArCl/ArI search uses the same OPLS parameters, even if CHPC
# changes its default module.
SCHRODINGER_MODULE="schrodinger/2024-4"
SCHRODINGER_PIN="/uufs/chpc.utah.edu/sys/installdir/schrodinger/2024-4"
export SCHRODINGER_LICENSE_RETRY=300m

# --- Conformer search and clustering (same as the ArBr library) ----------------
MAX_CONF=20          # cluster only if a molecule has MORE than this many conformers (SI: "greater than 20")
N_CLUSTERS=0         # 0 = number of clusters at the minimum Kelley penalty (SI). Any other n forces n clusters.

# --- Parallel MacroModel jobs -------------------------------------------------
# Defaults to the number of SLURM tasks. Each running job uses Schrödinger license tokens.
MAX_PARALLEL="${SLURM_NTASKS:-4}"

# --- Python (for reading the Kelley table and plotting in step 2, and steps 3 and 5) ----
# Step 3 needs RDKit; step 2 needs only the standard library (matplotlib optional, for plots).
PYTHON_MODULE="python/3.10"       # the old post-processing script used "ml python/3.10"; "" to skip
PYTHON="${PYTHON:-python3}"

load_python() {
    [[ -n "$PYTHON_MODULE" ]] && module load "$PYTHON_MODULE" 2>/dev/null || true
    command -v "$PYTHON" >/dev/null || { echo "ERROR: $PYTHON not found" >&2; exit 1; }
}

load_schrodinger() {
    module load "$SCHRODINGER_MODULE" 2>/dev/null || true
    if [[ -n "$SCHRODINGER_PIN" && -d "$SCHRODINGER_PIN" ]]; then
        export SCHRODINGER="$SCHRODINGER_PIN"
    elif [[ -n "$SCHRODINGER_PIN" ]]; then
        echo "WARNING: $SCHRODINGER_PIN not found; using the module default ($SCHRODINGER)." >&2
        echo "         Results may differ slightly from the ArBr library (different OPLS version)." >&2
    fi
    if [[ -z "${SCHRODINGER:-}" || ! -x "$SCHRODINGER/bmin" ]]; then
        echo "ERROR: Schrödinger not found (SCHRODINGER='${SCHRODINGER:-}')." >&2; exit 1
    fi
    echo "Using Schrödinger at $SCHRODINGER"
}
