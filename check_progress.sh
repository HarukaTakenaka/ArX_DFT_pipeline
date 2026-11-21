#!/bin/bash
# How far the Gaussian jobs in one or more gaussian/ folders have got. Read-only: moves nothing,
# so it is safe to run while jobs are still running (unlike 05_check_logs.py).
#
#   bash $PIPE_DIR/check_progress.sh                      # the current folder
#   bash $PIPE_DIR/check_progress.sh ~/arx_r2/*/gaussian  # several folders at once
#   bash $PIPE_DIR/check_progress.sh -v ArCl/gaussian     # also list failed jobs with their error link
#
# finished    = all EXPECTED "Normal termination" lines (opt, freq, NBO, NMR)
# failed      = "Error termination" in the log (fix with 06_restart_from_log.py after the folder is done)
# in progress = log exists, not finished, no error (running, or preempted and waiting to restart)
# not started = .com with no .log yet
# Logs already sorted by 05_check_logs.py (logs/, term_error/, imag/) are counted too.

EXPECTED=4
verbose=0
[[ "${1:-}" == "-v" ]] && { verbose=1; shift; }
(( $# )) || set -- .

printf "%-40s %6s %9s %7s %12s %12s\n" folder total finished failed "in progress" "not started"
for dir in "$@"; do
    [[ -d "$dir" ]] || { echo "$dir: not a folder"; continue; }
    shopt -s nullglob
    coms=( "$dir"/*.com "$dir"/term_error/*.com "$dir"/imag/*.com "$dir"/resubmit/*.com )
    declare -A seen=()
    for c in "${coms[@]}"; do seen[$(basename "${c%.com}")]=1; done
    for l in "$dir"/*.log "$dir"/logs/*.log "$dir"/term_error/*.log "$dir"/imag/*.log; do seen[$(basename "${l%.log}")]=1; done
    total=${#seen[@]}; done_=0; failed=0; running=0; failed_list=()
    for name in "${!seen[@]}"; do
        log=""
        for d in "$dir" "$dir/logs" "$dir/term_error" "$dir/imag"; do
            [[ -f "$d/$name.log" ]] && { log="$d/$name.log"; [[ "$d" == "$dir" ]] && break; }
        done
        [[ -z "$log" ]] && continue
        n=$(grep -c "Normal termination" "$log")
        if (( n >= EXPECTED )); then ((done_++))
        elif grep -q "Error termination" "$log"; then
            ((failed++))
            link=$(grep -o "Error termination via Lnk1e in [^ ]*/l[0-9]*" "$log" | tail -1 | grep -o "l[0-9]*$")
            failed_list+=( "$name ${link:-?}" )
        else ((running++)); fi
    done
    unset seen
    printf "%-40s %6d %9d %7d %12d %12d\n" "$dir" "$total" "$done_" "$failed" "$running" $((total - done_ - failed - running))
    if (( verbose && failed )); then printf "    failed: %s\n" "${failed_list[@]}"; fi
done

if command -v squeue >/dev/null; then
    echo
    echo "SLURM: $(squeue -u "$USER" -h -t RUNNING -r | wc -l) running, $(squeue -u "$USER" -h -t PENDING -r | wc -l) waiting (all your jobs)"
fi
