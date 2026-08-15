#!/bin/bash
# Is the ordering at p<1 a real transition, or a crossover with the only singular
# point at p=1?
#
# The test is size scaling in the window 0.55 < p < 1, which no existing data covers.
# If q_V(p) is L-independent throughout, the transition sits at p=1. If curves for
# different L separate or cross near some p, that locates a real critical point.
#
# Usage:   ./run_transition.sh            # phase 1: L = 64,128,256,512   (~30-60 min)
#          ./run_transition.sh --big      # also L = 1024                 (~2-3 h more)
#
# Safe to re-run: completed outputs are skipped, so an interrupted run resumes.

set -u
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$DIR/probabilisticNeighbour"
OUT="$DIR/outputs"
LOG="$DIR/transition_run.log"
STEPS=100000
SEEDS="1 2 3 4 5 6"
JOBS=28                     # concurrent simulations; lower this if the machine is busy

# p grid: the untested window, plus anchors at each end that already have data
PS="0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0"
LS="64 128 256 512"
[ "${1:-}" = "--big" ] && LS="$LS 1024"

mkdir -p "$OUT"
if [ ! -x "$BIN" ] || [ "$DIR/probabilisticNeighbour.cpp" -nt "$BIN" ]; then
    echo "building..."
    g++ -O3 -march=native -std=c++17 -o "$BIN" "$DIR/probabilisticNeighbour.cpp" || exit 1
fi

# a run counts as done if its snapshot file exists and is non-trivial
done_already() {
    local f="$OUT/snapshots_L_${1}_rho_0.2_p_${2}_seed_${3}.tsv"
    [ -s "$f" ] && [ "$(stat -c%s "$f")" -gt 10000 ]
}

started=0; skipped=0
echo "=== transition sweep started $(date) ===" | tee -a "$LOG"
# big L first so the long jobs are not left stranded at the end
for L in $(echo "$LS" | tr ' ' '\n' | sort -rn); do
    for p in $PS; do
        for seed in $SEEDS; do
            if done_already "$L" "$p" "$seed"; then
                skipped=$((skipped + 1)); continue
            fi
            while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
            echo "  L=$L p=$p seed=$seed" | tee -a "$LOG"
            "$BIN" "$L" 0.2 "$STEPS" "$p" "$seed" "$OUT" >/dev/null 2>&1 &
            started=$((started + 1))
        done
    done
done
wait
echo "=== done $(date): $started run, $skipped already present ===" | tee -a "$LOG"
echo
echo "Now make the figures:"
echo "    python $DIR/plots.py                 # regenerates everything"
echo "or just the two that answer the question:"
echo "    python -c 'import plots; plots.plot_qv_finite_size(); plots.plot_correlation_function()'"
