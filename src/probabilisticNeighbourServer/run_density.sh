#!/bin/bash
# Density sweep at L=512 on both sides of p_c, for fig4 of paperFigures.py. Runs locally.
#   SEED_LO=1 SEED_HI=8 ./run_density.sh     the first half (default)
#   SEED_LO=9 SEED_HI=16 ./run_density.sh    the second half, later
#   DRY_RUN=1 ./run_density.sh               print the plan only
# On the server, point OUTDIR at the main sweep so its rho=0.2 runs are skipped:
#   OUTDIR=/nbi/home/rpw391/cell-disk/solarFlares/outputs SEED_LO=9 SEED_HI=16 nohup ./run_density.sh &> density.out &
# Resumes like run_server.sh: a run counts as done once its spotSize_*.tsv.gz exists.

set -u

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${OUTDIR:-$DIR/outputs}"
SRC="$DIR/probabilisticNeighbourServer.cpp"
EXECUTABLE="$DIR/probabilisticNeighbourServer"
LOG="$DIR/run_density.log"

L=512
STEPS=100000
NSNAP=8
SEED_LO=${SEED_LO:-1}
SEED_HI=${SEED_HI:-8}
DRY_RUN=${DRY_RUN:-0}
MAX_JOBS=${MAX_JOBS:-$(nproc --ignore=2)}   # leave 2 cores free

RHOS="0.9 0.8 0.6 0.4 0.2 0.1 0.05"   # densest first, so no long job is left stranded at the end
PS="0.3 0.5 0.8 0.9"

if [ ! -x "$EXECUTABLE" ] || [ "$SRC" -nt "$EXECUTABLE" ]; then
    echo "building..."
    g++ -O3 -march=native -std=c++17 -o "$EXECUTABLE" "$SRC" || exit 1
fi

run_simulation() {
    local L=$1 rho=$2 steps=$3 p=$4 seed=$5 nsnap=$6
    local pf rf
    pf=$(printf "%g" "$p"); rf=$(printf "%g" "$rho")   # match C++ double formatting
    local stem="L_${L}_rho_${rf}_p_${pf}_seed_${seed}"
    if [ -s "$OUTDIR/spotSize_${stem}.tsv.gz" ]; then  # written last, so it means done
        return 0
    fi
    "$EXECUTABLE" "$L" "$rho" "$steps" "$p" "$seed" "$nsnap" "$OUTDIR"
    if [ $? -ne 0 ]; then
        echo "FAILED L=$L rho=$rho p=$p seed=$seed" >&2
        return 1
    fi
    gzip -6 -f "$OUTDIR/snapshots_${stem}.tsv" "$OUTDIR/spotSize_${stem}.tsv" \
              "$OUTDIR/emission_${stem}.tsv" "$OUTDIR/maxSpot_${stem}.tsv" 2>/dev/null
}
export -f run_simulation
export EXECUTABLE OUTDIR

combinations=()
for rho in $RHOS; do for p in $PS; do for seed in $(seq "$SEED_LO" "$SEED_HI"); do
    combinations+=("$L $rho $STEPS $p $seed $NSNAP")
done; done; done

echo "${#combinations[@]} runs, ${MAX_JOBS} parallel, seeds ${SEED_LO}-${SEED_HI}, output -> $OUTDIR"
if [ "$DRY_RUN" != "0" ]; then
    exit 0
fi

mkdir -p "$OUTDIR"
echo "=== started $(date) with ${#combinations[@]} runs ===" | tee -a "$LOG"
printf '%s\n' "${combinations[@]}" \
    | xargs -n 6 -P "$MAX_JOBS" bash -c 'run_simulation "$@"' _
echo "=== finished $(date) ===" | tee -a "$LOG"
