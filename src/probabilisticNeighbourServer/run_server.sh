#!/bin/bash
# Parameter sweep for the server. See readme.md for the grid, cost and disk estimates.
#   DRY_RUN=1 ./run_server.sh    print the plan only
#   SEED_MULT=2 ./run_server.sh  scale every seed count

set -u

OUTDIR="${OUTDIR:-/nbi/home/rpw391/cell-disk/solarFlares/outputs}"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$DIR/probabilisticNeighbourServer.cpp"
EXECUTABLE="$DIR/probabilisticNeighbourServer"
LOG="$DIR/run_server.log"

STEPS=100000
SEED_MULT=${SEED_MULT:-1}
DRY_RUN=${DRY_RUN:-0}
MAX_JOBS=$(nproc --ignore=2)   # leave 2 cores free

nsnap_for() { case "$1" in 128) echo 4;; 256) echo 6;; 512) echo 8;; *) echo 6;; esac; }
seeds() { seq 1 $(( $1 * SEED_MULT )); }

if [ "$DRY_RUN" = "0" ]; then
    mkdir -p "$OUTDIR" || { echo "cannot create $OUTDIR" >&2; exit 1; }
fi

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

P_FINE="0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0"
P_COARSE="0.0 0.2 0.4 0.5 0.6 0.7 0.8 0.85 0.9 0.95 1.0"
P_TRANS="0.8 0.9 0.95 1.0"

combinations=()   # largest L first, so no long job is left stranded at the end
if [ "${BIG:-0}" != "0" ]; then   # ~30 h per run, off by default
    for seed in $(seeds 3); do for p in $P_TRANS; do combinations+=("2048 0.2 $STEPS $p $seed $(nsnap_for 2048)"); done; done
fi
for seed in $(seeds 16); do for p in $P_COARSE; do combinations+=("1024 0.2 $STEPS $p $seed $(nsnap_for 1024)"); done; done
for seed in $(seeds 24); do for p in $P_FINE;   do combinations+=("512 0.2 $STEPS $p $seed $(nsnap_for 512)"); done; done
for seed in $(seeds 32); do for rho in 0.2 0.6; do for p in $P_FINE; do combinations+=("256 $rho $STEPS $p $seed $(nsnap_for 256)"); done; done; done
for seed in $(seeds 64); do for rho in 0.2 0.4 0.6 0.8; do for p in $P_FINE; do combinations+=("128 $rho $STEPS $p $seed $(nsnap_for 128)"); done; done; done

echo "${#combinations[@]} runs, ${MAX_JOBS} parallel, output -> $OUTDIR"
if [ "$DRY_RUN" != "0" ]; then
    printf '%s\n' "${combinations[@]}" | awk '{print $1}' | sort -n | uniq -c \
        | awk '{printf "  L=%-5s %5d runs\n", $2, $1}'
    exit 0
fi

echo "=== started $(date) with ${#combinations[@]} runs ===" | tee -a "$LOG"
printf '%s\n' "${combinations[@]}" \
    | xargs -n 6 -P "$MAX_JOBS" bash -c 'run_simulation "$@"' _
echo "=== finished $(date) ===" | tee -a "$LOG"
