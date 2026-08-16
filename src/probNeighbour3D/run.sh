#!/bin/bash
# 3D parameter sweep. Writes to outputs/ (regenerable, not tracked).
# L=32 is the main grid; L=48 at a few p is the size check.
# Safe to re-run: completed outputs are skipped.

set -u
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$DIR/probNeighbour3D"
OUT="$DIR/outputs"
STEPS=50000
SEEDS="1 2 3"
JOBS=28

mkdir -p "$OUT"
if [ ! -x "$BIN" ] || [ "$DIR/probNeighbour3D.cpp" -nt "$BIN" ]; then
    g++ -O3 -march=native -std=c++17 -o "$BIN" "$DIR/probNeighbour3D.cpp" || exit 1
fi

run() {
    local f="$OUT/snapshots_L_${1}_rho_${2}_p_${4}_seed_${5}.tsv"
    [ -s "$f" ] && [ "$(stat -c%s "$f")" -gt 10000 ] && return
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
    "$BIN" "$@" "$OUT" >/dev/null 2>&1 &
}

for seed in $SEEDS; do
    for L in 48 32; do
        for rho in 0.2 0.6; do
            for p in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
                [ "$L" = 48 ] && case "$p" in 0.0|0.5|0.8|1.0) ;; *) continue ;; esac
                run "$L" "$rho" "$STEPS" "$p" "$seed"
            done
        done
    done
done
wait
echo "3D sweep complete"
