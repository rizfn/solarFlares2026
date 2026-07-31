#!/bin/bash

# Generate all data for the figures. Writes to outputs/ (regenerable, not tracked).
# Grid A: density sweep at L=128 for histograms, snapshots, measured q.
# Grid B: finite-size check at rho=0.2 for the two endpoints and the middle.

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$DIR/addedBipole"
OUT="$DIR/outputs"
STEPS=100000

if [ ! -f "$BIN" ]; then
    g++ -O3 -march=native -std=c++17 -o "$BIN" "$DIR/addedBipole.cpp"
fi

for seed in 1 2 3; do
    for rho in 0.2 0.4 0.6 0.8; do
        for r in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
            "$BIN" 128 "$rho" "$STEPS" "$r" "$seed" "$OUT" &
        done
    done
    wait
    for L in 64 256 512; do
        for r in 0.0 0.5 1.0; do
            "$BIN" "$L" 0.2 "$STEPS" "$r" "$seed" "$OUT" &
        done
    done
    wait
done
echo "All simulations completed!"
