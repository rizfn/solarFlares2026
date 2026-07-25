#!/bin/bash

# Generate all data for the figures. Writes to outputs/ (regenerable, not tracked).
# Grid A: density sweep at L=128 for histograms, snapshots, correlation length.
# Grid B: finite-size scaling near the critical point at rho=0.2.

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$DIR/probabilisticNeighbour"
OUT="$DIR/outputs"
STEPS=100000

if [ ! -f "$BIN" ]; then
    g++ -O3 -march=native -std=c++17 -o "$BIN" "$DIR/probabilisticNeighbour.cpp"
fi

for seed in 1 2 3; do
    for rho in 0.2 0.4 0.6 0.8; do
        for p in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0; do
            "$BIN" 128 "$rho" "$STEPS" "$p" "$seed" "$OUT" &
        done
    done
    for L in 64 128 256; do
        for p in 0.20 0.25 0.30 0.33 0.36 0.40 0.45 0.50 0.55; do
            "$BIN" "$L" 0.2 "$STEPS" "$p" "$seed" "$OUT" &
        done
    done
    wait
done
echo "All simulations completed!"
