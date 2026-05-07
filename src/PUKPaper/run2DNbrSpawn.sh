#!/bin/bash

# Script to run 2D Neighbour Spawn simulations in parallel for multiple densities
# L=128, steps=1000000, densities: 0.2, 0.4, 0.6, 0.8

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXECUTABLE="$SCRIPT_DIR/2DNeighbourSpawn"
L=64
STEPS=1000000
DENSITIES=(0.2 0.4 0.6 0.8)

# Compile if executable doesn't exist
if [ ! -f "$EXECUTABLE" ]; then
    echo "Compiling 2DNeighbourSpawn.cpp..."
    g++ -O3 -std=c++17 -o "$EXECUTABLE" "$SCRIPT_DIR/2DNeighbourSpawn.cpp"
    if [ $? -ne 0 ]; then
        echo "Compilation failed!"
        exit 1
    fi
    echo "Compilation successful!"
fi

echo "Starting parallel simulations..."
echo "L=$L, Steps=$STEPS, Densities: ${DENSITIES[@]}"
echo ""

# Run simulations in parallel, one per core
for density in "${DENSITIES[@]}"; do
    echo "Starting simulation with density=$density..."
    "$EXECUTABLE" "$L" "$density" "$STEPS" &
done

# Wait for all background processes to complete
wait
echo ""
echo "All simulations completed!"
