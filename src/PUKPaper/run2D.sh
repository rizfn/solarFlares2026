#!/bin/bash

# Script to run 2D Neighbour Spawn and Random Spawn simulations in parallel for multiple densities
# L=128, steps=1000000, densities: 0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
L=128
STEPS=1000000
DENSITIES=(0.2 0.4 0.6 0.7 0.8 0.85 0.9)
MODELS=(2DNeighbourSpawn 2DRandomSpawn)

# Compile if executable doesn't exist
for model in "${MODELS[@]}"; do
    EXECUTABLE="$SCRIPT_DIR/$model"
    if [ ! -f "$EXECUTABLE" ]; then
        echo "Compiling $model.cpp..."
        g++ -O3 -std=c++17 -o "$EXECUTABLE" "$SCRIPT_DIR/$model.cpp"
        if [ $? -ne 0 ]; then
            echo "Compilation failed!"
            exit 1
        fi
        echo "Compilation successful!"
    fi
done

echo "Starting parallel simulations..."
echo "L=$L, Steps=$STEPS, Densities: ${DENSITIES[@]}"
echo ""

# Run simulations in parallel, one per core
for model in "${MODELS[@]}"; do
    for density in "${DENSITIES[@]}"; do
        echo "Starting $model simulation with density=$density..."
        "$SCRIPT_DIR/$model" "$L" "$density" "$STEPS" &
    done
done

# Wait for all background processes to complete
wait
echo ""
echo "All simulations completed!"
