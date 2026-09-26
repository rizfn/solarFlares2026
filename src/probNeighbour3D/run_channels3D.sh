#!/bin/bash
# Emission-channel sweep in 3D, to test whether psi (the rising exponent of 1-q)
# depends on dimension. L=64 has the same site count as the 2D L=512 runs.
set -u
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${OUT:-$DIR/outputs}"
EXE="$DIR/emissionChannels3D"
[ -x "$EXE" ] || g++ -O3 -march=native -std=c++17 -o "$EXE" "$DIR/emissionChannels3D.cpp" || exit 1
mkdir -p "$OUT"
PS="0.5 0.6 0.7 0.75 0.8 0.85 0.9 0.95 1.0"
jobs=()
for p in $PS; do for s in 1 2 3; do jobs+=("64 0.2 30000 $p $s"); done; done
for p in $PS; do for s in 1 2 3 4; do jobs+=("48 0.2 30000 $p $s"); done; done
run() {
    local stem="L_$1_rho_$(printf %g $2)_p_$(printf %g $4)_seed_$5"
    [ -s "$OUT/chQratio_${stem}.tsv" ] && return 0
    "$EXE" "$1" "$2" "$3" "$4" "$5" "$OUT"
}
export -f run; export EXE OUT
echo "${#jobs[@]} runs"
printf '%s\n' "${jobs[@]}" | xargs -n 5 -P 28 bash -c 'run "$@"' _
echo "done $(date)"
