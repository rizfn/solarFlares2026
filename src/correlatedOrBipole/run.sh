#!/bin/bash
# Sweep matching probabilisticNeighbour/run_channels.sh, so the two models can be
# compared p by p at the same L, rho and sweep count. The question is whether the
# U in 1 - q(s) survives when the uncorrelated channel injects a bipole instead of
# two independent monomers.
set -u
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${OUT:-$DIR/outputs}"
EXE="$DIR/correlatedOrBipole"
[ -x "$EXE" ] || g++ -O3 -march=native -std=c++17 -o "$EXE" "$DIR/correlatedOrBipole.cpp" || exit 1
mkdir -p "$OUT"

PS="0.0 0.3 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0"
# cheapest L first: the 256 block finishes in minutes and already answers the
# question, while the 1024 block only sharpens it
jobs=()
for p in $PS; do for s in 1 2 3 4 5 6; do jobs+=("256 0.2 30000 $p $s"); done; done
for p in $PS; do for s in 1 2 3 4; do jobs+=("512 0.2 30000 $p $s"); done; done
for p in $PS; do for s in 1 2 3; do jobs+=("1024 0.2 20000 $p $s"); done; done

run() {
    local stem="L_$1_rho_$(printf %g $2)_p_$(printf %g $4)_seed_$5"
    [ -s "$OUT/cbQratio_${stem}.tsv" ] && return 0   # already done
    "$EXE" "$1" "$2" "$3" "$4" "$5" "$OUT" 2>/dev/null
}
export -f run; export EXE OUT
echo "${#jobs[@]} runs"
printf '%s\n' "${jobs[@]}" | xargs -n 5 -P 30 bash -c 'run "$@"' _
echo "done $(date)"
