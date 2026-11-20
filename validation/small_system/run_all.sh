#!/bin/bash
# Runs 22 SEGC-WL and 88 1/t GC-WL runs of the small test system on 11 cores (about 30 min),
# then compares them. Needs GCWL.jl (package GCWL); set GCWL to its path (default: cloned next to this repository as ../GCWL.jl).
# Usage, from this folder:  bash run_all.sh
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"; SEGC="$(cd "$HERE/../.." && pwd)"; GC="$(cd "${GCWL:-$SEGC/../GCWL.jl}" && pwd)"
{ for i in $(seq 1 22); do echo "$SEGC run_segc.jl $((3000+i)) $HERE/runs/segc$i"; done
  for i in $(seq 1 88); do echo "$GC run_gcwl.jl $((4000+i)) $HERE/runs/gcwl$i"; done; } |
  xargs -P 11 -L 1 sh -c 'julia --project=$0 '"$HERE"'/$1 $2 $3 40 4.0 1.2 > /dev/null 2>&1'
for d in "$HERE"/runs/*; do find "$d" -type f ! -name lnQ.csv ! -name time.txt -delete; done
julia "$HERE/compare.jl" "$HERE/runs" | tee "$HERE/results.txt"
