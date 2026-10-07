#!/bin/bash
# A/B of the MkFinderV2p2 per-layer candidate policy (SESSIONS.md S4).
#
# Three switches, all default ON, each recovering the pre-2026-09-22 behaviour
# when set to 0:
#   --v2p2-wsr          set and act on the within-sensitive-region verdict
#   --v2p2-hole-limits  apply maxHolesPerCand / maxConsecHoles
#   --v2p2-stop-cuts    apply minPtCut and the looper stop at pull-in
#
# Reports the standard quality-val summary: found / pT10 / pT20 / no_mc_assoc
# and nH >= 80%. Same events and the same seeds in every configuration.
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-20}
SAMPLE=${2:-/foo/matevz/mic-dev/trackingNtuple_HLT_2026_March.bin}

run() {  # wsr holes stops label
  printf '@@CFG@@ %-34s ' "$4"
  echo .q | ./mkFit --geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE" \
      --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 \
      --num-iters-cmssw 1 --remove-dup --quality-val \
      --v2p2-wsr "$1" --v2p2-hole-limits "$2" --v2p2-stop-cuts "$3" 2>&1 \
    | sed -n '/Sum up of quality-val/,+2p' | tail -2 | tr '\n' ' '
  echo
}

run 0 0 0 "BASE__none"
run 1 0 0 "WSR__wsr"
run 0 1 0 "HOL__holes"
run 0 0 1 "STP__stops"
run 1 1 0 "WH___wsr+holes"
run 1 1 1 "ALL__wsr+holes+stops"
