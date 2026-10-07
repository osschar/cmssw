#!/bin/bash
# Per-hit surface reference (module normal) vs the MkBins cylinder scaffold.
#
# The per-hit version uses the HIT'S OWN module normal, so tilted and flat layers
# need no branching and TBPS is not over-widened. Expect it to be at least as
# good as the cylinder version, and better in the OT barrel.
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-100}
SAMPLE=/foo/matevz/mic-dev/trackingNtuple_HLT_2026_March.bin
T=../RecoTracker/MkFitCore/standalone/test

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command 'val_seeds(0, true)')

# mkbins_cyl  perhit  extra_dq  label
CFG=("0 0 3.0 BASE_____none______dq3"
     "1 0 3.0 CYL______cylinder__dq3"
     "0 1 3.0 HIT______modnorm___dq3"
     "0 1 2.0 HIT______modnorm___dq2"
     "0 1 1.5 HIT______modnorm___dq1.5"
     "0 1 1.0 HIT______modnorm___dq1"
     "1 1 1.5 BOTH_____cyl+mod___dq1.5")

for c in "${CFG[@]}"; do
  set -- $c
  CMD+=(--shell-command "val_surf_q($1)" --shell-command "val_surf_q_hit($2)"
        --shell-command "val_extra_dq($3)"
        --shell-command "printf(\"@@CFG@@ %s\\n\", \"$4\")")
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)" --shell-command 's.ProcessEventHlt()')
  done
done
echo .q | "${CMD[@]}"
