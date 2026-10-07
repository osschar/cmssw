#!/bin/bash
# Does the surface-referenced window let EXTRA_DQ come down?
#
# EXTRA_DQ = 3 was chosen when the dq window was missing its surface reference
# and was therefore up to 9x too small at |eta| > 2 -- i.e. the factor was
# compensating for that, not for genuine tails. The recorded scan showing
# "3 -> 1 costs 12 % of good tracks" was taken with the broken window, so it has
# to be redone with the reference on.
#
# Measures the PHYSICS, not the window: per-event found / nH>=80% / no_mc_assoc
# from the standard quality output, plus the chopped-pT5 recovery efficiency,
# which is truth-free (exact (layer,index) match).
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-10}
SAMPLE=/foo/matevz/mic-dev/trackingNtuple_HLT_2026_March.bin
T=../RecoTracker/MkFitCore/standalone/test

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command 'val_seeds(0, true)')

# surf_q, extra_dq, label
CFG=("0 3.0 BASE__surfOFF_dq3"
     "1 3.0 SURF__surfON__dq3"
     "1 2.0 SURF__surfON__dq2"
     "1 1.5 SURF__surfON__dq1.5"
     "1 1.0 SURF__surfON__dq1")

for c in "${CFG[@]}"; do
  set -- $c
  CMD+=(--shell-command "val_surf_q($1)" --shell-command "val_extra_dq($2)"
        --shell-command "printf(\"@@CFG@@ %s\\n\", \"$3\")")
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)" --shell-command 's.ProcessEventHlt()')
  done
done
echo .q | "${CMD[@]}"
