#!/bin/bash
# LST-into-pixels, CHOPPED pT5 seeds: does the in-layer combinatorial search put
# the chopped hits back, and does the inward hole penalty matter?
#
# Why this sample and this control. A pT5's pixel hits were FOUND by the upstream
# reconstruction, so chopping them off gives a denominator of hits known to be
# findable, and the comparison is an exact (layer, index) match -- NO truth
# matching, so it is immune both to the nH>=80% inflation (raw reco hits against
# sim layers, which the combinatorial games by construction) and to the mc_match
# similarity trap. It is also the only place the head/body asymmetry of the hole
# penalty can be tested at all: this is the INWARD search.
#
# Read the numbers as an UPPER REFERENCE, not as a control group -- these are
# tracks the upstream reconstruction already reconstructed. Plain T5 (kind 1) is
# the population the search exists to rescue; run it too and quote separately.
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-50}
SAMPLE=${2:-/foo/matevz/mic-dev/trackingNtuple_HLT_2026_March.bin}
T=../RecoTracker/MkFitCore/standalone/test

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command 'val_seeds(0, true)')

#     comb slot mode eps    label
CFG=("0    0    0    0.99   BASE__best_hit"
     "1    0    0    0.99   COMB__linear"
     "1    1    0    0.99   COMB__linear__hole_slot"
     "1    0    1    0.9999 COMB__loglh"
     "1    1    1    0.9999 COMB__loglh__hole_slot")

for c in "${CFG[@]}"; do
  set -- $c
  CMD+=(--shell-command "val_in_layer_comb($1)"
        --shell-command "val_reserve_hole_slot($2)"
        --shell-command "val_score_mode($3, $4)"
        --shell-command 'val_te_reset()')
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)"
          --shell-command 's.ProcessEventHlt()'
          --shell-command 'val_chop_ev(s.event())'
          --shell-command 'val_te_ev(s.event())')
  done
  CMD+=(--shell-command "val_chop_report(\"$5\")"
        --shell-command "val_te_report(\"$5\")")
done
echo .q | "${CMD[@]}"
