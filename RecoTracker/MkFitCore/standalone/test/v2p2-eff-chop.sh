#!/bin/bash
# RESOLVED chop-recovery for the in-layer combinatorial search: the INWARD
# search, chopped pT5 seeds, HLT March sample.
#
# Truth-free by construction -- an exact (layer, index) match against the pixel
# hits the chop removed, which the upstream reconstruction had already found.
# That is what makes this sample usable at all: it predates the split-cluster
# arbitration fix, so its rec->sim links are stale and anything truth-matched on
# it is biased. Nothing here reads a truth link.
#
#   usage:  v2p2-eff-chop.sh [n_events] [sample] [out_prefix]
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-20}
SAMPLE=${2:-/foo/matevz/mic-dev/trackingNtuple_HLT_2026_March.bin}
OUT=${3:-eff-chop}
T=../RecoTracker/MkFitCore/standalone/test

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command 'val_seeds(0, true)'
     --shell-command 'val_chopres_reset()'
     --shell-command 'val_chopres_ref("base_c6")')

#     comb cap mode eps     label
CFG=("0    6   0    0.99    base_c6"
     "1    6   0    0.99    comb_c6"
     "1    3   0    0.99    comb_c3"
     "1    3   1    0.99    LL99_c3"
     "1    3   1    0.9999  LL9999_c3")

for c in "${CFG[@]}"; do
  set -- $c
  CMD+=(--shell-command "val_in_layer_comb($1)"
        --shell-command "val_max_cands($2)"
        --shell-command "val_score_mode($3, $4)")
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)"
          --shell-command 's.ProcessEventHlt()'
          --shell-command "val_chopres_ev(s.event(), \"$5\")")
  done
done
CMD+=(--shell-command "val_chopres_report(\"$OUT\")")
echo .q | "${CMD[@]}"
