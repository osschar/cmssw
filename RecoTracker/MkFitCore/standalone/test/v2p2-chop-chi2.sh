#!/bin/bash
# Kalman hit acceptance on the INWARD search: chopped pT5 into the pixels.
# Companion to v2p2-eff-chi2.sh (forward). Production defaults otherwise; only
# Policy::hit_chi2_cut and Policy::chi2_trk_fac move.
# Metrics: val_chop_report, exact (layer, index) match against the chopped hits,
# no truth; val_te_report, truth-matched candidate hits only.
#
#   usage:  v2p2-chop-chi2.sh [n_events] [sample]
#   CFGS="cut:fac ..." overrides the list.
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-50}
SAMPLE=${2:-/foo/matevz/mic-dev/trackingNtuple_HLT_2026_March.bin}
T=../RecoTracker/MkFitCore/standalone/test
CFGS=${CFGS:-"30:1 13.8:2"}

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command 'val_seeds(0, true)')

for c in $CFGS; do
  cut=${c%%:*}; fac=${c##*:}
  CMD+=(--shell-command "val_hit_chi2($cut, $fac)")
  CMD+=(--shell-command 'val_te_reset()')
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)"
          --shell-command 's.ProcessEventHlt()'
          --shell-command 'val_chop_ev(s.event())'
          --shell-command 'val_te_ev(s.event())')
  done
  CMD+=(--shell-command "val_chop_report(\"c$c\")"
        --shell-command "val_te_report(\"c$c\")")
done
echo .q | "${CMD[@]}"
