#!/bin/bash
# A/B for the surface-referenced q window (MkBins::surface_reference_dq).
#
# Hypothesis under test: the pre-selection q window is built from a covariance
# transported to a fixed PATH LENGTH and never referenced to the layer surface,
# which makes sigma_q too small by ~1/sin^2(theta). If so, turning the correction
# on should FLATTEN the eta dependence of the track term.
#
# Scope: this touches the WINDOW AND THE dq CUT ONLY. The Kalman update already
# carries the term (jacCurv2Loc's cosz), so chi2 is expected NOT to move.
#
# Chopped pT5 seeds, i.e. the control where the pixel hits are known findable.
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
     --shell-command 'val_seeds(0, true)')      # pT5, chop the pixel hits

for MODE in 0 1; do
  CMD+=(--shell-command "val_surf_q($MODE)" --shell-command 'val_sr_reset()')
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)"
          --shell-command 's.ProcessEventHlt()'
          --shell-command "val_sr_ev(s.event(), $i)")
  done
  CMD+=(--shell-command "val_sr_write(\"$T/val-surfq-$MODE.root\")")
done

echo .q | "${CMD[@]}"
