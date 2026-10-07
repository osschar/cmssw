#!/bin/bash
# Kalman hit acceptance scan: Policy::hit_chi2_cut x Policy::chi2_trk_fac.
#
# The acceptance chi2 is the residual against chi2_trk_fac^2 * C_trk + C_hit in
# the module's local frame; the hit covariance is taken as it is. The reference
# is the plain chi2 < 30 (trk_fac 1). Production defaults otherwise.
# Paired: same events, same seeds, one val_hit_chi2() call apart.
#
#   usage:  v2p2-eff-chi2.sh [n_events] [sample] [out_prefix]
#   CFGS="cut:fac ..." overrides the grid; the first entry is the reference.
#   ASSOC=mtv associates by MTV's rule (> 75 % of all hits true) instead of
#   quality-val's (>= 50 % of the non-seed hits).
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-30}
SAMPLE=${2:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev-rt.bin}
OUT=${3:-eff-chi2}
T=../RecoTracker/MkFitCore/standalone/test
CFGS=${CFGS:-"30:1 20:1 13.8:1 9.21:1 30:1.5 20:1.5 13.8:1.5 9.21:1.5 20:2 13.8:2 9.21:2 13.8:3 9.21:3"}
REF=${CFGS%% *}

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --read-cmssw-tracks --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command "val_assoc_mtv($([ "${ASSOC:-}" = mtv ] && echo true || echo false))"
     --shell-command 'val_eff_reset()'
     --shell-command "val_eff_ref(\"c$REF\")")

for c in $CFGS; do
  cut=${c%%:*}; fac=${c##*:}
  CMD+=(--shell-command "val_hit_chi2($cut, $fac)")
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)"
          --shell-command 's.ProcessEventStd()'
          --shell-command "val_eff_ev(s.event(), \"c$c\")")
  done
done
CMD+=(--shell-command "val_eff_report(\"$OUT\")")
echo .q | "${CMD[@]}"
