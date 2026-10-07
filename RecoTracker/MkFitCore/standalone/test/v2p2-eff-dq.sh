#!/bin/bash
# EXTRA_DQ scan on the CURRENT forward configuration (in-layer search, cap 3).
#
# WHY THIS IS A BASELINE AND NOT AN ANSWER. g_v2p2_extra_dq multiplies BOTH terms
# of the pre-selection cut (MkFinderV2p2.cc:1329):
#
#   ddq < EXTRA_DQ * dq_trk + EXTRA_DQ * DDQ_PRESEL_FAC * hit_q_half_length
#
# The first is trust in the propagated covariance; the second is containment of a
# geometric extent that is known exactly. Scanning the single global moves them
# together, so it can only find the best compromise, not the right shape. What it
# IS good for: the curve any layer/step-dependent form has to beat, and a
# re-check of the recorded "2.0 is slightly better than 3.0" -- which was
# measured on the chopped-pT5 INWARD search in the best-hit era, i.e. a different
# direction, a different sample and a different algorithm.
#
# Paired: same events, same seeds, one val_extra_dq() call apart. The phi window
# is pinned to the flat constant it was recorded with, since the default moved.
#
#   usage:  v2p2-eff-dq.sh [n_events] [sample] [out_prefix]
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-30}
SAMPLE=${2:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev-rt.bin}
OUT=${3:-eff-dq}
T=../RecoTracker/MkFitCore/standalone/test

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --read-cmssw-tracks --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command 'val_phi_per_hit(false, 1.0)'
     --shell-command 'val_dphi(1.0, 0.0246, 1)'
     --shell-command 'val_eff_reset()'
     --shell-command "val_eff_ref(\"dq${DQREF:-3.0}\")")

for dq in ${DQLIST:-3.0 2.5 2.0 1.5}; do
  CMD+=(--shell-command "val_extra_dq($dq)")
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)"
          --shell-command 's.ProcessEventStd()'
          --shell-command "val_eff_ev(s.event(), \"dq$dq\")")
    if [ "$dq" = "3.0" ]; then
      CMD+=(--shell-command 'val_eff_cmssw_ev(s.event(), "cmssw_V1")')
    fi
  done
done
CMD+=(--shell-command "val_eff_report(\"$OUT\")")
echo .q | "${CMD[@]}"
