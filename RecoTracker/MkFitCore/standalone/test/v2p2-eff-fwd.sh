#!/bin/bash
# RESOLVED efficiency for the in-layer combinatorial search: where does the
# +4.5 % actually land?  Outward search, D121 PU200, offline CMSSW seeds.
#
# Everything measured for this work so far is a scalar summed over a run. This
# runs every configuration over the SAME events with the SAME seeds in ONE
# process and reports the per-sim-track efficiency binned in |eta|, pT and the
# sim track's own hits-per-layer, with the paired per-event error bar.
#
# The metric is quality-val's association rule (2*mccount >= nCandHits over the
# non-seed hits), turned round so the denominator is SIM tracks a seed points
# at.  nH >= 80 % is deliberately not used anywhere: it tests raw reco HITS
# against sim LAYERS and so rewards exactly what the in-layer search changes.
#
#   usage:  v2p2-eff-fwd.sh [n_events] [sample] [out_prefix]
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-30}
SAMPLE=${2:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev.bin}
OUT=${3:-eff-fwd}
T=../RecoTracker/MkFitCore/standalone/test

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command 'val_eff_reset()'
     --shell-command 'val_eff_ref("base_c6")')

#     comb cap mode eps     label
CFG=("0    6   0    0.99    base_c6"
     "0    3   0    0.99    base_c3"
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
          --shell-command 's.ProcessEventStd()'
          --shell-command "val_eff_ev(s.event(), \"$5\")")
  done
done
CMD+=(--shell-command "val_eff_report(\"$OUT\")")
echo .q | "${CMD[@]}"
