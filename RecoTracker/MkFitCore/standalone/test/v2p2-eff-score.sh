#!/bin/bash
# What each term of the log-likelihood score actually does, per region.
#
# WHY AN ABLATION AND NOT A COMPARISON AGAINST THE LINEAR SCORE. The two differ
# in eps, in the chi2 weight (0.5 against 1.0) and in ln(det V) all at once, so a
# regional difference between them names no single term. Worse, rho and eps enter
# the likelihood IDENTICALLY -- both multiply n_hits -- so replacing log_rho by a
# constant is exactly a global shift of eps. The only thing rho can do that eps
# cannot is VARY, across layers and eta.
#
# Hence: hold eps fixed and replace each term by ITS OWN MEASURED MEAN, which
# removes that term's variation and leaves the average hit-versus-hole balance
# untouched. The means come from val_score_term_stats(), measured on this sample:
#   mean ln(rho)   =  3.8482   (807190 hits, 3 events)
#   mean ln(det V) = -13.5346
# Re-measure them if the sample or the windows change.
#
#   usage:  v2p2-eff-score.sh [n_events] [sample] [out_prefix]
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-30}
SAMPLE=${2:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev.bin}
OUT=${3:-eff-score}
T=../RecoTracker/MkFitCore/standalone/test
RHO=3.8482
DETV=-13.5346

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command 'val_in_layer_comb(1)'
     --shell-command 'val_max_cands(3)'
     --shell-command 'val_eff_reset()'
     --shell-command 'val_eff_ref("lin")')

#     mode eps     rho   detv  label
CFG=("0    0.99    1     1     lin"
     "1    0.9999  1     1     LL_e9999"
     "1    0.9999  0     1     e9999_flatrho"
     "1    0.99    1     1     LL_e99"
     "1    0.99    0     1     e99_flatrho"
     "1    0.99    1     0     e99_flatdetv"
     "1    0.99    0     0     e99_flatboth")

for c in "${CFG[@]}"; do
  set -- $c
  CMD+=(--shell-command "val_score_mode($1, $2)"
        --shell-command "val_score_terms($3, $RHO, $4, $DETV)")
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)"
          --shell-command 's.ProcessEventStd()'
          --shell-command "val_eff_ev(s.event(), \"$5\")")
  done
done
CMD+=(--shell-command "val_eff_report(\"$OUT\")")
echo .q | "${CMD[@]}"
