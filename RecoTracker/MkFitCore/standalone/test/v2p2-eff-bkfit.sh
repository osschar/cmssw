#!/bin/bash
# Does the BACKWARD FIT close the resolution gap to production?
#
# The recorded comparison "production is 37 % narrower in d(pT)/pT" was not a
# result: production's tracks are fully fitted and selected, ours came straight
# out of one forward search with no backward fit at all (Shell::ProcessEventStd
# takes m_backward_fit from Config::backwardFit and no scan passed it). This
# sets that half up.
#
# BACKWARD FIT ONLY, NOT BACKWARD SEARCH. SetupBackwardSearch() picks the search
# up at TB2S OTLayer4 / TEDD2 / TFPX6, chosen for LST T5 seeds that start in the
# outer tracker. Iteration 0 here is seeded by pixel quadruplets, so a backward
# search from those layers re-searches ground the forward pass already covered,
# from a seed region the track did not come from. Measured on one event it costs
# 164 found tracks of 3446.
#
# Both passes run in ONE process over the same events and the same seeds, so the
# delta is the fit and nothing else. cmsswTracks_ ride along with the first pass.
#
#   usage:  v2p2-eff-bkfit.sh [n_events] [sample] [out_prefix]
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-30}
SAMPLE=${2:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev-rt.bin}
OUT=${3:-eff-bkfit}
T=../RecoTracker/MkFitCore/standalone/test

# No --backward-fit on the command line: both passes are driven from the shell
# setters instead, so the two configurations differ in exactly one call.
CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --read-cmssw-tracks --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command 'val_eff_reset()'
     --shell-command 'val_eff_ref("nofit")'
     --shell-command 's.SetBackwardSearch(false)')

#     bkfit  label
CFG=("false  nofit"
     "true   bkfit")

for c in "${CFG[@]}"; do
  set -- $c
  CMD+=(--shell-command "s.SetBackwardFit($1)")
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)"
          --shell-command 's.ProcessEventStd()'
          --shell-command "val_eff_ev(s.event(), \"$2\")")
    if [ "$2" = "nofit" ]; then
      CMD+=(--shell-command 'val_eff_cmssw_ev(s.event(), "cmssw_V1")')
    fi
  done
done
CMD+=(--shell-command "val_eff_report(\"$OUT\")")
echo .q | "${CMD[@]}"
