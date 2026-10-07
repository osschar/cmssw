#!/bin/bash
# Hits per OT layer on a track, for v2p2's final tracks, the production tracks
# in the sample and the sim tracks v2p2 found. Production v2p2 configuration
# (in-layer search on, maxCandsPerSeed 3), then a second pass with the search
# off at the same cap. See an-ot-nhits.C.
#
#   usage:  v2p2-ot-nhits.sh [n_events] [sample] [out.root]
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-30}
SAMPLE=${2:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev-rt.bin}
OUT=${3:-ot-nhits.root}
T=../RecoTracker/MkFitCore/standalone/test
CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --read-cmssw-tracks --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command "gROOT->ProcessLine(\".L $T/an-ot-nhits.C\")"
     --shell-command 'val_in_layer_comb(1)'
     --shell-command 'val_max_cands(3)'
     --shell-command 'ot_nh_reset()')
for ((i=1;i<=N;i++)); do
  CMD+=(--shell-command "s.GoToEvent($i)" --shell-command 's.ProcessEventStd()' --shell-command 'ot_nh_ev(s.event(), 0)')
done
# second pass: one best hit per layer, same cap
CMD+=(--shell-command 'val_in_layer_comb(0)')
for ((i=1;i<=N;i++)); do
  CMD+=(--shell-command "s.GoToEvent($i)" --shell-command 's.ProcessEventStd()' --shell-command 'ot_nh_ev(s.event(), 1)')
done
CMD+=(--shell-command "ot_nh_write(\"$OUT\")")
echo .q | "${CMD[@]}"
