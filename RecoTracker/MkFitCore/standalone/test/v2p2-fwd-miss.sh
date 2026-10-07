#!/bin/bash
# Where the FORWARD search loses the sim track's own hit -- the taxonomy that
# exists for the inward search and did not for this one.
#
# The instrument is val_search_event(); only its truth join needed fixing. It
# keyed on TrCandMeta::global_seed, which is a seed's LABEL, and a label equals
# its index only after relabelSeedTracksSequentially() -- which ProcessEventHlt()
# calls and ProcessEventStd deliberately does not. It keys on TrCandMeta::seed
# into currentSeedTracks() now, which is right on both paths.
#
# val_search_lite(1) skips the exact double-precision helix-plane solve, which is
# what makes the per-hit analysis ~45 s/event. All layers are kept (bit 2 would
# restrict to the pixel barrel).
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-20}
SAMPLE=${2:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev-rt.bin}
OUT=${3:-val-fwd-miss.root}
T=../RecoTracker/MkFitCore/standalone/test

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command 'val_in_layer_comb(1)'
     --shell-command 'val_max_cands(6)'
     --shell-command 'val_search_lite(1)'
     --shell-command 'val_sr_reset()')
for ((i=1;i<=N;i++)); do
  CMD+=(--shell-command "s.GoToEvent($i)"
        --shell-command 's.ProcessEventStd()'
        --shell-command "val_sr_ev(s.event(), $i)")
done
CMD+=(--shell-command "val_sr_write(\"$OUT\")")
echo .q | "${CMD[@]}"
