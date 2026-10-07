#!/bin/bash
# maxCandsPerSeed scan for the in-layer combinatorial search, resolved.
#
# 6 is the phase-2 default (CMS-phase2.cc SetupIterationParams); 3 is the
# operating point the scalar measurement preferred on time alone. This asks
# whether any REGION prefers a wider beam, which a summed counter cannot.
#
# The base row is the search OFF at the old default, so every delta is against
# what production does today.
#
#   usage:  v2p2-eff-cap.sh [n_events] [sample] [out_prefix]
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-30}
SAMPLE=${2:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev-rt.bin}
OUT=${3:-eff-cap}
T=../RecoTracker/MkFitCore/standalone/test

# --read-cmssw-tracks gives the PRODUCTION reference: whatever tracking the job
# that wrote the ntuple ran, which for these samples is mkFit V1 with
# prop-to-plane and selectHitIndicesV2.
#
# IT REQUIRES THE '-rt' SAMPLE, hence the default above. An earlier comment here
# claimed the flag costs nothing when the section is absent; it does not --
# DataFile::openRead() prints "Reading of CmsswTracks requested but data not
# available on file." and calls exit(1) (Event.cc:1320). The plain
# ttbar-PU200-D121-C22-100ev.bin carries Seeds + BeamSpot + SimHitStates and no
# CmsswTracks, so it aborts before the first event.
CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --read-cmssw-tracks --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command 'val_eff_reset()'
     --shell-command 'val_eff_ref("base_c6")')

#     comb cap  label
CFG=("0    6    base_c6"
     "1    3    comb_c3"
     "1    4    comb_c4"
     "1    5    comb_c5"
     "1    6    comb_c6")

for c in "${CFG[@]}"; do
  set -- $c
  CMD+=(--shell-command "val_in_layer_comb($1)"
        --shell-command "val_max_cands($2)")
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)"
          --shell-command 's.ProcessEventStd()'
          --shell-command "val_eff_ev(s.event(), \"$3\")")
    # cmsswTracks_ do not depend on the configuration, so they ride along with
    # the first pass and cost no extra finding. They do need ProcessEventStd to
    # have run: the seed-hit exclusion reads Event::currentSeedTracks().
    if [ "$3" = "base_c6" ]; then
      CMD+=(--shell-command 'val_eff_cmssw_ev(s.event(), "cmssw_V1")')
    fi
  done
done
CMD+=(--shell-command "val_eff_report(\"$OUT\")")
echo .q | "${CMD[@]}"
