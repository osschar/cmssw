#!/bin/bash
# PHI across the whole detector -- localising component 2, the eta-rising factor
# common to q and phi.
#
# phi is the only coordinate measurable everywhere: cluster spans are 2.2-3.3 in
# the precise direction in EVERY region (charge-interpolated, so Gaussian-ish)
# and the hit term is small, whereas q has 100 % single-cell clusters in every
# strip region and is therefore unmeasurable outside the pixels.
#
# THREE configurations, deliberately:
#   out  : real outward search (ProcessEventStd) -- real cmssw seeds, no invented
#          covariance. Covers the outer tracker. Does NOT scan the pixel barrel,
#          because that is where its seeds already are.
#   in   : chopped-pT5 inward search -- covers the pixel barrel.
#   sim1 : sim-seeded outward from layer 0 -- the only thing that covers the WHOLE
#   sim2   detector in one pass, but its seed covariance is INVENTED. Run at two
#          different fake priors so the prior-sensitivity can be shown rather
#          than assumed; deep layers should not care.
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-20}
SAMPLE=/foo/matevz/mic-dev/trackingNtuple_HLT_2026_March.bin
T=../RecoTracker/MkFitCore/standalone/test

base=(./mkFit --geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
      --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
      --shell-command 'gROOT->SetBatch(kTRUE)'
      --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
      --shell-command 'val_search_lite(1)'          # whole detector, skip exact solve
      --shell-command 'val_surf_q(true)')

run () {   # $1 = label, $2 = per-event shell command
  local CMD=("${base[@]}" --shell-command 'val_sr_reset()')
  [ -n "${3:-}" ] && CMD+=(--shell-command "$3")
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)" --shell-command "$2"
          --shell-command "val_sr_ev(s.event(), $i)")
  done
  CMD+=(--shell-command "val_sr_write(\"$T/val-phi-$1.root\")")
  echo .q | "${CMD[@]}"
}

run out  's.ProcessEventStd()'
run in   's.ProcessEventHlt()'          'val_seeds(0, true)'
run sim1 's.ProcessEventSimSeeded()'
