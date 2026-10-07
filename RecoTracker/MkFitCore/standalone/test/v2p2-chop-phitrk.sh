#!/bin/bash
# Per-hit phi extent on the INWARD search: chopped pT5 into the pixels.
#
# Companion to v2p2-eff-phitrk.sh, which measures the forward search only. The
# inward search is where windows have historically bound -- EXTRA_DQ = 1.0 cost
# it -92 tracks where the forward search shrugged -- because it runs into the
# pixels, which are covariance-dominated rather than containment-dominated. So a
# phi window carried by dphi_track alone is tested hardest here.
#
# Metrics are the two that cannot be gamed by taking more hits: val_chop_report
# is an exact (layer, index) match against the hits the chop removed, no truth;
# val_te_report counts only candidate hits whose mcTrackID is the sim label.
# Production configuration: in-layer search on, linear score, no hole slot.
#
#   usage:  v2p2-chop-phitrk.sh [n_events] [sample]
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-50}
SAMPLE=${2:-/foo/matevz/mic-dev/trackingNtuple_HLT_2026_March.bin}
T=../RecoTracker/MkFitCore/standalone/test

# label:per_hit_fac:trk_fac[:dq_hit_fac[:q_extra_bins[:phi_extra_bins[:precut]]]]
# -- per_hit_fac 0 means the flat 0.0246 rad constant; dq_hit_fac defaults to
# 1.8, with dq_trk_fac fixed at 1.5; both fetch margins default to one whole
# bin; precut is qphi (default, as in production), q, phi or 0 for the line pre-cuts
CFGS=${CFGS:-"flat:0:1.0 ph3_trk1.0:3:1.0 ph3_trk1.5:3:1.5 ph3_trk2.0:3:2.0 ph3_trk3.0:3:3.0"}

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command 'val_seeds(0, true)'
     --shell-command 'val_in_layer_comb(1)'
     --shell-command 'val_reserve_hole_slot(0)'
     --shell-command 'val_score_mode(0, 0.99)')

for cfg in $CFGS; do
  IFS=: read -r lab ph trk dqh qeb peb pc <<< "$cfg"
  case "${pc:-qphi}" in q) pcq=1 pcp=0 ;; phi) pcq=0 pcp=1 ;; qphi) pcq=1 pcp=1 ;; *) pcq=0 pcp=0 ;; esac
  CMD+=(--shell-command "val_precut($pcq, $pcp)")
  CMD+=(--shell-command "val_dq(1.5, ${dqh:-1.8}, ${qeb:-1})")
  if [ "$ph" = "0" ]; then
    CMD+=(--shell-command 'val_phi_per_hit(false, 1.0)')
  else
    CMD+=(--shell-command "val_phi_per_hit(true, $ph)")
  fi
  CMD+=(--shell-command "val_dphi($trk, 0.0246, ${peb:-1})"
        --shell-command 'val_te_reset()')
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)"
          --shell-command 's.ProcessEventHlt()'
          --shell-command 'val_chop_ev(s.event())'
          --shell-command 'val_te_ev(s.event())')
  done
  CMD+=(--shell-command "val_chop_report(\"$lab\")"
        --shell-command "val_te_report(\"$lab\")")
done
echo .q | "${CMD[@]}"
