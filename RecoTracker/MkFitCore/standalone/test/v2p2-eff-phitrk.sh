#!/bin/bash
# Per-hit phi extent: how much TRACK window does it need? Forward search.
#
# With Config::V2p2::Window::phi_per_hit on, the phi cut is
#
#   ddphi < trk_fac * dphi_track + hit_fac * hit_phi_half_extent
#
# and the hit term is ~4e-5 rad in TB2S against the flat 0.0246 it replaces, so
# the cut is carried by the track term alone. Measured so far at hit_fac 3:
# trk 1 costs -41 found tracks (4.2 sigma), trk 2 recovers it (+1), trk 3 and 5
# add nothing. This fills the 1 -> 2 interval, where the timing is expected to
# prefer the lower value. The eta-binned phi pull (1.25 below |eta| 0.8, 2.2-3.4
# above) says the answer should differ by region -- read the region rows.
#
# Paired: same events, same seeds, one val_dphi() call apart. Every configuration
# sets all window knobs itself, so the driver does not depend on the defaults.
# "flat" is the flat 0.0246 rad constant with dq_hit 1.8, the window before the
# per-hit extent, and is the reference; ph3_trk1.0 and ph3_trk2.0 reproduce
# recorded points.
#
#   usage:  v2p2-eff-phitrk.sh [n_events] [sample] [out_prefix]
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-30}
SAMPLE=${2:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev-rt.bin}
OUT=${3:-eff-phitrk}
T=../RecoTracker/MkFitCore/standalone/test

# label:per_hit_fac:trk_fac[:dq_hit_fac[:q_extra_bins[:phi_extra_bins[:precut]]]]
# -- per_hit_fac 0 means the flat 0.0246 rad constant; dq_hit_fac defaults to
# 1.8, with dq_trk_fac fixed at 1.5; both fetch margins default to one whole
# bin; precut is qphi (default, as in production), q, phi or 0 for the line pre-cuts
CFGS=${CFGS:-"flat:0:1.0 ph3_trk1.0:3:1.0 ph3_trk1.25:3:1.25 ph3_trk1.5:3:1.5 ph3_trk1.75:3:1.75 ph3_trk2.0:3:2.0"}

CMD=(./mkFit --geom CMS-phase2 --seed-input cmssw --read-cmssw-tracks --input-file "$SAMPLE"
     --num-events "$N" --num-thr 1 --build-mimi --build-mimi-v2p2 --shell
     --shell-command 'gROOT->SetBatch(kTRUE)'
     --shell-command "gROOT->ProcessLine(\".L $T/val-prop.C\")"
     --shell-command 'val_eff_reset()'
     --shell-command 'val_eff_ref("flat")')

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
  CMD+=(--shell-command "val_dphi($trk, 0.0246, ${peb:-1})")
  for ((i=1;i<=N;i++)); do
    CMD+=(--shell-command "s.GoToEvent($i)"
          --shell-command 's.ProcessEventStd()'
          --shell-command "val_eff_ev(s.event(), \"$lab\")")
  done
done
CMD+=(--shell-command "val_eff_report(\"$OUT\")")
echo .q | "${CMD[@]}"
