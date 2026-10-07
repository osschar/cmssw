#!/bin/bash
# Absolute build time of the pre-selection WINDOW configurations, production
# setup (in-layer search on, maxCandsPerSeed 3). Companion to the paired physics
# drivers v2p2-eff-dq.sh / v2p2-eff-phitrk.sh / v2p2-chop-phitrk.sh.
#
# The dq and phi timings recorded so far are PAIRED RELATIVE numbers taken in
# the trace build, valid as ratios only. This gives seconds comparable with
# v2p2-timing.sh's table: same sample, same inner/total definitions, minimum over
# repetitions, and the binary must be built with EVENT_RDF_TRACE OFF.
#
# Configurations are "label|mkFit options"; every option is a command-line one,
# so no shell and no ROOT is needed. Repetitions are INTERLEAVED (rep 1 of every
# config, then rep 2, ...) so slow drift in machine load lands on all of them.
#
#   usage:  v2p2-timing-win.sh [n_events] [n_reps] [sample]
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-100}
REPS=${2:-3}
SAMPLE=${3:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev.bin}
OUT=${TMPDIR:-/tmp}/mkfit-timing-win-$$

if ldd ./mkFit | grep -qi libCore; then
  echo "REFUSING: mkFit is linked against ROOT, so this is very likely a trace build." >&2
  exit 1
fi

BASE=(--geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
      --num-events "$N" --num-thr 1 --num-thr-ev 1 --num-iters-cmssw 1 --silent
      --build-mimi --build-mimi-v2p2 --v2p2-in-layer-comb 1 --max-cands-per-seed 3)

if [ -z "$CFGS" ]; then
  # The pre-2026-09-24 window, every knob pinned, against the current defaults.
  CFGS=$'flat_dq1.8|--v2p2-phi-per-hit 0 --v2p2-dphi-trk 1.0 --v2p2-dq-hit 1.8\ndef|'
fi

mkdir -p "$OUT"
echo "# window timing, $N events, $REPS reps interleaved, $(date -Is)"
echo "# production configuration: in-layer search on, maxCandsPerSeed 3"
echo
for ((r=1;r<=REPS;r++)); do
  while IFS='|' read -r lab opts; do
    [ -z "$lab" ] && continue
    ./mkFit "${BASE[@]}" $opts > "$OUT/$lab.$r.log" 2>&1
    inner=$(grep -m1 "Iteration 0 build time =" "$OUT/$lab.$r.log" | awk '{print $6}')
    total=$(grep -m1 "Total event loop time"    "$OUT/$lab.$r.log" | awk '{print $5}')
    echo "$lab rep$r inner=$inner total=$total"
  done <<< "$CFGS"
done

echo
printf "%-22s %10s %10s %10s %9s %s\n" config inner_s ms_per_ev total_s "vs first" options
REF=""
while IFS='|' read -r lab opts; do
  [ -z "$lab" ] && continue
  i=$(grep -h "Iteration 0 build time =" "$OUT/$lab".*.log | awk '{print $6}' | sort -g | head -1)
  t=$(grep -h "Total event loop time"    "$OUT/$lab".*.log | awk '{print $5}' | sort -g | head -1)
  [ -z "$REF" ] && REF=$i
  printf "%-22s %10.3f %10.1f %10.3f %8.1f%% %s\n" "$lab" "$i" \
         "$(echo "$i*1000/$N" | bc -l)" "$t" "$(echo "100*($i/$REF-1)" | bc -l)" "$opts"
done <<< "$CFGS"
echo
echo "# logs in $OUT"
