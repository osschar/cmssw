#!/bin/bash
# What EXTRA_DQ costs, at the production configuration (in-layer search, cap 3).
#
# The physics side is measured by v2p2-eff-dq.sh and is FLAT-to-slightly-better
# from 3.0 down to 1.5, so the whole case for narrowing rests on cost. This is
# the cost half.
#
# Same conventions as v2p2-timing.sh: inner is the forward search alone, total
# adds read + hit load + export, minimum over repetitions, and the binary must be
# built with EVENT_RDF_TRACE OFF.
#
#   usage:  v2p2-timing-dq.sh [n_events] [n_reps] [sample]
set -e
cd /foo/matevz/mic-dev/current/src/standalone
unset DISPLAY
export LD_LIBRARY_PATH=.
N=${1:-100}
REPS=${2:-3}
SAMPLE=${3:-/foo/matevz/mic-dev/ttbar-PU200-D121-C22-100ev.bin}
OUT=${TMPDIR:-/tmp}/mkfit-timing-dq-$$

if ldd ./mkFit | grep -qi libCore; then
  echo "REFUSING: mkFit is linked against ROOT, so this is very likely a trace build." >&2
  exit 1
fi

BASE=(--geom CMS-phase2 --seed-input cmssw --input-file "$SAMPLE"
      --num-events "$N" --num-thr 1 --num-thr-ev 1 --num-iters-cmssw 1 --silent
      --build-mimi --build-mimi-v2p2 --v2p2-in-layer-comb 1 --max-cands-per-seed 3
      --v2p2-phi-per-hit 0 --v2p2-dphi-trk 1.0)  # the phi window it was recorded with

mkdir -p "$OUT"
echo "# EXTRA_DQ timing, $N events, $REPS reps, $(date -Is)"
echo "# production configuration: in-layer search on, maxCandsPerSeed 3"
echo
for dq in ${DQLIST:-3.0 2.0 1.5 1.0}; do
  for ((r=1;r<=REPS;r++)); do
    ./mkFit "${BASE[@]}" --v2p2-extra-dq $dq > "$OUT/dq$dq.$r.log" 2>&1
    inner=$(grep -m1 "Iteration 0 build time =" "$OUT/dq$dq.$r.log" | awk '{print $6}')
    total=$(grep -m1 "Total event loop time"    "$OUT/dq$dq.$r.log" | awk '{print $5}')
    echo "dq=$dq rep$r inner=$inner total=$total"
  done
done

echo
printf "%-8s %10s %10s %10s %10s\n" EXTRA_DQ inner_s ms_per_ev total_s "vs dq3.0"
REF=""
for dq in ${DQLIST:-3.0 2.0 1.5 1.0}; do
  i=$(grep -h "Iteration 0 build time =" "$OUT/dq$dq".*.log | awk '{print $6}' | sort -g | head -1)
  t=$(grep -h "Total event loop time"    "$OUT/dq$dq".*.log | awk '{print $5}' | sort -g | head -1)
  [ -z "$REF" ] && REF=$i
  printf "%-8s %10.3f %10.1f %10.3f %9.1f%%\n" "$dq" "$i" \
         "$(echo "$i*1000/$N" | bc -l)" "$t" "$(echo "100*($i/$REF-1)" | bc -l)"
done
echo
echo "# logs in $OUT"
