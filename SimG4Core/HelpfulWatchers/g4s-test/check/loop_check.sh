#!/bin/bash
# loop_check.sh [LOOPS] (on uaf-4, after build.sh and check20.sh hgcal): re-simulates event 5 of check-hgcal-nosnitch,
# its engine state restored, LOOPS times (default 2) in one job with a looper (G4SNITCH_LOOPS, IterateNTimesLooper),
# step log on, in $G/loop<LOOPS>. Needs the patch to SimG4Core/Application that lets g4SimHits start a Geant4 run in
# each loop (without it the later loops simulate nothing: "Geometry not closed"). Checks:
# - each loop's hits, SimTracks and SimVertices (rerun.root, rerun001.root, ...) equal check-hgcal-rerun's (event 5
#   re-simulated once);
# - the snitch's step log of each loop equals the first loop's (the log cut at each "++ signal BeginOfRun");
# - the snitch's output of each loop (G4Snitch.root, then G4Snitch_run<k>.root) equals loop2-ttree/G4Snitch.root
#   (event 5 re-simulated once, without the patch), when that file is there.
G=/ceph/users/matevz/g4snitch-2026-11
C=$G/CMSSW_20_1_0_pre5/src/SimG4Core/HelpfulWatchers/g4s-test/check
N=${1:-2}
D=$G/loop$N
cd $G
rm -rf $D
G4SNITCH_INPUT=$G/check-hgcal-nosnitch/step1.root G4SNITCH_EVENTS=1:5 G4SNITCH_LOOPS=$N G4SNITCH_BFIELD=1 \
  timeout 1800 bash run20.sh $D rerunG4Snitch_D128.py > $D.log 2>&1
tail -1 $D/info.txt
echo "framework records: $(grep -c 'Begin processing' $D/e); snitch BeginOfRun: $(grep -c '++ signal BeginOfRun' $D/o)," \
     "BeginOfEvent: $(grep -c '++ signal BeginOfEvent' $D/o), EndOfRun: $(grep -c '++ signal EndOfRun' $D/o);" \
     "G4Exceptions: $(grep -c 'G4Exception-START' $D/e)"
grep 'Begin processing' $D/e | grep -o 'Event [0-9]*.*PST'
source env.sh > /dev/null 2>&1
for f in $(cd $D && ls rerun*.root); do
  echo "== $f against check-hgcal-rerun/rerun.root"
  python3 $C/g4s_hits_check.py $G/check-hgcal-rerun/rerun.root $D/$f 2>&1 | grep -v -i warning | tail -3
done
# The step log of each loop, from its BeginOfRun to the next.
awk -v d=$D '/^\+\+ signal BeginOfRun/ {n++} n > 0 {print > (d "/o.loop" n)}' $D/o
for k in $(seq 2 $N); do
  cmp $D/o.loop1 $D/o.loop$k && echo "step log of loop $k identical to loop 1's ($(wc -l < $D/o.loop1) lines)"
done
ls $D/G4Snitch*.root
if [ -f loop2-ttree/G4Snitch.root ]; then
  for f in $(cd $D && ls G4Snitch*.root); do
    echo "== the snitch's $f against loop2-ttree/G4Snitch.root (one pass, without the patch)"
    python3 $C/g4s_mt_check.py loop2-ttree/G4Snitch.root -- $D/$f 2>&1 | grep -v -i warning | tail -2
  done
fi
