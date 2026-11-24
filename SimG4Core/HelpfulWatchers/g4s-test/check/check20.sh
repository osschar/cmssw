#!/bin/bash
# check20.sh [GUN] (on uaf-4, after build.sh): checks of the snitch with the D128 test configuration and the gun GUN
# (default cms: into most of CMS; also vertex, or hgcal for the gun at the HGCal face), field at the path points on.
# Runs, in $G: check-GUN (snitch, DetIds on), check-GUN-0 (DetIds off), check-GUN-nosnitch (no snitch), and
# check-GUN-rerun (event 5 of the no-snitch sample re-simulated with the snitch on, the engine states restored, and
# only one particle's subtree written: the tracked non-primary of event 5 with the largest subtree). Checks:
# - the snitch's verbose logs are identical with the DetIds on and off;
# - all hits, SimTracks and SimVertices are identical with the snitch on and without it (the snitch is passive);
# - the re-run's hits, SimTracks and SimVertices equal the original's, and its subtree equals the same subtree in
#   check-GUN's snitch file (a single event of a sample can be re-run with the snitch);
# - each snitch DetId is among its event's hits, by detector, and the sensitive steps by detector and subdetector.
G=/ceph/users/matevz/g4snitch-2026-11
GUN=${1:-cms}; [ "$GUN" = hgcal ] && GUNV= || GUNV=$GUN
C=$G/CMSSW_20_1_0_pre5/src/SimG4Core/HelpfulWatchers/g4s-test/check
# Only this area's jobs count (other sessions run cmsRun under the same account elsewhere on uaf-4).
for p in $(pgrep -u $USER -x cmsRun); do
  case $(readlink /proc/$p/cwd) in $G/*) echo "cmsRun $p is running in $G; not starting"; exit 1 ;; esac
done
STDBUF=1 G4SNITCH_GUN=$GUNV G4SNITCH_BFIELD=1 bash $G/run20.sh $G/check-$GUN > $G/check-$GUN.log 2>&1 &
STDBUF=1 G4SNITCH_GUN=$GUNV G4SNITCH_BFIELD=1 G4SNITCH_DETID=0 bash $G/run20.sh $G/check-$GUN-0 > $G/check-$GUN-0.log 2>&1 &
G4SNITCH_GUN=$GUNV G4SNITCH_OFF=1 bash $G/run20.sh $G/check-$GUN-nosnitch > $G/check-$GUN-nosnitch.log 2>&1 &
wait
for d in check-$GUN check-$GUN-0 check-$GUN-nosnitch; do tail -1 $G/$d/info.txt; done
cmp $G/check-$GUN/o $G/check-$GUN-0/o && echo "verbose logs identical, DetIds on and off"
source $G/env.sh > /dev/null 2>&1
echo "== the snitch is passive: without it (A) against with it (B)"
python3 $C/g4s_hits_check.py $G/check-$GUN-nosnitch/step1.root $G/check-$GUN/step1.root 2>&1 | grep -v -i warning
ID=$(python3 $C/g4s_subtree_check.py pick $G/check-$GUN/G4Snitch.root 4 2>/dev/null | tail -1)
echo "== re-run of event 5 of the no-snitch sample with the snitch on, writing the subtree of Geant4 id $ID"
G4SNITCH_INPUT=$G/check-$GUN-nosnitch/step1.root G4SNITCH_EVENTS=1:5 G4SNITCH_SELECT_ID=$ID G4SNITCH_BFIELD=1 \
  bash $G/run20.sh $G/check-$GUN-rerun rerunG4Snitch_D128.py > $G/check-$GUN-rerun.log 2>&1
tail -1 $G/check-$GUN-rerun/info.txt
python3 $C/g4s_hits_check.py $G/check-$GUN-nosnitch/step1.root $G/check-$GUN-rerun/rerun.root 2>&1 | grep -v -i warning
python3 $C/g4s_subtree_check.py $G/check-$GUN/G4Snitch.root 4 $ID $G/check-$GUN-rerun/G4Snitch.root 2>&1 | grep -v -i warning
echo "== DetIds and coverage"
python3 $C/g4s_hits_check.py $G/check-$GUN-0/step1.root $G/check-$GUN/step1.root $G/check-$GUN/G4Snitch.root 2>&1 | grep -v -i warning
python3 $C/g4s_subdets.py $G/check-$GUN/G4Snitch.root 2>&1 | grep -v -i warning
