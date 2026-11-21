#!/bin/bash
# check20.sh [GUN] (on uaf-4, after build.sh): run the D128 test configuration with the gun GUN (default cms: into most
# of CMS; also vertex, or hgcal for the gun at the HGCal face) twice, with the DetIds on and off, field at the path
# points on, in $G/check-GUN and $G/check-GUN-0, and check that:
# - the verbose logs are identical (the snitch does not change the simulation),
# - every PCaloHit and PSimHit collection is identical, and each snitch DetId is among its event's hits, by detector,
# - and list the sensitive steps by detector and subdetector, with their |eta| range.
G=/ceph/users/matevz/g4snitch-2026-11
GUN=${1:-cms}; [ "$GUN" = hgcal ] && GUNV= || GUNV=$GUN
C=$G/CMSSW_20_1_0_pre5/src/SimG4Core/HelpfulWatchers/g4s-test/check
pgrep -u $USER -x cmsRun > /dev/null && { echo "cmsRun is running; not starting"; exit 1; }
STDBUF=1 G4SNITCH_GUN=$GUNV G4SNITCH_BFIELD=1 bash $G/run20.sh $G/check-$GUN > $G/check-$GUN.log 2>&1 &
STDBUF=1 G4SNITCH_GUN=$GUNV G4SNITCH_BFIELD=1 G4SNITCH_DETID=0 bash $G/run20.sh $G/check-$GUN-0 > $G/check-$GUN-0.log 2>&1 &
wait
tail -1 $G/check-$GUN/info.txt; tail -1 $G/check-$GUN-0/info.txt
cmp $G/check-$GUN/o $G/check-$GUN-0/o && echo "verbose logs identical, DetIds on and off"
source $G/env.sh > /dev/null 2>&1
python3 $C/g4s_hits_check.py $G/check-$GUN-0/step1.root $G/check-$GUN/step1.root $G/check-$GUN/G4Snitch.root 2>&1 | grep -v -i warning
python3 $C/g4s_subdets.py $G/check-$GUN/G4Snitch.root 2>&1 | grep -v -i warning
