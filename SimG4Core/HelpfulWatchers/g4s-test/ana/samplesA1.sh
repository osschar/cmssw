#!/bin/bash
# ceph-samplesA1.sh (copied to /ceph/users/matevz/g4snitch-2026-11/samplesA1.sh): the A1 sample matrix on uaf-4,
# one particle per event from the origin, eta 1.6-2.9, D128, field branch on, no step log: e-, pi+, gamma at pT 1, 5,
# 10, 20, 50, 100 GeV, 200 events each, mu- 100 events, in $G/A1/<name>_pt<pT>/. All 24 jobs in parallel.
# TTree output (G4SNITCH_FORMAT, default ttree here): the analysis scripts in ana/ read tree T as a TTree.
G=/ceph/users/matevz/g4snitch-2026-11
pgrep -u $USER -x cmsRun > /dev/null && { echo "cmsRun is running; not starting"; exit 1; }
mkdir -p $G/A1
for sp in e:11:200 pi:211:200 gamma:22:200 mu:13:100; do
  IFS=: read name pdg n <<< "$sp"
  for pt in 1 5 10 20 50 100; do
    d=$G/A1/${name}_pt$pt
    G4SNITCH_GUN=single G4SNITCH_PDG=$pdg G4SNITCH_PT=$pt G4SNITCH_NEVT=$n G4SNITCH_BFIELD=1 G4SNITCH_VERBOSE=0 G4SNITCH_FORMAT=${G4SNITCH_FORMAT:-ttree} \
      bash $G/run20.sh $d > $d.log 2>&1 &
  done
done
wait
for d in $G/A1/*/; do printf "%-14s %s  %s\n" $(basename $d) "$(tail -1 $d/info.txt)" "$(du -sh $d/G4Snitch.root | cut -f1)"; done
du -sh $G/A1
