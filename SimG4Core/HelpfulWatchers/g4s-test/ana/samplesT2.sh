#!/bin/bash
# samplesT2.sh (copied to /ceph/users/matevz/g4snitch-2026-11/samplesT2.sh): samples T2 on uaf-4, one particle per
# event from the origin, D128, field branch on, no step log, in $G/T2/<name>_pt<pT>_s<seed>/:
# - pi+ at pT 1, 3, 10 GeV, eta uniform in -3..3, 8 jobs of 250 events per pT (2000): nuclear interactions and kinks in
#   the tracker;
# - mu- at pT 2, 10 GeV, eta uniform in 2.5..4.5, 4 jobs of 250 events per pT (1000): the forward pixel discs (TFPX
#   inner rings, all of TEPX), beyond T1's |eta| < 3.
# Seeds distinct from T1's. All 32 jobs in parallel.
# TTree output (G4SNITCH_FORMAT, default ttree here): the analysis scripts in ana/ read tree T as a TTree.
G=/ceph/users/matevz/g4snitch-2026-11
pgrep -u $USER -x cmsRun > /dev/null && { echo "cmsRun is running; not starting"; exit 1; }
mkdir -p $G/T2
run() {  # name pdg pt etamin etamax nevt seed
  local d=$G/T2/$1_pt$3_s$7
  G4SNITCH_GUN=single G4SNITCH_PDG=$2 G4SNITCH_PT=$3 G4SNITCH_ETAMIN=$4 G4SNITCH_ETAMAX=$5 G4SNITCH_NEVT=$6 \
    G4SNITCH_SEED=$7 G4SNITCH_BFIELD=1 G4SNITCH_VERBOSE=0 G4SNITCH_FORMAT=${G4SNITCH_FORMAT:-ttree} bash $G/run20.sh $d > $d.log 2>&1
}
for pt in 1 3 10; do
  for i in 1 2 3 4 5 6 7 8; do run pi 211 $pt -3 3 250 $((50000 + 1000 * pt + i)) & done
done
for pt in 2 10; do
  for i in 1 2 3 4; do run mufwd 13 $pt 2.5 4.5 250 $((70000 + 1000 * pt + i)) & done
done
wait
for d in $G/T2/*/; do printf "%-18s %s  %s\n" $(basename $d) "$(tail -1 $d/info.txt)" "$(du -sh $d/G4Snitch.root | cut -f1)"; done
du -sh $G/T2
