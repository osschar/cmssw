#!/bin/bash
# samplesT1.sh (copied to /ceph/users/matevz/g4snitch-2026-11/samplesT1.sh): the tracker-wide sample T1 on uaf-4 for
# the field-integral and material comparisons with mkFit. One particle per event from the origin, eta uniform in
# -3..3, phi uniform, D128 (tracker T35 = D121's), field branch on, no step log. mu- at pT 1, 2, 5, 10, 100 GeV and
# pi+ at pT 1, 10 GeV; each configuration in 4 jobs of different G4SNITCH_SEED (250 mu or 125 pi events each), in
# $G/T1/<name>_pt<pT>_s<seed>/. All 28 jobs in parallel.
G=/ceph/users/matevz/g4snitch-2026-11
pgrep -u $USER -x cmsRun > /dev/null && { echo "cmsRun is running; not starting"; exit 1; }
mkdir -p $G/T1
for spec in mu:13:250:1,2,5,10,100 pi:211:125:1,10; do
  IFS=: read name pdg n pts <<< "$spec"
  for pt in ${pts//,/ }; do
    for seed in 1 2 3 4; do
      d=$G/T1/${name}_pt${pt}_s$seed
      G4SNITCH_GUN=single G4SNITCH_PDG=$pdg G4SNITCH_PT=$pt G4SNITCH_ETAMIN=-3 G4SNITCH_ETAMAX=3 \
        G4SNITCH_NEVT=$n G4SNITCH_SEED=$((1000 * pt + seed)) G4SNITCH_BFIELD=1 G4SNITCH_VERBOSE=0 \
        bash $G/run20.sh $d > $d.log 2>&1 &
    done
  done
done
wait
for d in $G/T1/*/; do printf "%-16s %s  %s\n" $(basename $d) "$(tail -1 $d/info.txt)" "$(du -sh $d/G4Snitch.root | cut -f1)"; done
du -sh $G/T1
