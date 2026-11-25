#!/bin/bash
# rnt_check.sh (on uaf-4, after build.sh): the snitch's output at 1 and 4 threads. Re-simulates the 20 events of
# check-hgcal-nosnitch with their engine states restored, field at the path points on, in $G/rnt1 (1 thread) and
# $G/rnt4 (4 threads, one file written by the four snitches together). References, kept in $G from the
# multi-threading test of 2026-11-23 (made by the TTree code of that day): rr1-off, rr4-off (no snitch, 1 and 4
# threads) and rr1-on (1 thread, the snitch's TTree output). Checks:
# - the snitch is passive: the hits, SimTracks and SimVertices equal the no-snitch run's at the same threads;
# - the content, event by event (g4s_mt_check.py): rnt1 and rnt4 against rr1-on, and against each other;
# - each file's entries, clusters and event order, its tables' sizes, and the keys of table C that describe it;
# - the wall times and the file sizes.
G=/ceph/users/matevz/g4snitch-2026-11
C=$G/CMSSW_20_1_0_pre5/src/SimG4Core/HelpfulWatchers/g4s-test/check
cd $G
for p in $(pgrep -u $USER -x cmsRun); do
  case $(readlink /proc/$p/cwd) in $G/*) echo "cmsRun $p is running in $G; not starting"; exit 1 ;; esac
done
IN=$G/check-hgcal-nosnitch/step1.root
R() {  # R DIR THREADS
  rm -rf $G/$1
  G4SNITCH_INPUT=$IN G4SNITCH_NTHREADS=$2 G4SNITCH_BFIELD=1 G4SNITCH_VERBOSE=0 \
    timeout 1800 bash run20.sh $G/$1 rerunG4Snitch_D128.py > $G/$1.log 2>&1
}
R rnt1 1
R rnt4 4
for d in rr1-on rnt1 rnt4; do
  echo "$d: $(head -1 $d/info.txt | cut -c7-25) .. $(tail -1 $d/info.txt | cut -c5-) first event $(grep -m1 'Begin processing the 1st' $d/e | grep -o '[0-9:]*\.[0-9]* PST')"
  grep -m2 -i "exception\|fatal" $d/e
done
ls -la rr1-on/G4Snitch.root rnt1/*.root rnt4/*.root | awk '{print $5, $9}'
source env.sh > /dev/null 2>&1
H() { python3 $C/g4s_hits_check.py $1 $2 2>&1 | grep -v -i warning | tail -1; }
echo "=== the snitch is passive: no snitch against rnt1 (1 thread), rnt4 (4 threads)"
H rr1-off/rerun.root rnt1/rerun.root
H rr4-off/rerun.root rnt4/rerun.root
M() { python3 $C/g4s_mt_check.py "$@" 2>&1 | grep -v -i warning | tail -4; }
echo "=== content: rr1-on (1 thread, TTree, 2026-11-23) against rnt1"
M rr1-on/G4Snitch.root -- rnt1/G4Snitch.root
echo "=== content: rr1-on against rnt4 (4 threads)"
M rr1-on/G4Snitch.root -- rnt4/G4Snitch.root
echo "=== content: rnt1 against rnt4"
M rnt1/G4Snitch.root -- rnt4/G4Snitch.root
python3 - <<'EOF' 2>&1 | grep -v -i warning
import ROOT
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
for fn in ["rr1-on/G4Snitch.root", "rnt1/G4Snitch.root", "rnt4/G4Snitch.root"]:
    sf = SnitchFile(fn)
    C = sf.config()
    n = {k: sum(1 for _ in sf.rows(k)) for k in "PVD"}
    evs = [int(e.ev) for e in sf]
    what = f"TTree {len(sf)} entries" if sf.is_ttree else \
        f"RNTuple {len(sf)} entries in {sf.r.GetDescriptor().GetNClusters()} clusters"
    kinds = {k: ("TTree" if isinstance(sf.f.Get(k), ROOT.TTree) else "RNTuple") for k in "PVDC"}
    print(f"{fn}: {what}, ev {evs}; tables P {n['P']} V {n['V']} D {n['D']} ({', '.join(f'{k} {v}' for k, v in kinds.items())});"
          f" C: output_name {C.get('output_name')}, snitches {C.get('snitches')}, run_in_job {C.get('run_in_job')}")
    print(f"    format: {C.get('format')}")
EOF
