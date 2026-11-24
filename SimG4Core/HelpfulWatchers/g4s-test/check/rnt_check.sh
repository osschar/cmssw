#!/bin/bash
# rnt_check.sh (on uaf-4, after build.sh and after the restored-state runs rr1-off, rr1-on, rr4-off, rr4-on of the
# multi-threading test): checks of the RNTuple output (output_format rntuple and both). Re-simulates the 20 events of
# check-hgcal-nosnitch with their engine states restored, field at the path points on, in $G:
#   rnt4   4 threads, output_format rntuple: one G4Snitch.root for the job;
#   both1  1 thread, output_format both: G4Snitch.root (TTree) and G4Snitch_rnt.root;
#   both4  4 threads, output_format both: G4Snitch_t<thread>.root (TTrees) and G4Snitch_rnt.root.
# Checks:
# - the snitch is passive: each run's hits, SimTracks and SimVertices equal the no-snitch run's at the same threads;
# - the content, event by event (g4s_mt_check.py): the RNTuple of rnt4 against rr1-on's TTree (1 thread, the snitch's
#   own tables), and in both1 and both4 the RNTuple against the TTrees written by the same job;
# - the RNTuple's entries, clusters and event order, the tables' sizes, and the keys of tree C that describe the file;
# - the wall times and the file sizes.
G=/ceph/users/matevz/g4snitch-2026-11
C=$G/CMSSW_20_1_0_pre5/src/SimG4Core/HelpfulWatchers/g4s-test/check
cd $G
for p in $(pgrep -u $USER -x cmsRun); do
  case $(readlink /proc/$p/cwd) in $G/*) echo "cmsRun $p is running in $G; not starting"; exit 1 ;; esac
done
IN=$G/check-hgcal-nosnitch/step1.root
R() {  # R DIR THREADS FORMAT
  rm -rf $G/$1
  G4SNITCH_INPUT=$IN G4SNITCH_NTHREADS=$2 G4SNITCH_FORMAT=$3 G4SNITCH_BFIELD=1 G4SNITCH_VERBOSE=0 \
    timeout 1800 bash run20.sh $G/$1 rerunG4Snitch_D128.py > $G/$1.log 2>&1
}
R rnt4 4 rntuple
R both1 1 both
R both4 4 both
for d in rr1-on rr4-on rnt4 both1 both4; do
  echo "$d: $(head -1 $d/info.txt | cut -c7-25) .. $(tail -1 $d/info.txt | cut -c5-) first event $(grep -m1 'Begin processing the 1st' $d/e | grep -o '[0-9:]*\.[0-9]* PST')"
  grep -m2 -i "exception\|fatal" $d/e
done
ls -la rr1-on/G4Snitch.root rr4-on/G4Snitch_t*.root rnt4/*.root both1/*.root both4/*.root | awk '{print $5, $9}'
source env.sh > /dev/null 2>&1
H() { python3 $C/g4s_hits_check.py $1 $2 2>&1 | grep -v -i warning | tail -1; }
echo "=== the snitch is passive: no snitch against rnt4 (4 threads), both1 (1 thread), both4 (4 threads)"
H rr4-off/rerun.root rnt4/rerun.root
H rr1-off/rerun.root both1/rerun.root
H rr4-off/rerun.root both4/rerun.root
M() { python3 $C/g4s_mt_check.py "$@" 2>&1 | grep -v -i warning | tail -4; }
echo "=== content: rr1-on (1 thread, TTree) against rnt4 (4 threads, RNTuple)"
M rr1-on/G4Snitch.root -- rnt4/G4Snitch.root
echo "=== content: both1, TTree against RNTuple"
M both1/G4Snitch.root -- both1/G4Snitch_rnt.root
echo "=== content: both4, the 4 TTree files against the RNTuple"
M both4/G4Snitch_t*.root -- both4/G4Snitch_rnt.root
python3 - <<'EOF' 2>&1 | grep -v -i warning
import ROOT
ROOT.gSystem.Load("libRootG4Snitch.so")
for fn in ["rr1-on/G4Snitch.root", "rnt4/G4Snitch.root", "both1/G4Snitch.root", "both1/G4Snitch_rnt.root",
           "both4/G4Snitch_t0.root", "both4/G4Snitch_rnt.root"]:
    f = ROOT.TFile.Open(fn)
    C = {str(x.key): str(x.value) for x in f.Get("C")}
    n = {k: f.Get(k).GetEntries() for k in "PVD"}
    t = f.Get("T")
    if isinstance(t, ROOT.TTree):
        what = f"TTree {t.GetEntries()} entries, ev {[int(e.ev) for e in t]}"
    else:
        r = ROOT.RNTupleReader.Open("T", fn)
        e = r.CreateEntry()
        evs = []
        for i in range(r.GetNEntries()):
            r.LoadEntry(i, e)
            evs.append(int(e["ev"]))
        what = f"RNTuple {r.GetNEntries()} entries in {r.GetDescriptor().GetNClusters()} clusters, ev {evs}"
    print(f"{fn}: {what}; tables P {n['P']} V {n['V']} D {n['D']}; C: output_format {C.get('output_format')}, "
          f"output_name {C.get('output_name')}, thread_id {C.get('thread_id')}, snitches {C.get('snitches')}")
    print(f"    format: {C.get('format')}")
    print(f"    tables: {C.get('tables')}")
EOF
