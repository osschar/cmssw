# G4Snitch checks (CMSSW_20_1_0_pre5, D128)

The Python scripts read either format through `g4s_read.py`; the ROOT macros (`.C`) read tree T as a TTree, so they need
output written with `G4SNITCH_FORMAT=ttree` (the default output is an RNTuple since 2026-11-24).

Run in the area's environment (`/ceph/users/matevz/g4snitch-2026-11/env.sh` on uaf-4); the ROOT macros need
`gSystem->Load("libRootG4Snitch.so"); gSystem->AddIncludePath("-I$CMSSW_BASE/src");` and are run with ACLiC (`+`).

- `check20.sh [cms|vertex|hgcal]`: the gun twice, DetIds on and off; logs, hits, DetIds, coverage (see its header).
- `rnt_check.sh`: the RNTuple output (rntuple at 4 threads, both at 1 and 4) against the TTree output and the runs
  without the snitch (see its header).
- `loop_check.sh [LOOPS]`: one event re-simulated LOOPS times in one job with a looper; each loop's hits and step log
  against the first (needs the area's patch to `SimG4Core/Application`).
- `g4s_read.py`: the Python reader the scripts here share; a snitch file's events T as a TTree or an RNTuple, read the
  same way (`SnitchFile(fn)`, then `for ev in sf: ev.ev, ev.p, ev.s, ...`).
- `g4s_mt_check.py A... -- B...`: two outputs of the same events compared event by event, matched by `ev`, whatever
  the format and the number of files.
- `g4s_subtree_check.py`: a re-run that writes one particle's subtree against the same subtree in a full file.
- `g4s_hits_check.py A/step1.root B/step1.root [B/G4Snitch.root]`: hit collections equal; snitch DetIds among the
  event's hits, by sensitive detector (HCAL by part: its ids are test-numbering indices).
- `g4s_subdets.py G4Snitch.root`: sensitive steps by detector class and subdetector, with their |eta| range.
- `g4s_attrib.C(file)`: steps filed under a particle with a smaller Geant4 track id than its own (misattributed).
- `g4s_aggr_check.C(off, on)`: inert sums equal the descendants' inert steps, per particle (aggregation off vs on).
- `g4s_sums.C(file)`: the inert sums by material and the volumes with the most energy.
- `g4s_zsplit.C(file)`: steps and path points by |z|.
