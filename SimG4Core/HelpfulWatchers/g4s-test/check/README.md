# G4Snitch checks (CMSSW_20_1_0_pre5, D128)

Run in the area's environment (`/ceph/users/matevz/g4snitch-2026-11/env.sh` on uaf-4); the ROOT macros need
`gSystem->Load("libRootG4Snitch.so"); gSystem->AddIncludePath("-I$CMSSW_BASE/src");` and are run with ACLiC (`+`).

- `check20.sh [cms|vertex|hgcal]`: the gun twice, DetIds on and off; logs, hits, DetIds, coverage (see its header).
- `g4s_hits_check.py A/step1.root B/step1.root [B/G4Snitch.root]`: hit collections equal; snitch DetIds among the
  event's hits, by sensitive detector (HCAL by part: its ids are test-numbering indices).
- `g4s_subdets.py G4Snitch.root`: sensitive steps by detector class and subdetector, with their |eta| range.
- `g4s_attrib.C(file)`: steps filed under a particle with a smaller Geant4 track id than its own (misattributed).
- `g4s_aggr_check.C(off, on)`: inert sums equal the descendants' inert steps, per particle (aggregation off vs on).
- `g4s_sums.C(file)`: the inert sums by material and the volumes with the most energy.
- `g4s_zsplit.C(file)`: steps and path points by |z|.
