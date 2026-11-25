# G4Snitch checks (CMSSW_20_1_0_pre5, D128)

The snitch writes RNTuples since 2026-11-24 (events T, tables P, V, D, C); older files hold TTrees. Everything here reads
both: the Python scripts through `SimG4Core/HelpfulWatchers/python/g4s_read.py` (`from SimG4Core.HelpfulWatchers.g4s_read
import SnitchFile`), the ROOT macros through `interface/G4SnitchReader.h` (class `G4SnitchFile`).

Run in the area's environment (`/ceph/users/matevz/g4snitch-2026-11/env.sh` on uaf-4); the ROOT macros need
`gSystem->Load("libRootG4Snitch.so"); gSystem->AddIncludePath("-I$CMSSW_BASE/src");` and are run with ACLiC (`+`).

- `check20.sh [cms|vertex|hgcal]`: the gun twice, DetIds on and off; logs, hits, DetIds, coverage (see its header).
- `rnt_check.sh`: the output at 1 and 4 threads (one file written by all threads) against the runs without the snitch
  and the TTree output of 2026-11-23 (see its header).
- `loop_check.sh [LOOPS]`: one event re-simulated LOOPS times in one job with a looper; each loop's hits and step log
  against the first (needs the area's patch to `SimG4Core/Application`).
- `g4s_primaries.py FILE`: per event the primaries and the deposits at z > 0 and z < 0.
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
