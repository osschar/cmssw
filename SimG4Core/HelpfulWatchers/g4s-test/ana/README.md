# Analyses of G4Snitch samples (the A steps of the HGCal shower-following plan)

Run in the area's environment on uaf-4 (`/ceph/users/matevz/g4snitch-2026-11/env.sh`).

- `samplesA1.sh`: the A1 matrix (one particle per event from the origin, eta 1.6-2.9, D128): e-, pi+, gamma at pT 1,
  5, 10, 20, 50, 100 GeV, 200 events each, mu- 100 events; in `/ceph/users/matevz/g4snitch-2026-11/A1/`, 90 GB
  (2026-11-21, `cfe4cc928b1`).
- `g4s_a2.py FILE...`: where the shower starts (hadron: first inelastic interaction; photon: conversion; electron:
  half of the energy radiated), in which region, at how many X0 and lambda_I from the origin, and the HGCal layers
  crossed before it.
- `g4s_front_material.py FILE [ZMIN ZMAX]`: the primary's material between two |z| per volume (for muons).
- `g4s_towers.py FIELDTABLE FILE...`: how much of a single particle's visible HGCal energy lands in the tower its
  track points to, on the pointing-tower grid (|eta| from 1.305; EM 0.0174 x 2pi/360, hadronic 0.087 x 2pi/72);
  pointing by the track's crossing of the HGCal face (RK4 through the field table) and by the vertex direction.
- `g4s_lowpt.py --field FIELDTABLE FILE...`: the follower's ceiling at low pT: visible HGCal energy in the 3x3 towers
  around the face crossing against the 3x3 along the particle's true path at each deposit's depth, and corridors of
  3, 5, 10 cm; who reaches HGCal and what stops the others (samples `A1lo/`: pi, mu, e at pT 0.3, 0.5, 2 GeV).
- `g4s_minbias.py`: what a min-bias interaction (`G4SNITCH_GUN=minbias`, Pythia8 CP5 at 14 TeV, HL-LHC vertex
  smearing; `G4SNITCH_SEED` per job) puts into HGCal: visible energy by entrant class (nearest ancestor born at
  |z| < 320 cm) and pT, and the charged energy outside the 3x3 around each entrant's face crossing (samples `MB/`,
  1000 events). `g4s_seedcheck.py`: equal seeds give equal events, different seeds different ones.

## Spinoffs for mkFit's tracker models

- `samplesT1.sh`: the tracker-wide sample T1, one particle per event from the origin, eta -3..3, D128 (tracker T35 =
  D121's), field branch on: mu- at pT 1, 2, 5, 10, 100 GeV (1000 events each) and pi+ at pT 1, 10 GeV (500 each), 4
  jobs per configuration with different `G4SNITCH_SEED`; in `/ceph/users/matevz/g4snitch-2026-11/T1/`, 5.4 GB
  (2026-11-22).
- `g4s_field.py OUT FILE...`: mkFit's parametrised field against the CMSSW field map, by the bending integral along
  the primary's true path between its first and last tracker hit (outgoing leg only): the constant field, track
  building's model (ParabolicMf constants, Bz only) and the final fit's (bFieldParams, with Br), each integrated
  along the path and as mkFit samples it per step (building: Bz at the step's start; final fit: Bz at the chord mid
  point plus the radial-field kick of `FinalFit.cc`, with 1, 2 or 4 sub-steps). Per track the relative error of the
  azimuthal bending, the pT bias of a position fit, and the residual left after it. `g4s_field_sum.py OUT...`: the
  tables per |eta| bin.
