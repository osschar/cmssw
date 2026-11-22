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
- `g4s_tkmat.py GEO OUT FILE...`: mkFit's tracker material against Geant4's, per step between consecutive hits of
  the primary: the true X0 and lambda_I from one hit's mid point to the next, against the material mkFit applies on
  arriving at the destination hit (the module's MediumProperties radl, and the (|z|, r) grid, times p/|p.n|). GEO is
  mkFit's geometry dumped by `../geo/dump_mkfit_geo.C` (Ceph `T1/mat/mkfit_geo_D121.root`). `g4s_tkmat_sum.py
  OUT...`: per |eta| bin and per arrival at a physical layer. `g4s_tkmat_vol.py ETAMIN ETAMAX FILE...`: the true X0
  between the first and the last hit per Geant4 volume (module parts grouped per sub-detector). `g4s_tkmat_once.py
  OUT...`: the per-track material with the stack value of PS/2S modules applied per sensor (as mkFit does) and once.
- `g4s_tgeo_cmp.py ENVELOPES DIFFUSE FILE...`: the TGeo material session's maps (`/baz/matevz/root-dev/geo-stuff/`,
  D123, tracker T35) against Geant4 along the true tracks: Geant4's X0 and lambda_I inside and outside that session's
  layer and disc envelopes (`../geo/envelopes_zpad05_D123.txt`, printed by `../geo/tgeo_envs.C`), its diffuse table
  `matdiffuse_rz_zpad05.root` integrated along the same paths, and Geant4's X0 per disc crossing against its disc
  budgets.
- `g4s_auxcount.py FILE...`: what auxiliary points every 10, 20, 30 degrees of field turn would add to the paths, per
  particle class, counted from existing files (field branch needed): the field turn of each path segment from the
  stored times, momenta and field; and how many large turns are field turns.
- `samplesT2.sh`: sample T2, pi+ at pT 1, 3, 10 GeV over |eta| < 3 (2000 each) and mu- at pT 2, 10 GeV over eta
  2.5-4.5 (1000 each), field branch on; in `/ceph/users/matevz/g4snitch-2026-11/T2/`, 9.6 GB (2026-11-22).
- `g4s_eloss.py GEO OUT FILE...`: mkFit's energy loss (applyMaterialEffects' mean Bethe-Bloch with the module's or
  the grid's bbxi, transcribed) against the true loss of the primary between consecutive hits. `g4s_eloss_sum.py
  OUT... [--ptmax P]`: per pT and |eta| bin, per arrival at a physical layer.
- `g4s_hadint.py OUT FILE...`: single charged hadrons in the tracker: inelastic interactions and decays (where,
  after how many hits, the leading charged daughter's momentum fraction, angle and hits), and the kinks from
  hadElastic and CoulombScat. `g4s_hadint_sum.py OUT...`: the tables.
- `g4s_tkmat_ref.py CSV MD FILE...`: the reference table of the true material per arrival at a physical layer and
  |eta| bin, with mkFit's values; `ref/tkmat_ref_D128.{csv,md}` from T1's muons and T2's forward muons.
- `g4s_disc_split.py GEO ENVELOPES FILE...`: where a forward pixel disc's material lies relative to its sensors: the
  true X0 between hits on adjacent discs, split into the source disc's envelope after its hit, the destination's
  before its hit, and elsewhere; "SUMS" lines add up over runs on several files.
