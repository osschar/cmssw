# Analyses of G4Snitch samples (the A steps of the HGCal shower-following plan)

Run in the area's environment on uaf-4 (`/ceph/users/matevz/g4snitch-2026-11/env.sh`).

- `samplesA1.sh`: the A1 matrix (one particle per event from the origin, eta 1.6-2.9, D128): e-, pi+, gamma at pT 1,
  5, 10, 20, 50, 100 GeV, 200 events each, mu- 100 events; in `/ceph/users/matevz/g4snitch-2026-11/A1/`, 90 GB
  (2026-11-21, `cfe4cc928b1`).
- `g4s_a2.py FILE...`: where the shower starts (hadron: first inelastic interaction; photon: conversion; electron:
  half of the energy radiated), in which region, at how many X0 and lambda_I from the origin, and the HGCal layers
  crossed before it.
- `g4s_front_material.py FILE [ZMIN ZMAX]`: the primary's material between two |z| per volume (for muons).
