# G4Snitch_cfi.py: the G4Snitch watcher for g4SimHits, configured from the environment (the test configurations in
# g4s-test/ share it). Use: process.g4SimHits.Watchers = cms.VPSet(g4snitch_watcher()), or none() for a run without.
#   G4SNITCH_OFF=1            no snitch at all (for the passivity check)
#   G4SNITCH_VERBOSE=0        no step-by-step log (production; the checks compare these logs)
#   G4SNITCH_BFIELD=1         the field at each path point, branch "b"
#   G4SNITCH_ZEROSTEPS=1      keep the sensitive steps with no deposit
#   G4SNITCH_AGGR=0           store the untracked descendants' inert steps one by one instead of their sums per volume
#   G4SNITCH_DETID=0          do not ask the sensitive detectors for DetIds
#   G4SNITCH_SELECT_ID=N      write only the subtree of the tracked particle with Geant4 track id N
#   G4SNITCH_SELECT_CONE=eta,phi,dR   write only the steps and inert sums inside the cone
#   G4SNITCH_TRACK_EKIN_CUT=MeV       kinetic energy from which a stable particle is tracked (default 10)
#   G4SNITCH_LAST_EKIN_CUT=MeV        kinetic energy below which a particle's last state is no longer updated (default 1)
# Every value in use, with the code version and these environment variables, goes into the output's tree "C".
import os
import FWCore.ParameterSet.Config as cms

def g4snitch_watcher(env=os.environ):
    w = cms.PSet(
        type = cms.string('G4Snitch'),
        verbose = cms.untracked.bool(True),
        verbose_stack_level = cms.untracked.bool(False),
        verbose_transport = cms.untracked.bool(False),
        verbose_skip = cms.untracked.bool(False),
        verbose_skip_with_ids = cms.untracked.bool(False),

        output_sensitive_steps = cms.untracked.bool(True),
        output_inert_steps = cms.untracked.bool(True),
        sensitive_step_ecut = cms.untracked.double(0.0),  # MeV
        output_sensitive_zero_steps = cms.untracked.bool(False),  # sensitive steps with no deposit, as path marks
        inert_step_ecut = cms.untracked.double(10e-6),  # MeV: 10 eV
        inert_aggregate_descendants = cms.untracked.bool(True),  # untracked descendants' inert steps summed per volume
        sensitive_det_id = cms.untracked.bool(True),  # the DetId of each stored sensitive step
        no_det_id_volumes = cms.untracked.vstring("HVQF", "OQUA"),  # HF and its readout box: no hits from steps there

        output_path = cms.untracked.bool(True),
        path_max_angle = cms.untracked.double(0.05),  # rad
        path_max_length = cms.untracked.double(5.0),  # cm
        output_path_bfield = cms.untracked.bool(False),  # the field at each path point, in branch "b"

        track_ekin_cut = cms.untracked.double(10.0),  # MeV: a stable particle from this kinetic energy is tracked
        last_state_ekin_cut = cms.untracked.double(1.0),  # MeV: below this the last state is no longer updated
    )
    if env.get("G4SNITCH_BFIELD", "0") == "1":
        w.output_path_bfield = True
    if env.get("G4SNITCH_ZEROSTEPS", "0") == "1":
        w.output_sensitive_zero_steps = True
    if env.get("G4SNITCH_AGGR", "1") == "0":
        w.inert_aggregate_descendants = False
    if env.get("G4SNITCH_VERBOSE", "1") == "0":
        w.verbose = False
    if env.get("G4SNITCH_DETID", "1") == "0":
        w.sensitive_det_id = False
    if env.get("G4SNITCH_SELECT_ID"):
        w.select_track_id = cms.untracked.int32(int(env["G4SNITCH_SELECT_ID"]))
    if env.get("G4SNITCH_TRACK_EKIN_CUT"):
        w.track_ekin_cut = float(env["G4SNITCH_TRACK_EKIN_CUT"])
    if env.get("G4SNITCH_LAST_EKIN_CUT"):
        w.last_state_ekin_cut = float(env["G4SNITCH_LAST_EKIN_CUT"])
    if env.get("G4SNITCH_SELECT_CONE"):
        w.select_cone = cms.untracked.vdouble([float(x) for x in env["G4SNITCH_SELECT_CONE"].split(",")])
    return w

def g4snitch_watchers(env=os.environ):
    """The VPSet for g4SimHits.Watchers: the snitch, or nothing with G4SNITCH_OFF=1."""
    if env.get("G4SNITCH_OFF", "0") == "1":
        return cms.VPSet()
    return cms.VPSet(g4snitch_watcher(env))
