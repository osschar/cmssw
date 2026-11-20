#ifndef HelpfulWatchers_G4Snitch_h
#define HelpfulWatchers_G4Snitch_h
// -*- C++ -*-
//
// Package:     HelpfulWatchers
// Class  :     SimTracer
//
/**\class SimTracer SimTracer.h SimG4Core/HelpfulWatchers/interface/SimTracer.h

 Description: Prints a message for each Oscar signal

 Usage:
    <usage>

*/
//
// Original Author:
//         Created:  Tue Nov 22 16:41:33 EST 2005
//

// user include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimG4Core/Watcher/interface/SimWatcher.h"
#include "SimG4Core/Notification/interface/Observer.h"

#include "SimG4Core/HelpfulWatchers/interface/G4SnitchDataFormat.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

class MagneticField;
class IdealMagneticFieldRecord;

class TFile;
class TTree;

class BeginOfJob;
class BeginOfRun;
class BeginOfEvent;
class BeginOfTrack;
class G4Track;
class G4VProcess;
class G4Step;
class G4StepPoint;
class EndOfTrack;
class EndOfEvent;
class EndOfRun;
class EndOfJob;

//==============================================================================

/* On the output format

  What kind of output structure do we want and how much post-processing are we
  willing to do in the loader / REve / RenderCore?
  Also, do we just want kinematics tree or also energy deposits -- Marco?


*/

//==============================================================================

#define OBSERVES(type) public Observer<const type *>
#define UPDATE(type) void update(const type *) override

class G4Snitch : public SimWatcher,
                 OBSERVES(BeginOfJob),
                 OBSERVES(BeginOfRun),
                 OBSERVES(BeginOfEvent),
                 OBSERVES(BeginOfTrack),
                 OBSERVES(G4Step),
                 OBSERVES(EndOfTrack),
                 OBSERVES(EndOfEvent),
                 OBSERVES(EndOfRun),
                 OBSERVES(EndOfJob) {
public:
  G4Snitch(const edm::ParameterSet &pSet);

  // The CMSSW magnetic field, read at begin of run, for output_path_bfield.
  void registerConsumes(edm::ConsumesCollector) override;
  void beginRun(edm::EventSetup const &) override;
  virtual ~G4Snitch();

  UPDATE(BeginOfJob);
  UPDATE(BeginOfRun);
  UPDATE(BeginOfEvent);
  UPDATE(BeginOfTrack);
  UPDATE(G4Step);
  UPDATE(EndOfTrack);
  UPDATE(EndOfEvent);
  UPDATE(EndOfRun);
  UPDATE(EndOfJob);

private:
  void open_file_tree();
  void write_tree_close_file();
  void reset_output_structs();

  bool filter(const G4Track *t) const;

  float expected_event_progress() const;

  bool m_verbose = false;
  bool m_verbose_stack_level = false;
  bool m_verbose_transport = false;
  bool m_verbose_skip = false;
  bool m_verbose_skip_with_ids = false;

  bool m_output_sensitive_steps = true;
  bool m_output_inert_steps     = true;
  double m_sensitive_step_ecut = 0.0;  // MeV
  // Store sensitive steps with no deposit too (a particle's crossing of a sensitive volume, a potential path mark).
  bool m_output_sensitive_zero_steps = false;
  double m_inert_step_ecut = 10e-6; // MeV (Geant4 units): 10 eV
  // Sum the untracked descendants' inert steps per tracked particle and volume (G4S_ParticleSteps::m_inert_sums).
  bool m_inert_aggregate_descendants = true;
  // The DetId of each stored sensitive step, from the step's CMS sensitive detector (G4S_Step::m_det_id).
  bool m_sensitive_det_id = true;

  // Thinned path of each tracked particle (G4S_ParticleSteps::m_path).
  bool m_output_path = true;
  double m_path_max_angle = 0.05;  // rad: store a point when the direction turned by more than this ...
  double m_path_max_length = 5.0;  // cm: ... or the path since the last point is longer than this
  double m_path_cos_max_angle = 1.0;
  // State of the path of the track being followed, in G4 units.
  double m_path_dir[3] = {0, 0, 0};
  double m_path_len_since = 0;
  void add_path_point(int id, const double *pos_mm, double t_ns, const double *mom_MeV, double e_MeV);
  // The field Geant4 used at each path point, in the separate branch "b" (off by default).
  bool m_output_path_bfield = false;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> m_bfield_token;
  const MagneticField *m_bfield = nullptr;

  TFile *m_file;
  TTree *m_tree;

  std::shared_ptr<G4S_Info> m_info;
  std::shared_ptr<std::vector<G4S_Particle>> m_part_vec;
  std::shared_ptr<std::vector<G4S_ParticleSteps>> m_psteps_vec;
  std::shared_ptr<std::vector<G4S_PathBField>> m_pbfield_vec;  // only with m_output_path_bfield

  G4S_Particle& particle(int i) { return (*m_part_vec)[i]; }
  G4S_ParticleSteps& particle_steps(int i) { return (*m_psteps_vec)[i]; }
  // The per-particle vectors parallel to m_part_vec: grow, reserve and clear them together.
  void steps_emplace_back();
  void steps_resize(int n);
  void steps_reserve(int n);
  void steps_clear();

  int m_event = 0;
  int m_step_n = -1;
  int m_non_tracking_gid = -1;
  int m_tracks_skipped, m_steps_skipped, m_daughters_skipped;
  int m_event_tracks_accepted, m_event_tracks_skipped;
  int m_total_tracks_accepted, m_total_tracks_skipped;
  bool m_tracking = true;
  bool m_starting_new_event = false; // for primary initalization

  int m_vec_size = 0, m_vec_capacity = 0;
  int m_num_accepted_daugters_for_track = -1;
  int m_num_total_daugters_for_track = -1;
  int m_g4_num_primaries = -1;
  int m_g4_id_current_primary = -1; // as it comes from N_prim downwards -- primaries can be filtered!
  int m_id_current_primary = -1; // sequential up from 1
  int m_id; // current particle vector id
  int m_primary_n_accepted_tracks = -1; // needed to estimate kine tree size for vec.reserve()

  std::map<const G4Track*, int> m_gtp2vid;

  // The process table of the file (tree "P"): G4S_Particle::m_birth_proc and m_death_proc index it.
  std::map<const G4VProcess*, int> m_proc_index;
  std::vector<std::string> m_proc_names;
  std::vector<int> m_proc_types, m_proc_subtypes;
  int process_index(const G4VProcess *proc);

  // The volume table of the file (tree "V"): G4S_InertSum::m_volume indexes it. A volume is a placement, keyed by a
  // hash of the (physical volume, copy number) pairs of its touchable's whole history.
  std::unordered_map<uint64_t, int> m_vol_index;
  std::vector<std::string> m_vol_paths, m_vol_materials;
  std::vector<double> m_vol_x0, m_vol_density, m_vol_x, m_vol_y, m_vol_z;
  int volume_index(const G4StepPoint *sp);

  // Per-event accumulators of the inert sums, keyed by (particle slot << 32 | volume index). Positions and times
  // are summed relative to the first step of the sum, which keeps the RMS precise far from the origin.
  struct InertAcc {
    int slot, vol, n;
    double e, niel, ref[4], sw[4], sw2[4];
  };
  std::unordered_map<uint64_t, int> m_inert_acc_index;
  std::vector<InertAcc> m_inert_acc;
  void add_inert_step(int slot, const G4Step *step, double tot, double niel);
  void flush_inert_sums();

  int m_last_g4_id = -1;

  std::vector<int> m_stack; // stack of g4 ids
  // The particle-vector slot of each stack entry; -1 for an untracked one. While tracking is on, every entry below
  // the top is tracked, so a track that pauses tracking files its steps under its parent's slot.
  std::vector<int> m_stack_vid;
  int stack_level() const { return (int) m_stack.size() - 1; }
  int stack_top() const { return m_stack.back(); }
};

#undef UPDATE
#undef OBSERVES

#endif
