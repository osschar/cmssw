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
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

class MagneticField;
class IdealMagneticFieldRecord;

class TFile;
class TTree;
namespace ROOT {
  class RNTupleFillContext;
  namespace Detail {
    class RRawPtrWriteEntry;
  }
}  // namespace ROOT

class BeginOfJob;
class BeginOfRun;
class BeginOfEvent;
class BeginOfTrack;
class G4Track;
class G4VProcess;
class G4Step;
class G4StepPoint;
class G4VSensitiveDetector;
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

// The tables of an output file: processes (tree "P"), volumes ("V") and sensitive detectors ("D").
// G4S_Particle::m_birth_proc and m_death_proc index the processes, G4S_Step::m_volume and G4S_InertSum::m_volume the
// volumes, and the volume table's column "sd" the sensitive detectors. A snitch writing its own file numbers one entry
// per Geant4 process and sensitive-detector object, keyed by pointer in the snitch. A table shared by the snitches of
// several threads (the RNTuple output) is keyed by content, since Geant4 has a process and a sensitive-detector object
// per thread: a process by (name, type, subtype), a sensitive detector by name, a volume by the hash of its
// touchable's history (the physical volumes are shared by the threads).
struct G4SnitchTables {
  std::vector<std::string> proc_names;
  std::vector<int> proc_types, proc_subtypes;
  std::vector<std::string> vol_paths, vol_materials;
  std::vector<double> vol_x0, vol_lambda, vol_density, vol_x, vol_y, vol_z;
  std::vector<int> vol_sd;
  std::vector<std::string> sd_names, sd_classes;
  // The keys by content, used only by a shared table.
  std::map<std::tuple<std::string, int, int>, int> proc_key;
  std::unordered_map<uint64_t, int> vol_key;
  std::map<std::string, int> sd_key;
};

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
  void write_tables(const G4SnitchTables &tab);
  void write_config_tree(bool shared_file);
  void reset_output_structs();
  void open_shared_rntuple();
  void close_shared_rntuple();

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
  // Volumes whose steps make no hits of their own: no DetId for steps under them. HF (HVQF): its hits come from CMS's
  // shower library and from fibre light; its readout box (OQUA): fibre-bundle and PMT hits are off by default.
  std::vector<std::string> m_no_det_id_volumes = {"HVQF", "OQUA"};

  // Thinned path of each tracked particle (G4S_ParticleSteps::m_path).
  bool m_output_path = true;
  double m_path_max_angle = 0.05;  // rad: store a point when the direction turned by more than this ...
  double m_path_max_length = 5.0;  // cm: ... or the path since the last point is longer than this
  double m_path_cos_max_angle = 1.0;
  // State of the path of the track being followed, in G4 units.
  double m_path_dir[3] = {0, 0, 0};
  double m_path_len_since = 0;
  // Material crossed by the current Geant4 track since its birth (every track, tracked or not; reset at BeginOfTrack).
  double m_mat_x0 = 0, m_mat_lambda = 0;
  void add_path_point(int id, const double *pos_mm, double t_ns, const double *mom_MeV, double e_MeV);
  // The field Geant4 used at each path point, in the separate branch "b" (off by default).
  bool m_output_path_bfield = false;
  // Kinetic energies (MeV): from track_ekin_cut at creation a stable particle is tracked as its own G4S_Particle (see
  // filter()); below last_state_ekin_cut its m_x_last / m_p_last are no longer updated.
  double m_track_ekin_cut = 10.0;
  double m_last_state_ekin_cut = 1.0;
  // The two rules of filter() before the energy cut: unstable particles other than neutrons are tracked at any
  // energy (track_unstable_always); positrons are cut on total, not kinetic, energy (positron_by_total_energy).
  bool m_track_unstable_always = true;
  bool m_positron_by_total_energy = true;
  // The output file. With several framework threads every thread has its own snitch, so output_per_thread inserts
  // "_t<Geant4 thread id>" before ".root" (G4Snitch_t0.root, ...); events are matched across files by branch "ev".
  std::string m_output_file = "G4Snitch.root";
  bool m_output_per_thread = false;
  // The format of the event data: "rntuple" (the default: RNTuple T in one file shared by the snitches of all
  // threads, written in parallel), "ttree" (tree T, a file per snitch), or "both" (the TTree files as with "ttree",
  // and the RNTuple file with "_rnt" inserted before ".root"; the same objects feed both). With an RNTuple the tables
  // are shared by all snitches of the job, in both outputs.
  std::string m_output_format = "rntuple";
  bool m_write_ttree = true, m_write_rntuple = false;
  std::shared_ptr<ROOT::RNTupleFillContext> m_rnt_fill;
  std::unique_ptr<ROOT::Detail::RRawPtrWriteEntry> m_rnt_entry;
  std::string m_output_name;
  int m_thread_id = -1;
  int m_ev_number = -1;  // the Geant4 event id, which CMS sets to the CMSSW event number; branch "ev" of tree T
  // The parameter set as given, written with the values in use into the configuration tree "C".
  std::string m_pset_dump;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> m_bfield_token;
  const MagneticField *m_bfield = nullptr;

  TFile *m_file = nullptr;
  TTree *m_tree = nullptr;

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

  // The tables (trees P, V, D): the snitch's own, or with an RNTuple the job's shared one, which m_tables_mutex guards.
  G4SnitchTables m_own_tables;
  G4SnitchTables *m_tables = &m_own_tables;
  std::mutex *m_tables_mutex = nullptr;
  // This snitch's index of each process object, volume and sensitive-detector object into the tables. A volume is a
  // placement, keyed by a hash of the (physical volume, copy number) pairs of its touchable's whole history.
  std::map<const G4VProcess *, int> m_proc_index;
  int process_index(const G4VProcess *proc);
  std::unordered_map<uint64_t, int> m_vol_index;
  std::vector<bool> m_vol_no_det_id;  // by volume index: a volume under one of m_no_det_id_volumes
  int volume_index(const G4StepPoint *sp);
  std::map<const G4VSensitiveDetector *, int> m_sd_index;
  int sd_index(const G4VSensitiveDetector *sd);

  // Per-event accumulators of the inert sums, keyed by (particle slot << 32 | volume index). Positions and times
  // are summed relative to the first step of the sum, which keeps the RMS precise far from the origin.
  struct InertAcc {
    int slot, vol, n;
    double e, niel, ref[4], s1[4], s2[6], s2t;  // s2: xx, yy, zz, xy, xz, yz
  };
  std::unordered_map<uint64_t, int> m_inert_acc_index;
  std::vector<InertAcc> m_inert_acc;
  void add_inert_step(int slot, const G4Step *step, double tot, double niel);
  void flush_inert_sums();

  // Output selection at the end of the event (the simulation is not touched): the subtree of the tracked particle
  // with this Geant4 track id (a primary or a SimTrack; 0 = all), and the steps and inert sums inside a cone from the
  // origin (eta, phi, dR; empty = all; the particle tree and the paths stay whole).
  int m_select_track_id = 0;
  std::vector<double> m_select_cone;
  void select_output();

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
