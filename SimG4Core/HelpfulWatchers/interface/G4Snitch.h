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

class TFile;
class TTree;

class BeginOfJob;
class BeginOfRun;
class BeginOfEvent;
class BeginOfTrack;
class G4Track;
class G4Step;
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

  bool m_verbose;
  bool m_verbose_stack_level;
  bool m_verbose_transport;
  bool m_verbose_skip;
  bool m_verbose_skip_with_ids;

  bool m_output_sensitive_edeps = true;
  bool m_output_inert_edeps     = true;

  TFile *m_file;
  TTree *m_tree;
  std::shared_ptr<G4S_Info> m_info;
  std::shared_ptr<std::vector<G4S_Particle>> m_part_vec;
  G4S_Particle& particle(int i) { return (*m_part_vec)[i]; }

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

  int m_last_g4_id = -1;

  std::vector<int> m_stack; // stack of g4 ids
  int stack_level() const { return (int) m_stack.size() - 1; }
  int stack_top() const { return m_stack.back(); }
};

#undef UPDATE
#undef OBSERVES

#endif
