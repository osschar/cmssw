#ifndef HelpfulWatchers_G4SnitchDataFormat_h

#include "ROOT/REveVector.hxx"

#include <algorithm>
#include <limits>

// Information about the event

struct G4S_Info {
  int    m_n_primaries {0};
  double m_min_time {0};
  double m_max_time {0};

  void reset() {
    m_n_primaries = 0;
    m_min_time = std::numeric_limits<double>::max();
    m_max_time = -m_min_time;
  }
  void update_min_time(double t) { m_min_time = std::min(m_min_time, t); }
  void update_max_time(double t) { m_max_time = std::max(m_max_time, t); }

  // G4S_Info() = default;
  // G4S_Info& operator=(const G4S_Info&) = default;
};

// A tracked G4 particle

struct G4S_Particle {
  using Vec4D = ROOT::Experimental::REveVector4D;
  Vec4D  m_x_beg;
  Vec4D  m_p_beg;
  Vec4D  m_x_end;
  Vec4D  m_p_end;
  double m_mass {0};
  int    m_pdg  {0};
  int    m_charge {0};
  int    m_parent {-1};
  int    m_daughters_begin {-1};
  int    m_daughters_end   {-1};
  int    m_n_daughters_skipped     {0};
  int    m_n_sub_daughters_skipped {0};
  int    m_g4_id       {-1};
  int    m_g4_level    {-1};
  bool   m_was_tracked {false};

  int n_daughters() const { return m_daughters_end - m_daughters_begin; }

  // G4S_Particle() = default;
  // G4S_Particle& operator=(const G4S_Particle&) = default;
};

// Structures for collecting energy deposits.
//
// For particles it possible to keep them in vector-form and still maintain
// the tree information by grouping all secondaries of a given paticle into
// sequential index range (because any particle is fully tracked from start to
// end and we get a chance to place all the secondaries in their rightful place.
//
// Now, we do have cuts on particles -- we do not store all of them. But energy
// deposits from untracked seconaries can be still relevant, especially if the
// cut on secondary energy is relatively high. G4Snitch::filter() tries to keep
// particles that can produce significant energy deposits later on (non-stable
// particles and postrons) ... but this selection might be shakey (and also produce
// more particles than we care about).
//
// We still get callbacks into G4Snitch for all begin/end track events and for
// all step events for all descendant particles, whether the snitch has decided to
// track them or not. So we can have another cut on the individual energy deposits
// from a step and store them along with the deposits from the last tracked particle.
//
// Clearly, we can not anticipate in advance how many such deposits there will be,
// so we need a different storage strategy to capture all of the deposits that we
// care about.
// It seems we need a vector of steps + some summary data for each tracked particle.
// These vectors are then put together into an event-level vector, parallel to
// the particle vector.

// Struct binding together types of E-deposits.

struct G4S_EDep {
  double m_total {0};
  double m_niel  {0}; // non-ionizing energy loss (for low energy ions only?)

  void set_edep(double t, double n) {
    m_total = t;
    m_niel  = n;
  }
  void add_edep(double t, double n) {
    m_total += t;
    m_niel  += n;
  }

  // G4S_EDep() = default;
  // G4S_EDep& operator=(const G4S_EDep&) = default;
};

// Individual G4Steps with EDep information.

struct G4S_Step {
  using Vec4D = ROOT::Experimental::REveVector4D;
  Vec4D    m_x_beg;
  Vec4D    m_x_end;
  G4S_EDep m_edep;
  // To add: process type at step-end
  // Also, do we need to filter by process type in the snitch?
  int      m_g4_pid;   // G4 particle id -- to distinguish this particle or untracked secondaries.
  int      m_g4_pv_id; // PhysicalVolume -- actually does not exist.
                       // TODO: figure out what we need here.

  // G4S_Step() = default;
  // G4S_Step& operator=(const G4S_Step&) = default;
};

// Step and EDep storage, one per Particle.

struct G4S_ParticleSteps {
  std::vector<G4S_Step> m_steps_sensitive;
  std::vector<G4S_Step> m_steps_inert;

  G4S_EDep m_edep_sum_all; // Summed-up E-deposits from all steps
  int      m_n_edeps_all {0};

  G4S_EDep m_edep_sum_sensitive; // Summed-up E-deposits from steps with E-dep below threshold
  int      m_n_edeps_sensitive {0};

  void add_edep(double t, double n, bool is_sensitive) {
    m_edep_sum_all.add_edep(t, n);
    ++m_n_edeps_all;
    if (is_sensitive) {
      m_edep_sum_sensitive.add_edep(t, n);
      ++m_n_edeps_sensitive;
    }
  }

  // G4S_ParticleSteps() = default;
  // G4S_ParticleSteps(G4S_ParticleSteps &&) = default;
  // G4S_ParticleSteps& operator=(const G4S_Particle&) = default;
};

// In G4Snitch, these are then registered into a branch of type std::vector<G4S_ParticleSteps>.
// There can be two vectors, one for sensitive and the other for inert material.

// In post-processing, one can then bubble up energy from the bottom and assign it to tracked
// particles --- as the total energy deposited by a given particle and all of its descendants.
// This probably most interesting for the energy deposited in the sensitive volumes, so one
// can then reason about the importance of a specific particle for the shower detection.

#endif
