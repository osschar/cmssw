#ifndef HelpfulWatchers_G4SnitchDataFormat_h
#define HelpfulWatchers_G4SnitchDataFormat_h

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

// A particle's life cycle, positions and times in cm and ns, four-momenta in GeV:
// - birth: where and when it was created, and its four-momentum then (Geant4's creation vertex);
// - last:  its state after the last step that left it more than 1 MeV of kinetic energy (E_final_cut in G4Snitch.cc),
//          kept so that a drawing has an end direction; equal to birth for a particle killed in its first step;
// - death: where and when Geant4 ended the track, and its four-momentum then. The kinetic energy is zero when the
//          particle stopped or a process absorbed it (conversion, photoelectric, annihilation, inelastic
//          interaction, decay at rest). It is not zero for a decay in flight, when the particle left the world, or
//          when CMS's SteppingAction killed it (time window, dead region, energy threshold of a volume, beyond the
//          central |z|, too many steps, low energy in the beam pipe); then the death process is the last step's,
//          for example Transportation or hadElastic. D128 test, 20 events: 283 of 67259 tracked particles die
//          with kinetic energy, 279 of them by decay in flight.
// m_birth_proc, m_death_proc: indices into the file's process table, tree "P" (idx, name, type, subtype); -1 for
// none: the birth of a primary, or a particle never tracked (m_was_tracked false, death members zero).
struct G4S_Particle {
  using Vec4D = ROOT::Experimental::REveVector4D;
  Vec4D  m_x_birth;
  Vec4D  m_p_birth;
  Vec4D  m_x_last;
  Vec4D  m_p_last;
  Vec4D  m_x_death;
  Vec4D  m_p_death;
  int    m_birth_proc {-1};
  int    m_death_proc {-1};
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
  Vec4D    m_p_end;  // momentum (GeV) and total energy at the step's end point, of the stepping track (see m_g4_pid)
  G4S_EDep m_edep;
  // To add: process type at step-end
  // Also, do we need to filter by process type in the snitch?
  int      m_g4_pid;   // G4 particle id -- to distinguish this particle or untracked secondaries.
  int      m_g4_pv_id; // PhysicalVolume -- actually does not exist.
                       // TODO: figure out what we need here.

  // G4S_Step() = default;
  // G4S_Step& operator=(const G4S_Step&) = default;
};

// A point on the path of a tracked particle: position (cm) and time (ns), momentum (GeV) and total energy.

struct G4S_PathPoint {
  using Vec4D = ROOT::Experimental::REveVector4D;
  Vec4D m_x;
  Vec4D m_p;
};

// The magnetic field at each point of a particle's path (tesla), parallel to G4S_ParticleSteps::m_path: the CMSSW
// MagneticField that Geant4's field wraps, asked at the point (Geant4 itself may use the value at a point within
// its caching tolerance). Written to the separate branch "b" only when the G4Snitch parameter
// output_path_bfield is on, so files without it keep the format of the other branches.

struct G4S_PathBField {
  using Vec3F = ROOT::Experimental::REveVectorF;
  std::vector<Vec3F> m_b;
};

// Step and EDep storage, one per Particle.
//
// m_path is the thinned path of the particle itself (never of its untracked descendants): its start, its end, and
// the end points of the steps after which its direction has turned by more than path_max_angle, or its path has
// grown by more than path_max_length, since the last stored point (G4Snitch parameters). A step is stored when its
// energy deposit is at least sensitive_step_ecut or inert_step_ecut, in MeV (0 and 1e-5 MeV = 10 eV in the test
// configurations, so every sensitive step, most with no deposit, and inert steps from 10 eV), so the steps do not describe a
// trajectory; the path does. Deposits (G4S_EDep) are in MeV, unlike the GeV of the momenta.

struct G4S_ParticleSteps {
  std::vector<G4S_Step> m_steps_sensitive;
  std::vector<G4S_Step> m_steps_inert;
  std::vector<G4S_PathPoint> m_path;

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
