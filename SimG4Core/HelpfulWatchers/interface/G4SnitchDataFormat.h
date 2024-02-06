#ifndef HelpfulWatchers_G4SnitchDataFormat_h

#include "ROOT/REveVector.hxx"

#include <algorithm>
#include <limits>

struct G4S_Info {
  int m_n_primaries {0};
  double m_min_time {0};
  double m_max_time {0};

  void reset() {
    m_n_primaries = 0;
    m_min_time = std::numeric_limits<double>::max();
    m_max_time = -m_min_time;
  }
  void update_min_time(double t) { m_min_time = std::min(m_min_time, t); }
  void update_max_time(double t) { m_max_time = std::max(m_max_time, t); }

  G4S_Info() = default;
  G4S_Info& operator=(const G4S_Info&) = default;
};

struct G4S_Particle {
  using Vec4D = ROOT::Experimental::REveVector4D;
  Vec4D m_x_beg;
  Vec4D m_p_beg;
  Vec4D m_x_end;
  Vec4D m_p_end;
  double m_mass {0};
  int m_pdg {0};
  int m_charge {0};
  int m_parent {-1};
  int m_daughters_begin {-1};
  int m_daughters_end {-1};
  int m_n_daughters_skipped {0};
  int m_n_sub_daughters_skipped {0};
  int m_g4_id {-1};
  int m_g4_level {-1};
  bool m_was_tracked {false};

  int n_daughters() const { return m_daughters_end - m_daughters_begin; }

  G4S_Particle() = default;
  G4S_Particle& operator=(const G4S_Particle&) = default;
};

#endif
