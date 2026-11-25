// g4s_zsplit.C(file): inert and sensitive steps and their deposits by |z| of the step's end point, all particles.
#include <cmath>
#include <cstdio>
#include <vector>
#include "SimG4Core/HelpfulWatchers/interface/G4SnitchReader.h"
void g4s_zsplit(const char *fn) {
  G4SnitchFile f(fn);  // TTree or RNTuple file
  std::vector<G4S_ParticleSteps> *&s = f.s;
  std::vector<G4S_Particle> *&p = f.p;
  const double edges[] = {0, 300, 520, 1100, 5000, 1e9};  // cm
  const char *names[] = {"tracker+ETL |z|<3 m", "HGCal 3-5.2 m", "5.2-11 m", "11-50 m", "beyond 50 m"};
  const int nb = 5;
  long n[2][nb] = {}; double e[2][nb] = {};
  long n_path[nb] = {};
  for (long ev = 0; ev < f.n_events(); ++ev) {
    f.load(ev);
    for (size_t i = 0; i < s->size(); ++i) {
      for (int k = 0; k < 2; ++k) {
        const std::vector<G4S_Step> &v = k == 0 ? (*s)[i].m_steps_inert : (*s)[i].m_steps_sensitive;
        for (const G4S_Step &st : v) {
          double z = std::abs(st.m_x_end.fZ);
          int b = 0; while (b < nb - 1 && z >= edges[b + 1]) ++b;
          ++n[k][b]; e[k][b] += st.m_edep.m_total;
        }
      }
      for (const G4S_PathPoint &pp : (*s)[i].m_path) {
        double z = std::abs(pp.m_x.fZ);
        int b = 0; while (b < nb - 1 && z >= edges[b + 1]) ++b;
        ++n_path[b];
      }
    }
  }
  printf("%s, %lld events\n%-22s %12s %12s %12s %12s %10s\n", fn, (long long)f.n_events(), "|z| of step end", "inert steps", "inert MeV", "sens steps", "sens MeV", "path pts");
  for (int b = 0; b < nb; ++b)
    printf("%-22s %12ld %12.1f %12ld %12.1f %10ld\n", names[b], n[0][b], e[0][b], n[1][b], e[1][b], n_path[b]);
}
