// mu_count2.C: muons per min-bias interaction whose own path crosses the HGCal face (|z| = 320 cm) at
// 1.5 < |eta| < 3.045, after being born at |z| < 320; by pT at the crossing; and their own visible HGCal energy
// (sensitive steps of the muon itself, DetId detector 8, 9, 10).
#include <cmath>
#include <cstdio>
#include <vector>
#include <glob.h>
#include <string>
#include "SimG4Core/HelpfulWatchers/interface/G4SnitchReader.h"
// The files of the min-bias sample, TTree or RNTuple (G4SnitchReader.h reads both).
static std::vector<std::string> mb_files() {
  glob_t g;
  glob("/ceph/users/matevz/g4snitch-2026-11/MB/job_*/G4Snitch.root", 0, nullptr, &g);
  std::vector<std::string> v(g.gl_pathv, g.gl_pathv + g.gl_pathc);
  globfree(&g);
  return v;
}
void mu_count2() {
  long n_ev = 0; double n_mu[4] = {0, 0, 0, 0}, e_mu[4] = {0, 0, 0, 0}, n_other_fwd = 0;
  for (const auto &fn : mb_files()) {
    G4SnitchFile f(fn);
    std::vector<G4S_Particle> *&p = f.p;
    std::vector<G4S_ParticleSteps> *&s = f.s;
    for (long e = 0; e < f.n_events(); ++e, ++n_ev) {
      f.load(e);
      for (size_t i = 0; i < p->size(); ++i) {
        const G4S_Particle &q = (*p)[i];
        if (std::abs(q.m_pdg) != 13 || !q.m_was_tracked || std::abs(q.m_x_birth.fZ) >= 320) continue;
        const auto &path = (*s)[i].m_path;
        int b = -1; double pt = 0;
        for (size_t j = 1; j < path.size(); ++j) {
          double z0 = std::abs(path[j - 1].m_x.fZ), z1 = std::abs(path[j].m_x.fZ);
          if (z0 < 320 && z1 >= 320) {
            double w = (320 - z0) / (z1 - z0);
            double x = path[j - 1].m_x.fX + w * (path[j].m_x.fX - path[j - 1].m_x.fX);
            double y = path[j - 1].m_x.fY + w * (path[j].m_x.fY - path[j - 1].m_x.fY);
            double eta = std::asinh(320 / std::hypot(x, y));
            pt = std::hypot(path[j].m_p.fX, path[j].m_p.fY);
            if (eta > 1.5 && eta < 3.045) b = pt < 0.5 ? 0 : pt < 1 ? 1 : pt < 2 ? 2 : 3;
            else ++n_other_fwd;
            break;
          }
        }
        if (b < 0) continue;
        n_mu[b] += 1;
        for (const G4S_Step &st : (*s)[i].m_steps_sensitive) {
          int det = st.m_det_id >> 28;
          if (st.m_g4_pid == q.m_g4_id && (det == 8 || det == 9 || det == 10) && st.m_edep.m_total > 0) e_mu[b] += st.m_edep.m_total;
        }
      }
    }
  }
  double tot = n_mu[0] + n_mu[1] + n_mu[2] + n_mu[3], etot = e_mu[0] + e_mu[1] + e_mu[2] + e_mu[3];
  printf("%ld interactions. Muons crossing the HGCal face at 1.5 < |eta| < 3.045: %.4f per interaction (x200: %.1f); "
         "by pT at the face <0.5 / 0.5-1 / 1-2 / >2 GeV: %.4f %.4f %.4f %.4f\n", n_ev, tot / n_ev, 200 * tot / n_ev,
         n_mu[0] / n_ev, n_mu[1] / n_ev, n_mu[2] / n_ev, n_mu[3] / n_ev);
  printf("their own visible HGCal energy: %.2f MeV per interaction (x200: %.2f GeV); crossings outside that eta range: %.3f per interaction\n",
         etot / n_ev, 0.2 * etot / n_ev, n_other_fwd / n_ev);
}
