// mu_origin.C: for the muons whose path crosses the HGCal face at 1.5 < |eta| < 3.045 (born at |z| < 320 cm): the
// PDG code and birth pT of their parent, and where they were born (|z|, r).
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
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
void mu_origin() {
  long n_ev = 0;
  std::map<int, long> parent; std::map<std::string, long> where, ppt; long n = 0, prim = 0;
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
        bool ok = false;
        for (size_t j = 1; j < path.size(); ++j) {
          double z0 = std::abs(path[j - 1].m_x.fZ), z1 = std::abs(path[j].m_x.fZ);
          if (z0 < 320 && z1 >= 320) {
            double w = (320 - z0) / (z1 - z0);
            double x = path[j - 1].m_x.fX + w * (path[j].m_x.fX - path[j - 1].m_x.fX), y = path[j - 1].m_x.fY + w * (path[j].m_x.fY - path[j - 1].m_x.fY);
            double eta = std::asinh(320 / std::hypot(x, y));
            ok = eta > 1.5 && eta < 3.045;
            break;
          }
        }
        if (!ok) continue;
        ++n;
        if (q.m_parent <= 0) { ++prim; continue; }
        const G4S_Particle &par = (*p)[q.m_parent];
        ++parent[par.m_pdg];
        double r = std::hypot(q.m_x_birth.fX, q.m_x_birth.fY), z = std::abs(q.m_x_birth.fZ);
        where[r < 5 ? "r<5cm" : r < 30 ? "r 5-30" : r < 120 ? "r 30-120" : "r>120"]++;
        where[z < 50 ? "|z|<50" : z < 150 ? "|z| 50-150" : z < 270 ? "|z| 150-270" : "|z| 270-320"]++;
        double pt = std::hypot(par.m_p_birth.fX, par.m_p_birth.fY);
        ppt[pt < 0.1 ? "parent pT<0.1" : pt < 0.3 ? "parent pT 0.1-0.3" : pt < 1 ? "parent pT 0.3-1" : "parent pT>1"]++;
      }
    }
  }
  printf("%lld interactions, %ld muons crossing the HGCal face in acceptance (%ld primaries)\n", (long long)n_ev, n, prim);
  printf("parent PDG:"); for (auto &x : parent) if (x.second * 50 > n) printf(" %d:%ld", x.first, x.second); printf("\n");
  for (auto &x : where) printf("  born %-12s %ld\n", x.first.c_str(), x.second);
  for (auto &x : ppt) printf("  %-18s %ld\n", x.first.c_str(), x.second);
}
