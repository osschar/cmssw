// g4s_sums.C(file): the inert sums by material, the 8 volumes with the most energy, and what the steps branch holds.
#include <algorithm>
#include <cstdio>
#include <map>
#include <string>
#include <vector>
#include "SimG4Core/HelpfulWatchers/interface/G4SnitchReader.h"
void g4s_sums(const char *fn) {
  G4SnitchFile f(fn);  // TTree or RNTuple file
  const auto vp = f.column<std::string>("V", "path"), vm = f.column<std::string>("V", "material");
  const auto vx0 = f.column<double>("V", "x0"), vz = f.column<double>("V", "z");
  std::vector<G4S_ParticleSteps> *&s = f.s;
  std::map<std::string, double> e_mat; std::map<std::string, long> n_mat;
  std::vector<double> e_vol(vp.size(), 0);
  long n_sums = 0, n_own_inert = 0, n_sens = 0, n_path = 0; double e_all = 0, rms_sum = 0;
  for (long ev = 0; ev < f.n_events(); ++ev) {
    f.load(ev);
    for (const G4S_ParticleSteps &ps : *s) {
      n_own_inert += ps.m_steps_inert.size(); n_sens += ps.m_steps_sensitive.size(); n_path += ps.m_path.size();
      for (const G4S_InertSum &is : ps.m_inert_sums) {
        ++n_sums; e_all += is.m_edep.m_total;
        e_mat[vm[is.m_volume]] += is.m_edep.m_total; ++n_mat[vm[is.m_volume]];
        e_vol[is.m_volume] += is.m_edep.m_total;
        rms_sum += is.m_edep.m_total * std::sqrt(is.m_cov[0] + is.m_cov[1] + is.m_cov[2]);
      }
    }
  }
  long ne = f.n_events();
  printf("%s: %ld events, %zu volumes; per event: %ld sums, %ld own inert steps, %ld sensitive steps, %ld path points\n",
         fn, ne, vp.size(), n_sums / ne, n_own_inert / ne, n_sens / ne, n_path / ne);
  printf("energy-weighted spatial RMS of a sum: %.2f cm\n", rms_sum / e_all);
  std::vector<std::pair<double, std::string>> m;
  for (auto &x : e_mat) m.emplace_back(x.second, x.first);
  std::sort(m.rbegin(), m.rend());
  printf("summed energy by material (share, sums):\n");
  for (size_t i = 0; i < m.size() && i < 8; ++i) printf("  %5.1f %%  %-28s %8ld\n", 100 * m[i].first / e_all, m[i].second.c_str(), n_mat[m[i].second]);
  std::vector<int> idx(vp.size()); for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&](int a, int b) { return e_vol[a] > e_vol[b]; });
  printf("volumes with the most summed energy:\n");
  for (int k = 0; k < 6; ++k) { int i = idx[k]; std::string p = vp[i]; if (p.size() > 150) p = "..." + p.substr(p.size() - 147);
    printf("  %5.2f %%  x0 %6.2f cm  z %7.1f cm  %s\n", 100 * e_vol[i] / e_all, vx0[i], vz[i], p.c_str()); }
}
