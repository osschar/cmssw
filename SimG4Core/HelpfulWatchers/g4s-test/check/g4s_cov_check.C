// g4s_cov_check.C(off, on): the inert sums of "on" (inert_aggregate_descendants on) against the same sums recomputed
// from the individual descendant steps of "off" (the same events with it off): deposit, mean, covariance and time RMS
// per (particle, volume). A charged step counts as its midpoint plus D D^T / 12, a neutral one as its end point; the
// charge from the step's PDG code (TDatabasePDG; ions charged).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <utility>
#include <vector>
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TParticlePDG.h"
#include "TTree.h"
#include "SimG4Core/HelpfulWatchers/interface/G4SnitchDataFormat.h"
static bool charged(int pdg) {
  if (std::abs(pdg) >= 1000000000) return true;
  TParticlePDG *pp = TDatabasePDG::Instance()->GetParticle(pdg);
  return pp && pp->Charge() != 0;
}
void g4s_cov_check(const char *foff, const char *fon) {
  TFile a(foff), b(fon);
  TTree *ta = (TTree *)a.Get("T"), *tb = (TTree *)b.Get("T");
  std::vector<G4S_Particle> *pa = nullptr;
  std::vector<G4S_ParticleSteps> *sa = nullptr, *sb = nullptr;
  ta->SetBranchAddress("p", &pa); ta->SetBranchAddress("s", &sa); tb->SetBranchAddress("s", &sb);
  long n_sums = 0, n_missing = 0, n_neutral = 0, n_steps = 0;
  double max_de = 0, max_dm = 0, max_dc = 0, max_dt = 0;
  static const int ci[6] = {0, 1, 2, 0, 0, 1}, cj[6] = {0, 1, 2, 1, 2, 2};
  for (long ev = 0; ev < ta->GetEntries(); ++ev) {
    ta->GetEntry(ev); tb->GetEntry(ev);
    for (size_t i = 1; i < pa->size(); ++i) {
      if (!(*pa)[i].m_was_tracked) continue;
      std::map<int, std::vector<double>> acc;  // volume -> w, s1[4], s2[6], s2t (raw, relative to the first point)
      std::map<int, std::vector<double>> ref;
      for (const G4S_Step &st : (*sa)[i].m_steps_inert) {
        if (st.m_g4_pid == (*pa)[i].m_g4_id) continue;
        const double w = st.m_edep.m_total;
        if (w <= 0) continue;
        ++n_steps;
        bool c = charged(st.m_pdg);
        if (!c) ++n_neutral;
        double x[4], dd[3] = {0, 0, 0};
        if (c) {
          x[0] = 0.5 * (st.m_x_beg.fX + st.m_x_end.fX); x[1] = 0.5 * (st.m_x_beg.fY + st.m_x_end.fY);
          x[2] = 0.5 * (st.m_x_beg.fZ + st.m_x_end.fZ); x[3] = 0.5 * (st.m_x_beg.fT + st.m_x_end.fT);
          dd[0] = st.m_x_end.fX - st.m_x_beg.fX; dd[1] = st.m_x_end.fY - st.m_x_beg.fY; dd[2] = st.m_x_end.fZ - st.m_x_beg.fZ;
        } else {
          x[0] = st.m_x_end.fX; x[1] = st.m_x_end.fY; x[2] = st.m_x_end.fZ; x[3] = st.m_x_end.fT;
        }
        auto &r = ref[st.m_volume];
        if (r.empty()) r.assign(x, x + 4);
        auto &v = acc[st.m_volume];
        if (v.empty()) v.assign(12, 0.0);
        double d[4];
        for (int k = 0; k < 4; ++k) { d[k] = x[k] - r[k]; v[1 + k] += w * d[k]; }
        for (int cc = 0; cc < 6; ++cc) v[5 + cc] += w * (d[ci[cc]] * d[cj[cc]] + dd[ci[cc]] * dd[cj[cc]] / 12.0);
        v[11] += w * d[3] * d[3];
        v[0] += w;
      }
      for (const G4S_InertSum &is : (*sb)[i].m_inert_sums) {
        ++n_sums;
        auto it = acc.find(is.m_volume);
        if (it == acc.end()) { ++n_missing; continue; }
        const auto &v = it->second; const auto &r = ref[is.m_volume];
        const double W = v[0];
        double mu[4];
        for (int k = 0; k < 4; ++k) mu[k] = v[1 + k] / W;
        max_de = std::max(max_de, std::abs(W - is.m_edep.m_total) / W);
        max_dm = std::max({max_dm, std::abs(r[0] + mu[0] - is.m_x_mean.fX), std::abs(r[1] + mu[1] - is.m_x_mean.fY),
                           std::abs(r[2] + mu[2] - is.m_x_mean.fZ)});
        double tr = 0, c[6];
        for (int cc = 0; cc < 6; ++cc) { c[cc] = v[5 + cc] / W - mu[ci[cc]] * mu[cj[cc]]; if (cc < 3) { c[cc] = std::max(0.0, c[cc]); tr += c[cc]; } }
        for (int cc = 0; cc < 6; ++cc) max_dc = std::max(max_dc, std::abs(c[cc] - is.m_cov[cc]) / (tr + 1e-12));
        const double trms = std::sqrt(std::max(0.0, v[11] / W - mu[3] * mu[3]));
        max_dt = std::max(max_dt, std::abs(trms - is.m_t_rms) / (trms + 1e-9));
      }
    }
  }
  printf("%lld events: %ld sums, %ld without steps in the other file; %ld descendant steps, %ld of neutral tracks\n",
         ta->GetEntries(), n_sums, n_missing, n_steps, n_neutral);
  printf("largest differences: deposit %.2g (relative), mean %.2g cm, covariance %.2g (relative to its trace), time RMS %.2g (relative)\n",
         max_de, max_dm, max_dc, max_dt);
}
