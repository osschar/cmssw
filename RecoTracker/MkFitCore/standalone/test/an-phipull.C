// The PHI PULL: is sigma_phi_trk the right SIZE?
//
// This is the phi twin of the q measurement that found the window covariance
// 1.4x too small centrally and 22x at |eta| > 1.6, and it exists to settle one
// question. Replacing the flat phi constant with the per-hit covariance extent
// costs tracks unless dphi_track is doubled, and the recovery saturates at 2x.
// Two readings fit that equally:
//
//   (a) sigma_phi_trk is ~2x too small  -> core/sigma comes out ~2
//   (b) sigma_phi_trk is right and a 3 sigma window is simply not enough
//       containment, 2x collecting the Gaussian tail -> core/sigma ~1
//
// Only the pull separates them, so do not tune the phi cut before reading it.
//
//   root -l -b -q 'an-phipull.C("phipull.root")'
//
// DISCIPLINE, inherited from the q analysis and not optional:
//  * use the SIGNED residual (dphi_s). The code abs()es dphi before the trace,
//    so a bias is invisible in that branch.
//  * quote a ROBUST width, (p84-p16)/2 -- never an RMS. 27 % of hits on a sim
//    track belong to another particle.
//  * CLAMP mc_match geometrically. It is similarity, not a binding: a hit of
//    the right label anywhere in the layer matches, including 600 cm away.
//  * report trk_frac = sigma_trk^2/sigma_tot^2. A row where the HIT term
//    dominates says nothing about the covariance, however pretty the ratio.

#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {
  float quant(std::vector<float> &v, double f) {
    if (v.empty()) return 0.f;
    size_t k = (size_t)(f * (v.size() - 1));
    std::nth_element(v.begin(), v.begin() + k, v.end());
    return v[k];
  }
  struct Grp { const char *name; int l0, l1; };
}

void an_phipull(const char *fn = "phipull.root") {
  TFile f(fn);
  TTree *t = (TTree *)f.Get("search");
  if (!t) { printf("no 'search' tree in %s\n", fn); return; }

  // 53 M entries at 1.2 GB: read ONLY the branches used, or this takes minutes.
  t->SetBranchStatus("*", 0);
  for (const char *b : {"h.layer", "h.mc_match", "h.is_barrel", "h.dphi_s", "h.dq_s",
                        "h.sigma_phi_trk", "h.hit_phi_sigma", "h.hit_q_half_len", "h.theta"})
    t->SetBranchStatus(b, 1);

  int layer; bool mc_match, is_barrel;
  float dphi_s, dq_s, sigma_phi_trk, hit_phi_sigma, hit_q_half_len, theta;
  t->SetBranchAddress("h.layer", &layer);
  t->SetBranchAddress("h.mc_match", &mc_match);
  t->SetBranchAddress("h.is_barrel", &is_barrel);
  t->SetBranchAddress("h.dphi_s", &dphi_s);
  t->SetBranchAddress("h.dq_s", &dq_s);
  t->SetBranchAddress("h.sigma_phi_trk", &sigma_phi_trk);
  t->SetBranchAddress("h.hit_phi_sigma", &hit_phi_sigma);
  t->SetBranchAddress("h.hit_q_half_len", &hit_q_half_len);
  t->SetBranchAddress("h.theta", &theta);

  // ETA BINS, matching the recorded surface-q control rows so the two can be
  // compared directly: that entry reports sigma_phi ratios 1.05 / 1.17 / 1.51 /
  // 1.83 over |eta| 0-0.8 / 0.8-1.6 / 1.6-2.0 / >2.0, i.e. a phi covariance that
  // narrows with eta, roughly as 1 + 0.17 eta^2. A REGIONAL split aliases this,
  // because the forward regions are where high eta lives -- so bin in eta first
  // and only then ask whether any region is anomalous on top of it.
  const double eta_edge[] = {0.0, 0.8, 1.6, 2.0, 99.0};
  const char  *eta_name[] = {"|eta| 0-0.8", "0.8-1.6", "1.6-2.0", "> 2.0"};
  const int NE = 4;
  std::vector<float> eres[NE], esig[NE], ehsig[NE];

  const Grp grps[] = {{"PixB 0-3", 0, 3},   {"TBPS-P 4,6,8", 4, 9},
                      {"TBPS-S 5,7,9", 4, 9}, {"TB2S 10-15", 10, 15},
                      {"fwd pix 16-27", 16, 27}, {"TEDD 28-37", 28, 37}};
  const int NG = sizeof(grps) / sizeof(grps[0]);
  std::vector<float> res[NG], sig[NG], hsig[NG];
  long n_clamped = 0, n_tot = 0;

  const Long64_t N = t->GetEntries();
  for (Long64_t i = 0; i < N; ++i) {
    t->GetEntry(i);
    if (!mc_match || sigma_phi_trk <= 0.f || dphi_s < -900.f) continue;
    ++n_tot;
    // GEOMETRIC CLAMP: a "matched" hit further in q than the hit's own extent
    // allows by a wide margin is a label coincidence, not a measurement.
    if (hit_q_half_len > 0.f && std::fabs(dq_s) > 10.f * hit_q_half_len) { ++n_clamped; continue; }
    if (theta > -900.f && theta > 0.f && theta < M_PI) {
      const double et = std::fabs(-std::log(std::tan(0.5 * (double)theta)));
      for (int e = 0; e < NE; ++e)
        if (et >= eta_edge[e] && et < eta_edge[e + 1]) {
          eres[e].push_back(dphi_s); esig[e].push_back(sigma_phi_trk);
          ehsig[e].push_back(hit_phi_sigma); break;
        }
    }
    for (int g = 0; g < NG; ++g) {
      if (layer < grps[g].l0 || layer > grps[g].l1) continue;
      if (g == 1 && (layer % 2) != 0) continue;   // TBPS-P: even
      if (g == 2 && (layer % 2) != 1) continue;   // TBPS-S: odd
      res[g].push_back(dphi_s); sig[g].push_back(sigma_phi_trk);
      hsig[g].push_back(hit_phi_sigma);
    }
  }
  printf("entries %lld, mc-matched with a track sigma %ld, clamped away %ld (%.1f %%)\n\n",
         (long long)N, n_tot, n_clamped, n_tot ? 100.0 * n_clamped / n_tot : 0.0);
  printf("%-16s %8s %11s %11s %11s %8s %9s\n",
         "region", "n", "bias[urad]", "core[urad]", "sig_trk", "RATIO", "trk_frac");
  for (int g = 0; g < NG; ++g) {
    if (res[g].size() < 50) continue;
    std::vector<float> a = res[g];
    const float p16 = quant(a, 0.16), p50 = quant(a, 0.50), p84 = quant(a, 0.84);
    const float core = 0.5f * (p84 - p16);
    std::vector<float> b = sig[g];  const float st = quant(b, 0.50);
    std::vector<float> c = hsig[g]; const float sh = quant(c, 0.50);
    const float tot2 = st * st + sh * sh;
    printf("%-16s %8zu %11.2f %11.2f %11.2e %8.2f %9.3f\n",
           grps[g].name, res[g].size(), 1e6 * p50, 1e6 * core, st,
           core / st, tot2 > 0 ? st * st / tot2 : 0.f);
  }
  printf("\n--- BINNED IN |eta|, the axis the regional split aliases ---\n");
  printf("%-16s %8s %11s %11s %11s %8s %9s\n",
         "eta bin", "n", "bias[urad]", "core[urad]", "sig_trk", "RATIO", "trk_frac");
  for (int e = 0; e < NE; ++e) {
    if (eres[e].size() < 50) continue;
    std::vector<float> a = eres[e];
    const float p16 = quant(a, 0.16), p50 = quant(a, 0.50), p84 = quant(a, 0.84);
    const float core = 0.5f * (p84 - p16);
    std::vector<float> b = esig[e];  const float st = quant(b, 0.50);
    std::vector<float> c = ehsig[e]; const float sh = quant(c, 0.50);
    const float tot2 = st * st + sh * sh;
    printf("%-16s %8zu %11.2f %11.2f %11.2e %8.2f %9.3f\n",
           eta_name[e], eres[e].size(), 1e6 * p50, 1e6 * core, st,
           core / st, tot2 > 0 ? st * st / tot2 : 0.f);
  }
  printf("  recorded for comparison (surface-q control): 1.05 / 1.17 / 1.51 / 1.83\n");

  printf("\nRATIO is the measured core over the QUOTED 1-sigma track error.\n"
         "  ~1 -> the covariance is right; the 2x that recovers the per-hit cut is\n"
         "        containment, so take the 13 %% and widen to ~6 sigma.\n"
         "  ~2 -> the covariance is SHORT, and that 2x is a symptom to fix, not a\n"
         "        setting to keep -- the phi analogue of the surface reference.\n"
         "trk_frac near 0 means the HIT term dominates: that row says nothing.\n");
}

void an_phipull() { an_phipull("phipull.root"); }
