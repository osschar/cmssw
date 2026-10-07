#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// Forward search, per LAYER-SEARCH, for the population CMSSW MTV calls fakes
// at: candidate pT < pt_cut and eta_lo <= |eta| < eta_hi, against the same
// layers at pT >= pt_cut. Input is val_search_write()'s "miss" tree.
//
// Two questions, per layer group:
//  - where the sim track has NO hit in the layer (verdict 0), how often a wrong
//    hit passes the acceptance cut and is kept, and at what chi2 -- in the
//    linear layer-step score a hole beats a hit above chi2 = (hit_bonus +
//    miss_penalty) / chi2_weight, 11 at the defaults;
//  - where it HAS one, where that hit is lost (verdicts 1-7), and how often a
//    wrong hit is kept in its place.
// All candidates are counted, not only those that end on a final track.
//   root -l -b -q 'an-fwdfake.C("fwdfake-b3-20ev.root")'
static double afq(std::vector<double> u, double p) {
  if (u.empty()) return -1;
  std::sort(u.begin(), u.end());
  return u[(size_t)(p / 100. * (u.size() - 1))];
}

void an_fwdfake(const char *fn, float eta_lo = 1.7f, float eta_hi = 2.7f, float pt_cut = 0.9f,
                float hole_chi2 = 11.f) {
  TFile f(fn);
  TTree *t = (TTree *)f.Get("miss");
  if (!t) { printf("no 'miss' tree in %s\n", fn); return; }
  ValSearchMiss *m = nullptr;
  t->SetBranchAddress("m", &m);
  // layer groups
  const int NG = 5;
  const char *gn[NG] = {"pixel barrel 0-3", "TBPS 4-9", "TB2S 10-15", "FPix 16-27,38-49", "TEDD 28-37,50-59"};
  auto grp = [](int l) {
    if (l < 4) return 0;
    if (l < 10) return 1;
    if (l < 16) return 2;
    if ((l >= 28 && l < 38) || l >= 50) return 4;
    return 3;
  };
  // [population 0 = low pT forward, 1 = high pT same eta][group]
  long n[2][NG] = {}, v[2][NG][8] = {};
  long h0_pass[2][NG] = {}, h0_kept[2][NG] = {}, h0_kept_below[2][NG] = {};
  long h1_kept_wrong[2][NG] = {}, h1_kept_wrong_only[2][NG] = {};
  std::vector<double> c2w[2][NG], c2mc[2][NG], dphin[2][NG], dqn[2][NG];
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (m->sim_label < 0 || m->layer < 0 || m->verdict < 0 || m->verdict > 7) continue;
    if (m->wsr == 2) continue;  // declined, scanned nothing
    const float ae = std::fabs(m->eta);
    if (ae < eta_lo || ae >= eta_hi) continue;
    const int p = m->pt < pt_cut ? 0 : 1, g = grp(m->layer);
    ++n[p][g];
    ++v[p][g][m->verdict];
    if (m->verdict == 0) {
      if (m->n_pass_wrong > 0) ++h0_pass[p][g];
      if (m->n_kept_wrong > 0) {
        ++h0_kept[p][g];
        if (m->best_wrong_chi2 < hole_chi2) ++h0_kept_below[p][g];
      }
      if (m->n_pass_wrong > 0) c2w[p][g].push_back(m->best_wrong_chi2);
    } else {
      if (m->n_kept_wrong > 0) {
        ++h1_kept_wrong[p][g];
        if (m->verdict != 7) ++h1_kept_wrong_only[p][g];
      }
      if (m->verdict >= 5 && m->mc_chi2 > -900) c2mc[p][g].push_back(m->mc_chi2);
      if (m->verdict == 1) {
        dphin[p][g].push_back(m->sim_dphi_norm);
        dqn[p][g].push_back(std::fabs(m->sim_dq_norm));
      }
    }
  }
  const char *pn[2] = {"LOW pT", "HIGH pT"};
  printf("\nFORWARD SEARCH, layer-searches with a sim label and WSR not outside,\n");
  printf("candidate %.1f <= |eta| < %.1f; LOW pT < %.1f, HIGH pT >= %.1f. File %s\n", eta_lo, eta_hi, pt_cut, pt_cut, fn);
  for (int g = 0; g < NG; ++g) {
    if (!n[0][g] && !n[1][g]) continue;
    printf("\n  %s\n", gn[g]);
    printf("    %-44s %14s %14s\n", "", pn[0], pn[1]);
    printf("    %-44s %14ld %14ld\n", "searches", n[0][g], n[1][g]);
    auto row = [&](const char *lab, long a0, long d0, long a1, long d1) {
      printf("    %-44s %13.1f%% %13.1f%%\n", lab, d0 ? 100. * a0 / d0 : 0., d1 ? 100. * a1 / d1 : 0.);
    };
    row("sim track has NO hit in the layer (v0)", v[0][g][0], n[0][g], v[1][g][0], n[1][g]);
    row("  of those: a wrong hit passes the cut", h0_pass[0][g], v[0][g][0], h0_pass[1][g], v[1][g][0]);
    row("  of those: a wrong hit is KEPT", h0_kept[0][g], v[0][g][0], h0_kept[1][g], v[1][g][0]);
    char lab[64];
    snprintf(lab, sizeof lab, "  of those: kept, chi2 < %.0f (beats a hole)", hole_chi2);
    row(lab, h0_kept_below[0][g], v[0][g][0], h0_kept_below[1][g], v[1][g][0]);
    printf("    %-44s %14.2f %14.2f\n", "  best passing wrong chi2, median", afq(c2w[0][g], 50), afq(c2w[1][g], 50));
    const long s0 = n[0][g] - v[0][g][0], s1 = n[1][g] - v[1][g][0];
    printf("    %-44s %14ld %14ld\n", "sim track HAS a hit in the layer", s0, s1);
    const char *vl[8] = {"", "  1 outside the window", "  2 in window, not scanned", "  3 scanned, failed presel",
                         "  4 preselected, evicted", "  5 Kalman, failed the cut", "  6 passed, not kept",
                         "  7 kept"};
    for (int k = 1; k < 8; ++k) row(vl[k], v[0][g][k], s0, v[1][g][k], s1);
    row("  a wrong hit kept", h1_kept_wrong[0][g], s0, h1_kept_wrong[1][g], s1);
    row("  a wrong hit kept, the true one not", h1_kept_wrong_only[0][g], s0, h1_kept_wrong_only[1][g], s1);
    printf("    %-44s %14.2f %14.2f\n", "  true hit chi2 (v5-7), median", afq(c2mc[0][g], 50), afq(c2mc[1][g], 50));
    printf("    %-44s %14.2f %14.2f\n", "  v1: |dphi| / window half-width, median", afq(dphin[0][g], 50), afq(dphin[1][g], 50));
    printf("    %-44s %14.2f %14.2f\n", "  v1: |dq| / window half-width, median", afq(dqn[0][g], 50), afq(dqn[1][g], 50));
  }
}
