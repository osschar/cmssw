#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"

// A/B of MkBins::surface_reference_dq(): does referencing the pre-selection q
// window to the LAYER SURFACE (instead of to the fixed path length the
// covariance is transported to) flatten its eta dependence?
//
// Prediction, if the mechanism is right:
//   1. ratio_trk = (measured core width) / (quoted sigma_q) flattens in |eta|;
//   2. the flat level it lands on is the SEPARATE eta-independent factor (2-3.3),
//      which this change cannot touch;
//   3. the effect is EVEN in eta -- the correction lies in the (r,z) plane;
//   4. sigma_phi is unchanged -- nothing here touches the transverse direction;
//   5. pre-selection efficiency for TRUE hits rises, most at high |eta|.
// Any of 1-4 failing falsifies it.

static double qtl(std::vector<double> u, double p) {
  if (u.empty()) return 0;
  std::sort(u.begin(), u.end());
  return u[(size_t)(p / 100. * (u.size() - 1))];
}
// IQR-based sigma: insensitive to the 27% delta-ray contamination in the tails.
static double core(std::vector<double> &u) {
  return u.size() < 20 ? -1 : (qtl(u, 75) - qtl(u, 25)) / 1.349;
}

struct Bin {
  std::vector<double> dq_s, dphi_s, sig_q, sig_phi;
  long n_true = 0, n_true_pass = 0;
};

static void fill(const char *fn, Bin *b, const double *eb, int NE, bool split_sign,
                 long &n_all, long &n_pixb) {
  TFile f(fn);
  TTree *t = (TTree *)f.Get("search");
  ValSearchHit *h = nullptr;
  t->SetBranchAddress("h", &h);
  n_all = t->GetEntries();
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (h->layer < 0 || h->layer >= 4) continue;   // pixel barrel only
    ++n_pixb;
    if (!h->mc_match) continue;
    if (h->sigma_q_trk <= 0 || h->dq_s < -900 || h->eta < -900) continue;
    double e = split_sign ? h->eta : std::fabs(h->eta);
    int k = -1;
    for (int j = 0; j < NE; ++j) if (e >= eb[j] && e < eb[j + 1]) { k = j; break; }
    if (k < 0) continue;
    b[k].dq_s.push_back(h->dq_s);
    b[k].sig_q.push_back(h->sigma_q_trk);
    if (h->dphi_s > -900 && h->sigma_phi_trk > 0) {
      b[k].dphi_s.push_back(h->dphi_s);
      b[k].sig_phi.push_back(h->sigma_phi_trk);
    }
    ++b[k].n_true;
    if (h->passed_preselect) ++b[k].n_true_pass;
  }
}

void an_surfq(const char *f0 = "val-surfq-0.root", const char *f1 = "val-surfq-1.root") {
  const double eb[] = {0, 0.8, 1.6, 2.0, 9};
  const int NE = 4;
  const char *en[NE] = {"< 0.8", "0.8-1.6", "1.6-2.0", "> 2.0"};
  Bin B[2][NE];
  long n_all[2] = {0, 0}, n_pixb[2] = {0, 0};
  fill(f0, B[0], eb, NE, false, n_all[0], n_pixb[0]);
  fill(f1, B[1], eb, NE, false, n_all[1], n_pixb[1]);

  printf("\n================ SURFACE-REFERENCED q WINDOW: A/B ================\n");
  printf("Pixel barrel (layers 0-3), MC-matched hits. Core width = IQR/1.349 of the\n");
  printf("SIGNED residual; quoted sigma = median sigma_q_trk. ratio = core/quoted,\n");
  printf("so 1.0 is a correctly sized window and >1 is too tight.\n\n");
  printf("  %-9s %7s | %9s %9s %7s | %9s %9s %7s | %7s\n", "|eta|", "n",
         "core OFF", "quoted", "ratio", "core ON", "quoted", "ratio", "gain");
  for (int k = 0; k < NE; ++k) {
    if (B[0][k].dq_s.size() < 20 || B[1][k].dq_s.size() < 20) continue;
    double c0 = core(B[0][k].dq_s) * 1e4, q0 = qtl(B[0][k].sig_q, 50) * 1e4;
    double c1 = core(B[1][k].dq_s) * 1e4, q1 = qtl(B[1][k].sig_q, 50) * 1e4;
    printf("  %-9s %7zu | %9.1f %9.1f %7.2f | %9.1f %9.1f %7.2f | %7.2f\n",
           en[k], B[0][k].dq_s.size(), c0, q0, c0 / q0, c1, q1, c1 / q1, q1 / q0);
  }
  printf("  (core and quoted in um; 'gain' = how much the window grew)\n");

  printf("\n  CONTROL -- sigma_phi must NOT move (the correction is in the r-z plane):\n");
  printf("  %-9s | %9s %9s %7s | %9s %9s %7s\n", "|eta|",
         "core OFF", "quoted", "ratio", "core ON", "quoted", "ratio");
  for (int k = 0; k < NE; ++k) {
    if (B[0][k].dphi_s.size() < 20 || B[1][k].dphi_s.size() < 20) continue;
    double c0 = core(B[0][k].dphi_s) * 1e3, q0 = qtl(B[0][k].sig_phi, 50) * 1e3;
    double c1 = core(B[1][k].dphi_s) * 1e3, q1 = qtl(B[1][k].sig_phi, 50) * 1e3;
    printf("  %-9s | %9.3f %9.3f %7.2f | %9.3f %9.3f %7.2f\n",
           en[k], c0, q0, c0 / q0, c1, q1, c1 / q1);
  }
  printf("  (mrad)\n");

  printf("\n  PAYOFF -- pre-selection efficiency for TRUE hits:\n");
  printf("  %-9s | %8s %8s | %8s %8s | %7s\n", "|eta|", "n OFF", "pass OFF", "n ON", "pass ON", "delta");
  for (int k = 0; k < NE; ++k) {
    if (B[0][k].n_true < 20) continue;
    double p0 = 100. * B[0][k].n_true_pass / B[0][k].n_true;
    double p1 = 100. * B[1][k].n_true_pass / B[1][k].n_true;
    printf("  %-9s | %8ld %7.1f%% | %8ld %7.1f%% | %+6.1f\n",
           en[k], B[0][k].n_true, p0, B[1][k].n_true, p1, p1 - p0);
  }

  printf("\n  COST: hits scanned, all layers  %ld -> %ld  (%+.1f%%)\n",
         n_all[0], n_all[1], 100. * (n_all[1] - n_all[0]) / n_all[0]);
  printf("        hits scanned, pixel barrel %ld -> %ld  (%+.1f%%)\n",
         n_pixb[0], n_pixb[1], 100. * (n_pixb[1] - n_pixb[0]) / n_pixb[0]);

  // --- falsifier: the correction is geometric in (r,z), so it must be EVEN in eta.
  const double sb[] = {-9, -2.0, -1.6, -0.8, 0, 0.8, 1.6, 2.0, 9};
  const int NS = 8;
  Bin S[2][NS];
  long a, b;
  fill(f0, S[0], sb, NS, true, a, b);
  fill(f1, S[1], sb, NS, true, a, b);
  printf("\n  FALSIFIER -- signed eta. A geometric (r,z) correction is EVEN in eta,\n");
  printf("  so the OFF and ON ratios must each mirror about 0.\n");
  printf("  %-12s %7s | %7s %7s\n", "eta", "n", "OFF", "ON");
  for (int k = 0; k < NS; ++k) {
    if (S[0][k].dq_s.size() < 20) continue;
    double c0 = core(S[0][k].dq_s), q0 = qtl(S[0][k].sig_q, 50);
    double c1 = core(S[1][k].dq_s), q1 = qtl(S[1][k].sig_q, 50);
    printf("  %+5.1f..%+5.1f %7zu | %7.2f %7.2f\n",
           sb[k], sb[k + 1], S[0][k].dq_s.size(), c0 / q0, c1 / q1);
  }
  printf("\n");
}
