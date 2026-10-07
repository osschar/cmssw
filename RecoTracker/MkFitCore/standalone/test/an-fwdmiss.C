#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// Where the forward search loses the sim track's own hit. One row per
// LAYER-SEARCH, so it answers "where were hits missed", not "why was this track
// not found" -- the two differ, since a track survives losing a hit.
//
// Verdicts are mutually exclusive and ordered most-upstream first, as
// val_search_event() assigns them.
static double qq(std::vector<double> u, double p) {
  if (u.empty()) return 0;
  std::sort(u.begin(), u.end());
  return u[(size_t)(p / 100. * (u.size() - 1))];
}
void an_fwdmiss(const char *fn = "val-fwd-miss.root") {
  TFile f(fn);
  TTree *t = (TTree*) f.Get("miss");
  if (!t) { printf("no 'miss' tree in %s\n", fn); return; }
  ValSearchMiss *m = nullptr; t->SetBranchAddress("m", &m);
  const char *vn[8] = {"0 track never hits layer", "1 outside the window",
                       "2 in window, not scanned", "3 scanned, failed presel",
                       "4 preselected, evicted",   "5 Kalman, chi2 >= 30",
                       "6 Kalman, outranked",      "7 WON its layer"};
  // regions by the searched layer: pixel barrel, OT barrel, forward
  const char *rn[3] = {"pixel barrel 0-3", "OT barrel 4-15", "forward >= 16"};
  long V[3][8] = {}, tot[3] = {}, nolbl = 0, nohere[3] = {};
  // Split on the WSR verdict the search itself acted on. With Config::V2p2::Policy::use_wsr
  // on -- the default -- a WSR_Outside candidate scans NO hits, so its row is a
  // DECLINED search, not a failed one. Counting it makes verdict 0 look like a
  // plan defect and verdict 1 look like a window too small; it is neither.
  long W[3][3] = {}, wund[3] = {}, wout_hashit[3] = {}, ngap[3] = {};
  std::vector<double> D[3], PH[3], QQ[3];
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (m->sim_label < 0) { ++nolbl; continue; }
    if (m->layer < 0 || m->verdict < 0 || m->verdict > 7) continue;
    const int r = m->layer < 4 ? 0 : (m->layer < 16 ? 1 : 2);
    if (m->wsr >= 0 && m->wsr <= 2) ++W[r][m->wsr]; else ++wund[r];
    if (m->wsr_in_gap) ++ngap[r];
    if (m->wsr == 2) {                       // declined: it scanned nothing
      if (m->verdict != 0) ++wout_hashit[r]; // ... and the sim track HAD a hit
      continue;
    }
    ++V[r][m->verdict];
    if (m->verdict == 0) { ++nohere[r]; if (m->near_d3d > -900) D[r].push_back(m->near_d3d); }
    else {
      ++tot[r];                       // the sim track HAS a hit in this layer
      if (m->verdict == 1) { PH[r].push_back(m->sim_dphi_norm);
                             QQ[r].push_back(std::fabs(m->sim_dq_norm)); }
    }
  }
  printf("\nWHERE THE FORWARD SEARCH LOSES THE SIM TRACK'S OWN HIT\n");
  printf("%lld layer-searches, %ld with no sim label on the seed.\n", t->GetEntries(), nolbl);
  printf("\nWSR FIRST -- a WSR_Outside candidate is DECLINED before any hit is scanned,\n");
  printf("so those rows are excluded from everything below.\n");
  printf("  %-26s", "WSR of searches w/ a label");
  for (int r = 0; r < 3; ++r) printf(" | %18s", rn[r]);
  printf("\n");
  const char *wnm[3] = {"inside", "edge", "OUTSIDE -- declined"};
  for (int w = 0; w < 3; ++w) {
    printf("  %-26s", wnm[w]);
    for (int r = 0; r < 3; ++r) {
      const long n = W[r][0] + W[r][1] + W[r][2];
      printf(" | %9ld %7.2f%%", W[r][w], n ? 100.*W[r][w]/n : 0.);
    }
    printf("\n");
  }
  printf("  %-26s", "of those, sim HAD a hit");
  for (int r = 0; r < 3; ++r)
    printf(" | %9ld %7.2f%%", wout_hashit[r], W[r][2] ? 100.*wout_hashit[r]/W[r][2] : 0.);
  printf("   <- WSR false miss\n");
  printf("  %-26s", "in_gap (disc r-hole)");
  for (int r = 0; r < 3; ++r) printf(" | %9ld %8s", ngap[r], "");
  printf("\n  %-26s", "no WSR recorded");
  for (int r = 0; r < 3; ++r) printf(" | %9ld %8s", wund[r], "");
  printf("\n");
  printf("Percentages are over searches where the sim track ACTUALLY HAS a hit in the\n");
  printf("searched layer -- verdict 0 is reported separately because the layer plans are\n");
  printf("deliberately inclusive and a search on a layer the track never crosses is\n");
  printf("expected, not a defect -- and the WSR has already declined most of them.\n\n");
  printf("  %-26s", "verdict");
  for (int r = 0; r < 3; ++r) printf(" | %18s", rn[r]);
  printf("\n");
  for (int v = 1; v < 8; ++v) {
    printf("  %-26s", vn[v]);
    for (int r = 0; r < 3; ++r)
      printf(" | %9ld %7.2f%%", V[r][v], tot[r] ? 100.*V[r][v]/tot[r] : 0.);
    printf("\n");
  }
  printf("  %-26s", "TOTAL with a hit here");
  for (int r = 0; r < 3; ++r) printf(" | %9ld %8s", tot[r], "");
  printf("\n  %-26s", vn[0]);
  for (int r = 0; r < 3; ++r)
    printf(" | %9ld %7.2f%%", nohere[r],
           (tot[r]+nohere[r]) ? 100.*nohere[r]/(tot[r]+nohere[r]) : 0.);
  printf("   <- of ALL searches\n");
  for (int r = 0; r < 3; ++r) {
    if (D[r].empty() && PH[r].empty()) continue;
    printf("\n  %s:\n", rn[r]);
    if (!D[r].empty())
      printf("    verdict 0 -- nearest sim hit %.2f cm away (med), p90 %.2f\n",
             qq(D[r],50), qq(D[r],90));
    if (!PH[r].empty())
      printf("    verdict 1 -- |dphi|/phi_delta med %.2f p90 %.2f ; |dq| in window\n"
             "      half-widths med %.2f p90 %.2f (n=%zu)\n",
             qq(PH[r],50), qq(PH[r],90), qq(QQ[r],50), qq(QQ[r],90), PH[r].size());
  }
  printf("\n");
}
