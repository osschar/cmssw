#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// Does ANY candidate follow a low-pT forward particle into the outer tracker?
// Input is val_search_write()'s "miss" tree (val_sr_reset / val_sr_ev / val_sr_write).
//
// an-fwdfake.C counts layer-searches over all candidates. This macro groups them by
// PARTICLE (event, sim label) and, per OT layer (TBPS 4-15 barrel strips, TEDD
// 28-37 / 50-59) in which the particle has a sim hit, takes the BEST verdict over all
// of the particle's candidates that searched that layer:
//   7 some surviving candidate holds the true hit
//   6 the true hit passed the cut, but no surviving candidate holds it
//   5 / 4 / 3 / 2 / 1 lost upstream (cut, eviction, pre-selection, not scanned, window)
//   - no candidate searched the layer at all (the beam stopped before it, or the
//     WSR declined it)
// The OT layers with a sim hit are taken in the order the particle's candidates first
// searched them (smallest step), so "1st" is the first OT layer the beam reached.
//   root -l -b -q 'an-fwdfollow.C("fwdfake-h2-20ev.root")'
void an_fwdfollow(const char *fn, float eta_lo = 1.7f, float eta_hi = 2.7f, float pt_lo = 0.0f,
                  float pt_hi = 0.9f) {
  TFile f(fn);
  TTree *t = (TTree *)f.Get("miss");
  if (!t) { printf("no 'miss' tree in %s\n", fn); return; }
  ValSearchMiss *m = nullptr;
  t->SetBranchAddress("m", &m);
  auto is_ot = [](int l) { return (l >= 4 && l < 16) || (l >= 28 && l < 38) || l >= 50; };

  struct LayerBest { int verdict = -1; int first_step = 1 << 30; bool tedd = false; };
  // (event, label) -> layer -> best
  std::map<std::pair<int, int>, std::map<int, LayerBest>> part;
  std::map<std::pair<int, int>, int> part_n_ot_sim;  // OT layers with a sim hit, as seen by any search
  const Long64_t N = t->GetEntries();
  for (Long64_t i = 0; i < N; ++i) {
    t->GetEntry(i);
    if (m->sim_label < 0 || !is_ot(m->layer)) continue;
    if (std::abs(m->eta) < eta_lo || std::abs(m->eta) >= eta_hi) continue;
    if (m->pt < pt_lo || m->pt >= pt_hi) continue;
    if (m->wsr == 2) continue;  // declined search, scanned nothing
    if (!m->has_sim_here) continue;
    auto &lb = part[{m->event, m->sim_label}][m->layer];
    lb.verdict = std::max(lb.verdict, m->verdict);
    lb.first_step = std::min(lb.first_step, m->step);
    lb.tedd = m->layer >= 28;
  }
  // Order each particle's OT layers by first search step.
  const int K = 4;
  long n_k[K] = {}, v_k[K][8] = {};
  long n_part = 0, n_part_any7 = 0, n_part_first7 = 0;
  long n_layers = 0, v_all[8] = {};
  for (auto &[key, layers] : part) {
    std::vector<std::pair<int, int>> ord;  // (first_step, layer)
    for (auto &[l, lb] : layers) ord.push_back({lb.first_step, l});
    std::sort(ord.begin(), ord.end());
    ++n_part;
    bool any7 = false;
    for (size_t k = 0; k < ord.size(); ++k) {
      int v = layers[ord[k].second].verdict;
      if (v < 0 || v > 7) continue;
      ++n_layers; ++v_all[v];
      if (v == 7) any7 = true;
      if (k < (size_t)K) { ++n_k[k]; ++v_k[k][v]; }
    }
    if (any7) ++n_part_any7;
    if (!ord.empty() && layers[ord[0].second].verdict == 7) ++n_part_first7;
  }
  printf("\nPER PARTICLE, OT layers with a sim hit that at least one candidate searched\n");
  printf("candidate %.1f <= |eta| < %.1f, %.2f <= pT < %.2f. File %s\n", eta_lo, eta_hi, pt_lo, pt_hi, fn);
  printf("  particles %ld; some candidate holds a true OT hit in >= 1 layer: %.1f %%; in the first: %.1f %%\n",
         n_part, 100. * n_part_any7 / std::max(1L, n_part), 100. * n_part_first7 / std::max(1L, n_part));
  const char *vn[8] = {"0 no sim hit", "1 outside the window", "2 in window, not scanned", "3 failed presel",
                       "4 evicted", "5 failed the cut", "6 passed, no survivor holds it", "7 a survivor holds it"};
  printf("  %-34s %9s %9s %9s %9s %9s\n", "best verdict over the candidates", "all", "1st OT", "2nd OT", "3rd OT", "4th OT");
  printf("  %-34s %9ld %9ld %9ld %9ld %9ld\n", "layers", n_layers, n_k[0], n_k[1], n_k[2], n_k[3]);
  for (int v = 1; v < 8; ++v) {
    printf("  %-34s %8.1f%%", vn[v], 100. * v_all[v] / std::max(1L, n_layers));
    for (int k = 0; k < K; ++k) printf(" %8.1f%%", 100. * v_k[k][v] / std::max(1L, n_k[k]));
    printf("\n");
  }
}
