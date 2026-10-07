#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
#include "TFile.h"
#include "TTree.h"
#include <algorithm>
#include <cstdio>
#include <map>
#include <vector>
#include <cmath>
// CAN ANYTHING SEE A DIVERGED CANDIDATE, WITHOUT TRUTH, OVER SEVERAL LAYERS?
//
// Three hedges against the wrong-hit failure were built and all three are null
// (reserved hole slot forward and inward, best-short). One cause: an alternative
// pays only if the diverged branch is eventually RECOGNISED, and nothing
// recognises it -- the first wrong hit sits at chi2 3.75, inside the range
// correct hits occupy.
//
// So the question is not "keep more alternatives" but "is the information there
// at all". By construction it cannot be per-hit: one residual is exactly what
// was shown to be uninformative. It has to read the SEQUENCE of a candidate's
// accepted hits. Four quantities, all computable at run time with no truth:
//
//   chi2_mean3   running mean chi2 over the last three accepted hits. The
//                reference -- if this separated, a tighter cut would have worked.
//   sign_run     fraction of consecutive accepted hits whose SIGNED q residual
//                keeps the same sign. chi2 squares the residual and throws the
//                sign away, so this is information chi2 cannot carry. A track
//                following a wrong trajectory should drift, giving correlated
//                signs; a correct one should alternate at random, 0.5.
//   bias3        |mean of the last three signed q pulls|. The same idea as a
//                magnitude rather than a count.
//   shrink       ln of how fast sigma_q_trk has fallen over the last three
//                steps. A candidate absorbing wrong hits shrinks its covariance
//                around a wrong trajectory.
//
// Truth is used ONLY to label, never to compute. A candidate is DIVERGED from
// its first accepted non-MC hit onward; hits before that, and all hits of a
// candidate that never took a wrong one, are PURE.
static double q(std::vector<double> u, double p) {
  if (u.empty()) return 0;
  std::sort(u.begin(), u.end());
  return u[(size_t)(p / 100. * (u.size() - 1))];
}
struct H { int step, layer; bool mc, acc; float chi2, dq_s, sig_q, hit_hl; };

void an_diverge(const char *fn = "/tmp/claude-1000/s8/val-diverge.root") {
  std::map<long long, std::vector<H>> C;
  TFile f(fn);
  TTree *t = (TTree *)f.Get("search");
  if (!t) { printf("no 'search' tree in %s\n", fn); return; }
  ValSearchHit *h = nullptr;
  t->SetBranchAddress("h", &h);
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (!h->had_kalman || h->chi2 < -900 || h->global_seed < 0 || h->event < 0) continue;
    C[(long long)h->event * 1000000LL + h->global_seed].push_back(
        {h->step, h->layer, h->mc_match, h->accepted, h->chi2, h->dq_s, h->sigma_q_trk, h->hit_q_half_len});
  }

  // pure[k] / div[k] for the four quantities
  std::vector<double> c2_p, c2_d, sr_p, sr_d, bi_p, bi_d, sh_p, sh_d;
  long n_cand = 0, n_div = 0, n_pure = 0;

  for (auto &kv : C) {
    auto v = kv.second;
    std::sort(v.begin(), v.end(), [](const H &a, const H &b) { return a.step < b.step; });
    // The accepted chain, in order, COLLAPSED TO ONE ENTRY PER LAYER.
    //
    // Necessary, not cosmetic. The in-layer combinatorial takes up to four hits
    // in one layer, so consecutive accepted hits are frequently the SAME
    // crossing seen twice -- their residuals share a propagation and are
    // correlated by construction. A sign run computed over them measures how
    // often the search took an overlap, not whether the track is drifting. The
    // first version of this macro did exactly that and reported pure tracks at a
    // sign run of 1.000, which is the artefact and not a signal.
    //
    // The question is about a trend over LAYERS, so one entry per layer it is.
    std::vector<H> a;
    {
      int last_layer = -1;
      for (auto &x : v) {
        if ( ! (x.acc && x.chi2 < 30.f)) continue;
        if (x.layer == last_layer) continue;   // keep the first of a layer
        last_layer = x.layer;
        a.push_back(x);
      }
    }
    if ((int)a.size() < 5) continue;   // need a window to look back over
    ++n_cand;
    int first_wrong = -1;
    for (size_t i = 0; i < a.size(); ++i) if (!a[i].mc) { first_wrong = (int)i; break; }
    if (first_wrong >= 0) ++n_div; else ++n_pure;

    for (size_t i = 3; i < a.size(); ++i) {
      // label THIS hit's history window: diverged if the wrong hit is already behind it
      const bool diverged = (first_wrong >= 0 && (int)i > first_wrong);
      // only compare like with like: a pure window from a candidate that never diverged
      if (!diverged && first_wrong >= 0) continue;

      double m = 0; for (int k = 0; k < 3; ++k) m += a[i - k].chi2; m /= 3;
      int same = 0; for (int k = 0; k < 3; ++k)
        same += (a[i - k].dq_s >= 0) == (a[i - k - 1].dq_s >= 0);
      double sr = same / 3.0;
      double b = 0; int nb = 0;
      for (int k = 0; k < 3; ++k) {
        const double sg = std::sqrt(a[i-k].sig_q * a[i-k].sig_q +
                                    a[i-k].hit_hl * a[i-k].hit_hl / 3.0);
        if (sg > 0) { b += a[i-k].dq_s / sg; ++nb; }
      }
      b = nb ? std::fabs(b / nb) : 0;
      double sh = (a[i-3].sig_q > 0 && a[i].sig_q > 0) ? std::log(a[i-3].sig_q / a[i].sig_q) : 0;

      if (diverged) { c2_d.push_back(m); sr_d.push_back(sr); bi_d.push_back(b); sh_d.push_back(sh); }
      else          { c2_p.push_back(m); sr_p.push_back(sr); bi_p.push_back(b); sh_p.push_back(sh); }
    }
  }

  printf("\n===== CAN ANYTHING SEE A DIVERGED CANDIDATE, WITHOUT TRUTH? =====\n");
  printf("%ld candidates with >=5 accepted hits: %ld pure, %ld diverged.\n",
         n_cand, n_pure, n_div);
  printf("Windows: %zu pure, %zu post-divergence. Truth labels, never computes.\n\n",
         c2_p.size(), c2_d.size());
  printf("  %-30s %9s %9s %9s | %9s %9s %9s | %7s\n",
         "quantity (last 3 acc. hits)", "pure p50", "p75", "p90", "div p50", "p75", "p90", "sep");
  auto row = [&](const char *nm, std::vector<double> &p, std::vector<double> &d) {
    // separation: shift in medians over the pure spread, a crude effect size
    const double s = q(p, 75) - q(p, 25);
    printf("  %-30s %9.3f %9.3f %9.3f | %9.3f %9.3f %9.3f | %7.2f\n", nm,
           q(p,50), q(p,75), q(p,90), q(d,50), q(d,75), q(d,90),
           s > 0 ? (q(d,50) - q(p,50)) / s : 0.0);
  };
  row("chi2, mean of last 3", c2_p, c2_d);
  row("signed-dq sign run [0..1]", sr_p, sr_d);
  row("|mean signed q pull|", bi_p, bi_d);
  row("ln(sigma_q shrink over 3)", sh_p, sh_d);
  printf("\n  sep = (diverged median - pure median) / pure IQR. A quantity that cannot\n");
  printf("  separate sits near 0. chi2 is the reference: it is KNOWN not to separate\n");
  printf("  per hit, and the question is whether a WINDOW of it does any better.\n");
  printf("  A sign run near 0.5 is a coin flip, i.e. no drift; above it means the\n");
  printf("  residuals are correlated in sign, which is what following a wrong\n");
  printf("  trajectory should look like and what chi2 discards.\n\n");
}
