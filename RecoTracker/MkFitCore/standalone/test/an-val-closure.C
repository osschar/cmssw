// Analysis of the propagation-closure sample. Needs NO mkFit build: only the
// dictionary library, so this runs in a bare root.exe.
//
//   . /home/matevz/root7.env ; unset DISPLAY
//   cd /foo/matevz/mic-dev/current/src/standalone
//   LD_LIBRARY_PATH=. root.exe -l -b -q \
//     '../RecoTracker/MkFitCore/standalone/test/an-val-closure.C("val-closure.root", 1)'
//
// Writes <prefix>.root (histograms) and <prefix>.txt (the numbers behind them).
//
// WHAT IS MEASURED.  Material OFF and a uniform field make the propagation a
// diffeomorphism, so the round trip A -> B -> A' must return the input exactly.
// |d pos| = |r(A') - r(A)| in cm; the EXPECTATION IS ZERO and every non-zero
// value is float arithmetic.
//
// The floor is computable rather than guessed: single precision has a 24-bit
// significand, so at a coordinate of magnitude |x| adjacent representable
// floats differ by 2^floor(log2|x|) * 2^-23 -- one ULP. A handful of operations
// lands a few ULP. So the second histogram, |d pos| / ULP, is the one with a
// clean expectation: O(1-10), dimensionless, independent of where in the
// detector the trial sat.

void an_val_closure(const char *in_file = "val-closure.root",
                    int cfg = 1,
                    const char *prefix = "110-pos-closure",
                    double pt_min = 0.3) {
  gSystem->Load("libMkFitRootDataFormats.so");

  FILE *log = fopen(Form("%s.txt", prefix), "w");
  auto P = [&](const char *fmt, ...) {
    va_list a; va_start(a, fmt); vprintf(fmt, a); va_end(a);
    va_start(a, fmt); vfprintf(log, fmt, a); va_end(a);
  };

  // ---- the configuration legend, so the plot can label itself
  float r_a = 0, r_b = 0;
  {
    ROOT::RDataFrame g("cfgs", in_file);
    auto rows = g.Filter(Form("g.cfg == %d", cfg)).Take<float>("g.r_a");
    auto rowb = g.Filter(Form("g.cfg == %d", cfg)).Take<float>("g.r_b");
    if (rows->empty()) { printf("*** no such cfg %d\n", cfg); return; }
    r_a = rows->at(0); r_b = rowb->at(0);
  }

  ROOT::RDataFrame df("closure", in_file);
  const auto n_all = *df.Count();

  auto d = df.Filter(Form("c.cfg == %d && c.pt >= %g", cfg, pt_min), "cfg and pT")
             .Filter("c.fail_oa == 0 && c.fail_ab == 0 && c.fail_ba == 0", "all three legs clean")
             // |d pos| : the closure residual, cm. Expectation 0.
             .Define("dpos",
                     [](const ROOT::RVec<float> &a, const ROOT::RVec<float> &c) {
                       const double dx = c[0] - a[0], dy = c[1] - a[1], dz = c[2] - a[2];
                       return std::sqrt(dx * dx + dy * dy + dz * dz);
                     },
                     {"c.par_a[6]", "c.par_c[6]"})
             // One ULP at the coordinate magnitude of this trial, cm.
             .Define("ulp",
                     [](const ROOT::RVec<float> &a) {
                       const double m = std::max({std::fabs(a[0]), std::fabs(a[1]), std::fabs(a[2])});
                       return m > 0 ? std::ldexp(1.0, std::ilogb(m) - 23) : 0.0;
                     },
                     {"c.par_a[6]"})
             .Define("dpos_ulp", "ulp > 0 ? dpos / ulp : 0.0")
             .Define("log10pt", "std::log10((double) c.pt)")
             .Define("abseta", "std::fabs((double) c.eta)")
             .Define("log10_ulp", "dpos_ulp > 0 ? std::log10(dpos_ulp) : -3.0")
             // Fairness check on the test itself: a track from the origin reaches
             // at most 2*R_c, so near tang_b = 1 plane B is TANGENT to the
             // trajectory and the solve is ill-posed by geometry, not by code.
             // R_c [cm] = pT / (0.01 * sol * B) = 87.515 * pT at B = 3.8112 T.
             .Define("tang_b", "2.0 * 87.515 * c.pt / c.r_b");

  const auto n_cfg = *df.Filter(Form("c.cfg == %d && c.pt >= %g", cfg, pt_min)).Count();
  const auto n_use = *d.Count();

  // Log-binned axes.
  auto logbins = [](int n, double lo, double hi, std::vector<double> &e) {
    e.resize(n + 1);
    for (int i = 0; i <= n; ++i) e[i] = std::pow(10.0, lo + (hi - lo) * i / n);
  };
  (void) logbins;
  auto dl = d.Define("log10_dpos", "dpos > 0 ? std::log10(dpos) : -9.0");

  auto h_dpos = dl.Histo1D({"dpos",
                            Form("closure residual, round trip r=%.0f #rightarrow %.0f #rightarrow %.0f cm;"
                                 "log_{10}( |#Deltapos| / 1 cm );trials", r_a, r_b, r_a),
                            140, -8.0, 4.0}, "log10_dpos");
  auto h_ulp  = dl.Histo1D({"dpos_ulp",
                            Form("closure residual in float32 ULP  (r=%.0f #rightarrow %.0f cm);"
                                 "log_{10}( |#Deltapos| / ULP );trials", r_a, r_b),
                            160, -1.0, 8.0}, "log10_ulp");


  // ---- WHERE. Declared expansion, spent because the tail above fired it.
  auto map = d.Profile2D({"map",
      Form("median-ish log_{10}( |#Deltapos| / ULP ), r=%.0f #rightarrow %.0f cm;"
           "log_{10} p_{T} [GeV];#eta;#LTlog_{10}(|#Deltapos|/ULP)#GT", r_a, r_b),
      28, -0.523, 2.301, 32, -4.0, 4.0},
      "log10pt", "c.eta", "log10_ulp");
  auto inc = d.Profile1D({"inc",
      "log_{10}( |#Deltapos| / ULP ) vs incidence at A;|cos(incidence)| at A;"
      "#LTlog_{10}(|#Deltapos|/ULP)#GT", 50, 0.0, 1.0},
      "c.cos_inc_a", "log10_ulp");

  auto q = d.Take<double>("dpos_ulp");
  auto qc = d.Take<double>("dpos");

  // ---- report
  P("### %s -- propagateHelixToPlaneMPlex round-trip closure\n", prefix);
  P("###\n");
  P("### definition : |d pos| = |r(A') - r(A)| after O->A->B->A', A and B planes at\n");
  P("###              r = %.0f and %.0f cm, normals radial. Material OFF, uniform B.\n", r_a, r_b);
  P("### unit       : cm  (and dimensionless when divided by one float32 ULP)\n");
  P("### expectation: exactly 0; the arithmetic floor is a few ULP\n");
  P("### population : pT >= %g GeV; %lld trials in this cfg of %lld in the file,\n"
    "###              %lld with all three legs clean\n",
    pt_min, (long long) n_cfg, (long long) n_all, (long long) n_use);
  P("###\n");

  std::vector<double> v = *q, w = *qc;
  std::sort(v.begin(), v.end());
  std::sort(w.begin(), w.end());
  auto pct = [](std::vector<double> &s, double p) {
    if (s.empty()) return 0.0;
    double x = 0.01 * p * (s.size() - 1);
    size_t i = (size_t) x;
    return i + 1 >= s.size() ? s.back() : s[i] * (1 - (x - i)) + s[i + 1] * (x - i);
  };

  P("%-22s %10s %10s %10s %10s %10s\n", "", "median", "p90", "p99", "p99.9", "max");
  P("%-22s %10.3g %10.3g %10.3g %10.3g %10.3g\n", "|d pos|  [cm]",
    pct(w, 50), pct(w, 90), pct(w, 99), pct(w, 99.9), w.empty() ? 0 : w.back());
  P("%-22s %10.3g %10.3g %10.3g %10.3g %10.3g\n", "|d pos| / ULP",
    pct(v, 50), pct(v, 90), pct(v, 99), pct(v, 99.9), v.empty() ? 0 : v.back());

  long n10 = 0, n100 = 0, n1e4 = 0;
  for (double x : v) { if (x > 10) ++n10; if (x > 100) ++n100; if (x > 1e4) ++n1e4; }
  P("\nfraction above  10 ULP : %7.3f %%\n", 100.0 * n10 / v.size());
  P("fraction above 100 ULP : %7.3f %%\n", 100.0 * n100 / v.size());
  P("fraction above 1e4 ULP : %7.3f %%\n", 100.0 * n1e4 / v.size());


  // ---- banded tables: which variable actually sorts the tail
  {
    auto band = [&](const char *what, const char *col, std::vector<double> edges) {
      P("\n  median and p99 of |d pos|/ULP, banded in %s\n", what);
      P("  %-18s %9s %12s %12s %12s\n", what, "n", "median", "p99", "max");
      for (size_t i = 0; i + 1 < edges.size(); ++i) {
        auto sel = d.Filter(Form("%s >= %g && %s < %g", col, edges[i], col, edges[i + 1]));
        auto tk = sel.Take<double>("dpos_ulp");
        std::vector<double> u = *tk;
        if (u.empty()) continue;
        std::sort(u.begin(), u.end());
        P("  %7.3g .. %-8.3g %9zu %12.4g %12.4g %12.4g\n",
          edges[i], edges[i + 1], u.size(), pct(u, 50), pct(u, 99), u.back());
      }
    };
    band("|eta|", "abseta", {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0});
    band("pT [GeV]", "c.pt", {0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 200.0});
    band("|cos incidence| at A", "c.cos_inc_a", {0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.01});
    band("2*R_c / r_B  (1 = tangent)", "tang_b", {1.0, 1.05, 1.2, 1.5, 2.0, 4.0, 1e9});
    P("\n  and pT bands again, but only for tang_b > 1.5 (plane B comfortably reachable):\n");
    {
      auto dd = d.Filter("tang_b > 1.5");
      std::vector<double> edges = {0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 200.0};
      P("  %-18s %9s %12s %12s %12s\n", "pT [GeV]", "n", "median", "p99", "max");
      for (size_t i = 0; i + 1 < edges.size(); ++i) {
        auto tk = dd.Filter(Form("c.pt >= %g && c.pt < %g", edges[i], edges[i+1])).Take<double>("dpos_ulp");
        std::vector<double> u = *tk;
        if (u.empty()) continue;
        std::sort(u.begin(), u.end());
        P("  %7.3g .. %-8.3g %9zu %12.4g %12.4g %12.4g\n",
          edges[i], edges[i+1], u.size(), pct(u, 50), pct(u, 99), u.back());
      }
    }
  }

  TFile f(Form("%s.root", prefix), "RECREATE");
  h_dpos->Write();
  h_ulp->Write();
  map->Write();
  inc->Write();
  f.Close();
  P("\nwrote %s.root (4 histograms) and %s.txt\n", prefix, prefix);
  fclose(log);
}

void an_val_closure_all() {
  for (int i = 0; i < 7; ++i) an_val_closure("val-closure.root", i, Form("cfg%d-closure", i));
}
