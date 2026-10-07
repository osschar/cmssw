// Does 2*R_c/r_B control the closure residual, or does pT?
//
// In ONE radius pair the two are the same variable (R_c is proportional to pT
// and r_B is fixed), so they cannot be separated. Across the seven pairs they
// can: a cell at pT 1 with r_B = 110 has the same 2*R_c/r_B as a cell at
// pT 0.5 with r_B = 50. If the tangency ratio is the controlling variable, the
// pooled table binned in it must be tight while the table binned in pT is not.
//
//   LD_LIBRARY_PATH=. root.exe -l -b -q \
//     '../RecoTracker/MkFitCore/standalone/test/an-val-closure-scan.C("val-closure.root")'

void an_val_closure_scan(const char *in_file = "val-closure.root",
                         const char *prefix = "120-closure-scan") {
  gSystem->Load("libMkFitRootDataFormats.so");
  FILE *log = fopen(Form("%s.txt", prefix), "w");
  auto P = [&](const char *fmt, ...) {
    va_list a; va_start(a, fmt); vprintf(fmt, a); va_end(a);
    va_start(a, fmt); vfprintf(log, fmt, a); va_end(a);
  };

  ROOT::RDataFrame df("closure", in_file);
  auto d = df.Filter("c.fail_oa == 0 && c.fail_ab == 0 && c.fail_ba == 0")
             .Define("dpos", [](const ROOT::RVec<float> &a, const ROOT::RVec<float> &c) {
                       const double dx = c[0]-a[0], dy = c[1]-a[1], dz = c[2]-a[2];
                       return std::sqrt(dx*dx + dy*dy + dz*dz); }, {"c.par_a[6]", "c.par_c[6]"})
             .Define("ulp", [](const ROOT::RVec<float> &a) {
                       const double m = std::max({std::fabs(a[0]), std::fabs(a[1]), std::fabs(a[2])});
                       return m > 0 ? std::ldexp(1.0, std::ilogb(m) - 23) : 0.0; }, {"c.par_a[6]"})
             .Define("dpos_ulp", "ulp > 0 ? dpos / ulp : 0.0")
             // R_c [cm] = pT / (0.01 * sol * B) = 87.515 * pT at B = 3.8112 T.
             .Define("tang_b", "2.0 * 87.515 * c.pt / c.r_b")
             .Define("log10pt_", "std::log10((double) c.pt)");

  auto med = [&](ROOT::RDF::RNode n) {
    auto tk = n.Take<double>("dpos_ulp");
    std::vector<double> u = *tk;
    if (u.empty()) return std::make_pair((size_t)0, 0.0);
    std::sort(u.begin(), u.end());
    return std::make_pair(u.size(), u[u.size() / 2]);
  };

  P("### %s -- is the closure residual controlled by pT or by tangency?\n", prefix);
  P("###\n");
  P("### quantity   : median |d pos| / (one float32 ULP), dimensionless\n");
  P("### expectation: a few; the propagation is exact with material off and uniform B\n");
  P("### 2*R_c/r_B  : how far plane B is from the track's turning circle. 1 = tangent.\n");
  P("###\n");

  const float rb[7]   = {40, 50, 70, 100, 80, 90, 110};
  const float ra[7]   = {30, 30, 30,  30, 70, 70,  70};
  std::vector<double> ptE = {0.3, 0.6, 1.0, 2.0, 4.0, 10.0, 40.0, 200.0};

  P("median |d pos|/ULP, rows = pT [GeV], cols = plane pair r_A -> r_B [cm]\n");
  P("  %-14s", "pT");
  for (int i = 0; i < 7; ++i) P("%12s", Form("%.0f->%.0f", ra[i], rb[i]));
  P("\n");
  for (size_t j = 0; j + 1 < ptE.size(); ++j) {
    P("  %5.3g .. %-6.3g", ptE[j], ptE[j+1]);
    for (int i = 0; i < 7; ++i) {
      auto r = med(d.Filter(Form("c.cfg == %d && c.pt >= %g && c.pt < %g", i, ptE[j], ptE[j+1])));
      if (r.first < 20) P("%12s", "-"); else P("%12.4g", r.second);
    }
    P("\n");
  }

  P("\npooled over all seven pairs, binned in 2*R_c/r_B:\n");
  P("  %-18s %9s %12s %12s\n", "2*R_c / r_B", "n", "median", "p99");
  std::vector<double> tE = {1.0, 1.1, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0, 1e9};
  for (size_t j = 0; j + 1 < tE.size(); ++j) {
    auto tk = d.Filter(Form("tang_b >= %g && tang_b < %g", tE[j], tE[j+1])).Take<double>("dpos_ulp");
    std::vector<double> u = *tk;
    if (u.size() < 20) continue;
    std::sort(u.begin(), u.end());
    P("  %7.3g .. %-8.3g %9zu %12.4g %12.4g\n", tE[j], tE[j+1], u.size(),
      u[u.size()/2], u[(size_t)(0.99*(u.size()-1))]);
  }

  P("\nsame pooled bins, but split by plane pair -- if tangency is the variable,\n"
    "each ROW must be flat across the columns:\n");
  P("  %-18s", "2*R_c / r_B");
  for (int i = 0; i < 7; ++i) P("%12s", Form("%.0f->%.0f", ra[i], rb[i]));
  P("\n");
  for (size_t j = 0; j + 1 < tE.size(); ++j) {
    P("  %7.3g .. %-8.3g", tE[j], tE[j+1]);
    for (int i = 0; i < 7; ++i) {
      auto r = med(d.Filter(Form("c.cfg == %d && tang_b >= %g && tang_b < %g", i, tE[j], tE[j+1])));
      if (r.first < 20) P("%12s", "-"); else P("%12.4g", r.second);
    }
    P("\n");
  }

  // The plot behind the tables: closure residual against the tangency ratio,
  // with one colour per plane pair. If tangency is the controlling variable the
  // seven curves lie on top of each other.
  TFile f(Form("%s.root", prefix), "RECREATE");
  auto dd = d.Define("log10_ulp", "dpos_ulp > 0 ? std::log10(dpos_ulp) : -3.0")
              .Define("log10_tang", "std::log10(tang_b)");
  auto all = dd.Profile1D({"tang_all",
      "closure vs tangency, all plane pairs;log_{10}( 2R_{c} / r_{B} );"
      "#LTlog_{10}( |#Deltapos| / ULP )#GT", 40, 0.0, 2.0}, "log10_tang", "log10_ulp");
  all->Write();
  for (int i = 0; i < 7; ++i) {
    auto h = dd.Filter(Form("c.cfg == %d", i))
               .Profile1D({Form("tang_cfg%d", i), Form("r %.0f #rightarrow %.0f cm;"
                          "log_{10}( 2R_{c} / r_{B} );#LTlog_{10}( |#Deltapos| / ULP )#GT",
                          ra[i], rb[i]), 40, 0.0, 2.0}, "log10_tang", "log10_ulp");
    h->Write();
  }
  auto pt = dd.Profile1D({"pt_all",
      "the same residual against p_{T} -- NOT the controlling variable;"
      "log_{10} p_{T} [GeV];#LTlog_{10}( |#Deltapos| / ULP )#GT", 40, -0.523, 2.301},
      "log10pt_", "log10_ulp");
  pt->Write();
  f.Close();
  P("\nwrote %s.root and %s.txt\n", prefix, prefix);
  fclose(log);
}
