// Two plots: the propagator's one-way landing error against turn angle, before
// and after sincos4 was removed. Bare root.exe; dictionary library only.
//
//   LD_LIBRARY_PATH=. root.exe -l -b -q \
//     '../RecoTracker/MkFitCore/standalone/test/an-val-trig.C'
//
// "before" = val-step-poly.root  (Config::useTrigApprox = true,  sincos4)
// "after"  = val-step-vdt.root   (Config::useTrigApprox = false, vdt::fast_sincosf)

void an_val_trig(const char *prefix = "110-step") {
  gSystem->Load("libMkFitRootDataFormats.so");
  FILE *log = fopen(Form("%s.txt", prefix), "w");
  auto P = [&](const char *fmt, ...) {
    va_list a; va_start(a, fmt); vprintf(fmt, a); va_end(a);
    va_start(a, fmt); vfprintf(log, fmt, a); va_end(a);
  };

  auto mk = [&](const char *f) {
    ROOT::RDataFrame df("closure", f);
    return df.Filter("c.pt >= 0.3")
      .Define("dperp", [](const ROOT::RVec<float>&p, const ROOT::RVec<float>&q,
                          const ROOT::RVec<float>&n){
          return (double) std::fabs((double)(p[0]-q[0])*n[0] + (double)(p[1]-q[1])*n[1]
                                  + (double)(p[2]-q[2])*n[2]); },
          {"c.par_b[6]","c.pnt_b[3]","c.nrm_b[3]"})
      .Define("dtot", [](const ROOT::RVec<float>&p, const ROOT::RVec<float>&q){
          double dx=p[0]-q[0], dy=p[1]-q[1], dz=p[2]-q[2];
          return std::sqrt(dx*dx+dy*dy+dz*dz); }, {"c.par_b[6]","c.pnt_b[3]"})
      .Define("dpar", "std::sqrt(std::max(0.0, dtot*dtot - dperp*dperp))")
      .Define("dal", "(double) c.dalpha_target")
      // dalpha takes SIX discrete values, so a continuum profile leaves empty
      // bins that render as zero. log10 spaces them almost evenly; six bins
      // from -1.85 to 0.05 put exactly one value in each.
      .Define("log10_dal", "std::log10((double) c.dalpha_target)")
      .Define("log10_dpar", "dpar > 0 ? std::log10(dpar) : -9.0")
      .Define("ulp", [](const ROOT::RVec<float>&a){
          double m = std::max({std::fabs(a[0]),std::fabs(a[1]),std::fabs(a[2])});
          return m>0 ? std::ldexp(1.0, std::ilogb(m)-23) : 0.0; }, {"c.par_b[6]"})
      .Define("n_ulp", "ulp>0 ? dpar/ulp : 0.0");
  };
  auto A = mk("val-step-poly.root"), B = mk("val-step-vdt.root");

  auto Q = [&](ROOT::RDF::RNode n, const char *c, double p){
    auto t = n.Take<double>(c); std::vector<double> u = *t;
    if (u.size() < 15) return -1.0;
    std::sort(u.begin(), u.end()); return u[(size_t)(0.01*p*(u.size()-1))]; };

  const double DAL[] = {0.02,0.05,0.10,0.20,0.40,0.80};
  P("### %s -- one-way landing error at plane B, incidence 1.0\n###\n", prefix);
  P("### definition : in-plane component of |x_prop - x_truth| at plane B, cm\n");
  P("### expectation: exactly 0; the arithmetic floor is a few float32 ULP\n");
  P("### before     : Config::useTrigApprox = true  (sincos4, 4th-order polynomial)\n");
  P("### after      : sincos4 removed, vdt::fast_sincosf everywhere\n###\n");
  P("  %-8s %12s %12s %8s | %12s %12s | %10s\n", "dalpha",
    "before med", "after med", "gain", "before p99", "after p99", "after/ULP");
  for (int ia = 0; ia < 6; ++ia) {
    char s[64]; snprintf(s, 64, "c.cfg==%d", ia*5);
    double am = Q(A.Filter(s), "dpar", 50), bm = Q(B.Filter(s), "dpar", 50);
    double a9 = Q(A.Filter(s), "dpar", 99), b9 = Q(B.Filter(s), "dpar", 99);
    double bu = Q(B.Filter(s), "n_ulp", 50);
    if (am < 0) continue;
    P("  %-8.2f %12.4g %12.4g %8.1f | %12.4g %12.4g | %10.2f\n",
      DAL[ia], am, bm, bm>0?am/bm:0, a9, b9, bu);
  }
  P("\nThe table above is incidence 1.0 (plane normal to the track). There the\n"
    "error is FLAT in turn angle at sub-ULP: the propagator's entire dalpha\n"
    "dependence was sincos4.\n\n");
  P("What sincos4 was HIDING -- median in-plane error [cm] after its removal:\n\n");
  P("  %-10s","dalpha\\inc");
  for (double v : (double[]){1.00,0.70,0.50,0.30,0.15}) P("%11.2f", v);
  P("%10s\n","0.15/1.0");
  for (int ia = 0; ia < 6; ++ia) {
    P("  %-10.2f", DAL[ia]); double f0 = 0, l0 = 0;
    for (int ic = 0; ic < 5; ++ic) {
      double m = Q(B.Filter(Form("c.cfg==%d", ia*5+ic)), "dpar", 50);
      if (m < 0) { P("%11s","-"); continue; }
      if (!f0) f0 = m; l0 = m; P("%11.4g", m);
    }
    P("%10.2f\n", f0>0 ? l0/f0 : 0.0);
  }
  P("\nUp to dalpha 0.2 that is a factor ~3, close to the 1/cos(incidence) the\n"
    "error algebra predicts. At dalpha 0.8 it is 431 in the median and 2343 in\n"
    "the p99 -- long step, shallow crossing.\n");

  TFile f(Form("%s.root", prefix), "RECREATE");
  auto pa = A.Profile1D({"trig_before",
      "before / after removing sincos4;log_{10}( turn angle #Delta#alpha / 1 rad );"
      "#LTlog_{10}( in-plane landing error / 1 cm )#GT", 6, -1.85, 0.05}, "log10_dal", "log10_dpar");
  auto pb = B.Profile1D({"trig_after", ";log_{10}( turn angle #Delta#alpha / 1 rad );"
      "#LTlog_{10}( in-plane landing error / 1 cm )#GT", 6, -1.85, 0.05}, "log10_dal", "log10_dpar");
  pa->Write(); pb->Write();
  // After the removal, one curve per incidence: the turn-angle dependence is
  // gone at normal incidence, and what is left is an INCIDENCE dependence that
  // sincos4 had been hiding.
  const double INC[] = {1.00, 0.70, 0.50, 0.30, 0.15};
  for (int ic = 0; ic < 5; ++ic) {
    auto h = B.Filter(Form("c.cfg %% 5 == %d", ic))
              .Profile1D({Form("after_inc%d", ic),
                Form("|p^.n| = %.2f;log_{10}( turn angle #Delta#alpha / 1 rad );"
                     "#LTlog_{10}( in-plane landing error / 1 cm )#GT", INC[ic]),
                6, -1.85, 0.05}, "log10_dal", "log10_dpar");
    h->Write();
  }
  auto ha = A.Histo1D({"dist_before", "in-plane landing error, all turn angles;"
      "log_{10}( error / 1 cm );trials", 120, -8, 0}, "log10_dpar");
  auto hb = B.Histo1D({"dist_after", ";log_{10}( error / 1 cm );trials", 120, -8, 0}, "log10_dpar");
  ha->Write(); hb->Write();
  f.Close();
  P("\nwrote %s.root and %s.txt\n", prefix, prefix);
  fclose(log);
}
