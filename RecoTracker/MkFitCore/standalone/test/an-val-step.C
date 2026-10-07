// One step there and back: propagation accuracy on its two real controls --
// how far the step turns, and how squarely it meets the target plane.
//
//   LD_LIBRARY_PATH=. root.exe -l -b -q \
//     '../RecoTracker/MkFitCore/standalone/test/an-val-step.C("val-step.root")'
//
// WHAT IS MEASURED.  Material off, uniform field: the propagation is a
// diffeomorphism, so O -> A -> B -> A' must return the input EXACTLY. Plane A
// sits at the production point with its normal along the momentum, so the
// return leg is perfectly conditioned and contributes no error of its own;
// everything below is a property of the outbound step A -> B.
//
// Unit: cm, and dimensionless once divided by one float32 ULP -- the gap
// between adjacent representable floats at that coordinate magnitude,
// 2^floor(log2|x|) * 2^-23. Expectation: a few ULP. Radius is never a control.

void an_val_step(const char *in_file = "val-step.root",
                 const char *prefix = "200-step", double pt_min = 0.3) {
  gSystem->Load("libMkFitRootDataFormats.so");
  FILE *log = fopen(Form("%s.txt", prefix), "w");
  auto P = [&](const char *fmt, ...) {
    va_list a; va_start(a, fmt); vprintf(fmt, a); va_end(a);
    va_start(a, fmt); vfprintf(log, fmt, a); va_end(a);
  };

  ROOT::RDataFrame df("closure", in_file);
  auto d = df.Filter(Form("c.pt >= %g", pt_min))
    .Define("dpos", [](const ROOT::RVec<float>&a, const ROOT::RVec<float>&c){
        double dx=c[0]-a[0], dy=c[1]-a[1], dz=c[2]-a[2];
        return std::sqrt(dx*dx+dy*dy+dz*dz); }, {"c.par_a[6]","c.par_c[6]"})
    // Compare against the EXACT crossing, kept in double: plane B is built to
    // pass through the truth helix's own crossing, so this is an absolute
    // reference. Storing it as float would make the measurement report the
    // reference's rounding -- at r ~ 40 cm that is 3.8e-6 cm, the same size as
    // the error.
    .Define("e_ref", [](const ROOT::RVec<float>&p, const ROOT::RVec<double>&q){
        double dx=(double)p[0]-q[0], dy=(double)p[1]-q[1], dz=(double)p[2]-q[2];
        return std::sqrt(dx*dx+dy*dy+dz*dz); }, {"c.par_b[6]","c.ref_pnt_b[3]"})
    .Define("ulp", [](const ROOT::RVec<double>&a){
        double m = std::max({std::fabs(a[0]),std::fabs(a[1]),std::fabs(a[2])});
        return m>0 ? std::ldexp(1.0, std::ilogb(m)-23) : 0.0; }, {"c.ref_pnt_b[3]"})
    .Define("n_ref", "ulp>0 ? e_ref/ulp : 0.0")
    .Define("r", "ulp>0 ? dpos/ulp : 0.0")
    // Landing at B, split: perpendicular to the plane is what the solver
    // solves for; in-plane is where along the plane it ended up.
    .Define("dperp", [](const ROOT::RVec<float>&p, const ROOT::RVec<float>&q,
                        const ROOT::RVec<float>&n){
        return (double) std::fabs((double)(p[0]-q[0])*n[0] + (double)(p[1]-q[1])*n[1]
                                + (double)(p[2]-q[2])*n[2]); },
        {"c.par_b[6]","c.pnt_b[3]","c.nrm_b[3]"})
    .Define("dtot", [](const ROOT::RVec<float>&p, const ROOT::RVec<float>&q){
        double dx=p[0]-q[0], dy=p[1]-q[1], dz=p[2]-q[2];
        return std::sqrt(dx*dx+dy*dy+dz*dz); }, {"c.par_b[6]","c.pnt_b[3]"})
    .Define("dpar", "std::sqrt(std::max(0.0, dtot*dtot - dperp*dperp))")
    // The conditioning leg is a ZERO-LENGTH step onto a plane the state is
    // already on. That is the same degenerate step the backward fit takes on
    // the first hit of essentially every track, so it is worth watching.
    .Define("condmiss", [](const ROOT::RVec<float>&p, const ROOT::RVec<float>&q){
        double dx=p[0]-q[0], dy=p[1]-q[1], dz=p[2]-q[2];
        return std::sqrt(dx*dx+dy*dy+dz*dz); }, {"c.par_a[6]","c.pnt_a[3]"})
    .Define("dal", "(double) c.dalpha_target")
    .Define("inc", "(double) c.cosinc_target")
    .Filter("(bool) std::isfinite(r)");

  auto med = [&](ROOT::RDF::RNode n, const char *col){
    auto t = n.Take<double>(col); std::vector<double> u = *t;
    if (u.size() < 15) return std::array<double,4>{0,0,0,0};
    std::sort(u.begin(), u.end());
    auto q=[&](double p){ return u[(size_t)(0.01*p*(u.size()-1))]; };
    return std::array<double,4>{(double)u.size(), q(50), q(99), u.back()};
  };

  const double DAL[] = {0.02,0.05,0.10,0.20,0.40,0.80};
  const double INC[] = {1.00,0.70,0.50,0.30,0.15};
  const int NI = 5;
  // Select by cfg index, never by float equality on a control value: the stored
  // value is a float, so (double) 0.02f != 0.02 and every such Filter is empty.
  auto CFG = [&](int ia, int ic){ return ia * NI + ic; };

  P("### %s -- one step there and back\n###\n", prefix);
  P("### model      : the propagator solves for a path length s with n.(x(s)-p)=0.\n");
  P("###              It does NOT propagate to a radius. The two controls are how\n");
  P("###              far the step turns (dalpha) and how squarely it meets the\n");
  P("###              target plane (|p^.n|). Radius is an output.\n");
  P("### definition : |d pos| = |r(A') - r(A)| after O->A->B->A'\n");
  P("### unit       : cm, and dimensionless in units of one float32 ULP at B\n");
  P("### expectation: exactly 0; the arithmetic floor is a few ULP\n");
  P("### population : pT >= %g GeV, %lld trials; configurations landing outside\n", pt_min, (long long)*d.Count());
  P("###              r < 120 / |z| < 300 cm are not generated\n###\n");

  P("(1) median |d pos| / ULP.  rows = turn angle, columns = incidence at B\n\n");
  P("  %-12s", "dalpha \\ inc");
  for (double v : INC) P("%11.2f", v);
  P("%11s\n", "ratio");
  for (int ia = 0; ia < 6; ++ia) {
    P("  %-12.2f", DAL[ia]);
    double first = 0, last = 0;
    for (int ic = 0; ic < NI; ++ic) {
      auto m = med(d.Filter(Form("c.cfg==%d", CFG(ia, ic))), "r");
      if (!m[0]) { P("%11s", "-"); continue; }
      if (!first) first = m[1];
      last = m[1];
      P("%11.4g", m[1]);
    }
    P("%11.3g\n", first > 0 ? last/first : 0.0);
  }
  P("\n  last column: how much incidence 0.15 costs relative to incidence 1,\n");
  P("  at that turn angle. If the two controls factorise it is constant.\n");

  P("\n(2) does it land ON plane B?  perpendicular vs in-plane, cm\n\n");
  P("  %-10s %-8s %8s %12s %12s %12s %12s\n", "dalpha","inc","n","med perp","p99 perp","med in-pl","p99 in-pl");
  for (int ia = 0; ia < 6; ++ia)
    for (int ic : {0, 3, 4}) {
      auto n = d.Filter(Form("c.cfg==%d", CFG(ia, ic)));
      auto mp = med(n, "dperp"), ml = med(n, "dpar");
      if (!mp[0]) continue;
      P("  %-10.2f %-8.2f %8.0f %12.4g %12.4g %12.4g %12.4g\n",
        DAL[ia], INC[ic], mp[0], mp[1], mp[2], ml[1], ml[2]);
    }

  P("\n(3) the zero-length conditioning leg O->A (the state is already on plane A)\n\n");
  auto cm = med(d, "condmiss");
  P("  |pos_A - plane-A point| : med %.4g  p99 %.4g  max %.4g cm\n", cm[1], cm[2], cm[3]);
  P("  fail_oa set : %lld of %lld\n",
    (long long)*d.Filter("c.fail_oa != 0").Count(), (long long)*d.Count());
  P("  any leg flagged : %lld\n",
    (long long)*d.Filter("c.fail_oa || c.fail_ab || c.fail_ba").Count());

  P("\n(4) where these steps actually are, for orientation (not a control)\n\n");
  P("  %-10s %10s %10s %10s %10s\n", "dalpha", "med r_B", "p99 r_B", "med |z_B|", "p99 |z_B|");
  for (int ia = 0; ia < 6; ++ia) {
    auto n = d.Filter(Form("c.cfg >= %d && c.cfg < %d", ia*NI, (ia+1)*NI))
              .Define("az","std::fabs((double)c.z_b)").Define("rb","(double)c.r_b");
    auto mr = med(n, "rb"), mz = med(n, "az");
    if (!mr[0]) continue;
    P("  %-10.2f %10.4g %10.4g %10.4g %10.4g\n", DAL[ia], mr[1], mr[2], mz[1], mz[2]);
  }

  TFile f(Form("%s.root", prefix), "RECREATE");
  auto dl = d.Define("log10_r", "r > 0 ? std::log10(r) : -1.0");
  for (int i = 0; i < 5; ++i) {
    auto h = dl.Filter(Form("c.cfg %% 5 == %d", i))
               .Profile1D({Form("dal_inc%d", i),
                 Form("|p^.n| = %.2f;turn angle #Delta#alpha  [rad];"
                      "#LTlog_{10}( |#Deltapos| / ULP )#GT", INC[i]), 40, 0, 0.85},
                 "dal", "log10_r");
    h->Write();
  }
  for (int i = 0; i < 6; ++i) {
    auto h = dl.Filter(Form("c.cfg / 5 == %d", i))
               .Profile1D({Form("inc_dal%d", i),
                 Form("#Delta#alpha = %.2f rad;incidence |p^.n| at B;"
                      "#LTlog_{10}( |#Deltapos| / ULP )#GT", DAL[i]), 20, 0, 1.05},
                 "inc", "log10_r");
    h->Write();
  }
  f.Close();
  P("\nwrote %s.root and %s.txt\n", prefix, prefix);
  fclose(log);
}
