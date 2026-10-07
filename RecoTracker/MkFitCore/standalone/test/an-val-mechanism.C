// What controls the closure residual: the tangency scan, the pixel-radius
// reproduction, the track-perpendicular plane control, and what survives it.
// Bare root.exe; only the dictionary library is needed.
//
//   LD_LIBRARY_PATH=. root.exe -l -b -q \
//     '../RecoTracker/MkFitCore/standalone/test/an-val-mechanism.C("val-closure.root")'
//
// cfg map, from TTree "cfgs":  0-6 outer radii, radial normals
//                              7-13 outer radii, plane PERPENDICULAR to the track
//                             14-17 pixel radii, radial normals
//                             18-21 pixel radii, perpendicular

void an_val_mechanism(const char *in_file = "val-closure.root",
                      const char *prefix = "130-mechanism",
                      double pt_min = 0.0) {
  gSystem->Load("libMkFitRootDataFormats.so");
  FILE *log = fopen(Form("%s.txt", prefix), "w");
  auto P = [&](const char *fmt, ...) {
    va_list a; va_start(a, fmt); vprintf(fmt, a); va_end(a);
    va_start(a, fmt); vfprintf(log, fmt, a); va_end(a);
  };

  ROOT::RDataFrame df("closure", in_file);
  auto d = df.Define("dpos", [](const ROOT::RVec<float>&a, const ROOT::RVec<float>&c){
              double dx=c[0]-a[0], dy=c[1]-a[1], dz=c[2]-a[2];
              return std::sqrt(dx*dx+dy*dy+dz*dz); }, {"c.par_a[6]","c.par_c[6]"})
           .Define("ulp", [](const ROOT::RVec<float>&a){
              double m = std::max({std::fabs(a[0]),std::fabs(a[1]),std::fabs(a[2])});
              return m>0 ? std::ldexp(1.0, std::ilogb(m)-23) : 0.0; }, {"c.par_a[6]"})
           .Define("r", "ulp>0 ? dpos/ulp : 0.0")
           .Filter("(bool) std::isfinite(r)")
           .Filter(Form("c.pt >= %g", pt_min))
           // R_c [cm] = pT / (0.01*sol*B) = 87.515 * pT at B = 3.8112 T
           .Define("tang_b", "2.0 * 87.515 * c.pt / c.r_b")
           .Define("log10_tang", "std::log10(tang_b)")
           .Define("log10_r", "r > 0 ? std::log10(r) : -1.0")
           .Define("dal", "(double) c.dalpha_ab");

  auto med = [&](ROOT::RDF::RNode n){
    auto t = n.Take<double>("r"); std::vector<double> u = *t;
    if (u.size() < 20) return std::array<double,3>{0,0,0};
    std::sort(u.begin(), u.end());
    return std::array<double,3>{(double)u.size(), u[u.size()/2], u[(size_t)(0.99*(u.size()-1))]};
  };

  const char *OUT_R = "c.cfg < 7", *OUT_P = "c.cfg >= 7 && c.cfg < 14";
  const char *PIX_R = "c.cfg >= 14 && c.cfg < 18", *PIX_P = "c.cfg >= 18";

  P("### %s -- what controls the round-trip closure residual\n###\n", prefix);
  P("### quantity   : median |d pos| / (one float32 ULP), dimensionless\n");
  P("### 2*R_c/r_B  : turning-circle DIAMETER over the radius of plane B.\n");
  P("###              1 = the track grazes plane B and turns back.\n");
  P("### pT cut     : >= %g GeV (0 = none; the pixel pairs need pT < 0.3 to\n", pt_min);
  P("###              reach tangency at all, since 2R_c/r_B = 1 at pT = r_B/175)\n###\n");

  std::vector<double> e = {1.0, 1.1, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0, 1e9};

  P("(1) THE POOLING TEST -- outer radii (30-110 cm) against PIXEL radii (4-16 cm).\n");
  P("    Radial normals. pT and coordinates differ by ~10x between the columns;\n");
  P("    if 2R_c/r_B is the controlling variable they must still agree.\n\n");
  P("  %-16s %8s %12s | %8s %12s\n", "2R_c/r_B", "n outer", "med outer", "n pixel", "med pixel");
  for (size_t j = 0; j + 1 < e.size(); ++j) {
    auto A = med(d.Filter(Form("%s && tang_b>=%g && tang_b<%g", OUT_R, e[j], e[j+1])));
    auto B = med(d.Filter(Form("%s && tang_b>=%g && tang_b<%g", PIX_R, e[j], e[j+1])));
    if (A[0] < 20 && B[0] < 20) continue;
    P("  %5.3g .. %-8.3g %8.0f %12.4g | %8.0f %12.4g\n", e[j], e[j+1], A[0], A[1], B[0], B[1]);
  }

  P("\n(2) THE PLANE CONTROL -- same truth crossings, but the plane is turned to be\n");
  P("    PERPENDICULAR TO THE TRACK. The plane equation n.(x(s)-pnt)=0 then has\n");
  P("    derivative n.dx/ds = 1, the best conditioning available.\n\n");
  P("  %-16s %12s %12s | %12s %12s\n", "2R_c/r_B", "radial", "perp", "pixel rad", "pixel perp");
  for (size_t j = 0; j + 1 < e.size(); ++j) {
    auto A = med(d.Filter(Form("%s && tang_b>=%g && tang_b<%g", OUT_R, e[j], e[j+1])));
    auto B = med(d.Filter(Form("%s && tang_b>=%g && tang_b<%g", OUT_P, e[j], e[j+1])));
    auto C = med(d.Filter(Form("%s && tang_b>=%g && tang_b<%g", PIX_R, e[j], e[j+1])));
    auto E = med(d.Filter(Form("%s && tang_b>=%g && tang_b<%g", PIX_P, e[j], e[j+1])));
    if (A[0] < 20 && B[0] < 20) continue;
    P("  %5.3g .. %-8.3g %12.4g %12.4g | %12.4g %12.4g\n", e[j], e[j+1], A[1], B[1], C[1], E[1]);
  }

  P("\n(3) WHAT SURVIVES IT -- banded in |dalpha|, the helix angle turned from A to B.\n");
  P("    Both radius scales pooled in each column.\n\n");
  P("  %-16s %8s %12s | %8s %12s\n", "|dalpha| [rad]", "n radial", "med radial", "n perp", "med perp");
  std::vector<double> a = {0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 2.5, 3.2};
  for (size_t j = 0; j + 1 < a.size(); ++j) {
    auto A = med(d.Filter(Form("(%s || %s) && dal>=%g && dal<%g", OUT_R, PIX_R, a[j], a[j+1])));
    auto B = med(d.Filter(Form("(%s || %s) && dal>=%g && dal<%g", OUT_P, PIX_P, a[j], a[j+1])));
    if (A[0] < 20 && B[0] < 20) continue;
    P("  %5.3g .. %-8.3g %8.0f %12.4g | %8.0f %12.4g\n", a[j], a[j+1], A[0], A[1], B[0], B[1]);
  }

  P("\n(4) non-finite residuals: %lld of %lld trials in the file.\n",
    (long long) *df.Define("dp", [](const ROOT::RVec<float>&x, const ROOT::RVec<float>&y){
        double dx=y[0]-x[0],dy=y[1]-x[1],dz=y[2]-x[2];
        return std::sqrt(dx*dx+dy*dy+dz*dz); }, {"c.par_a[6]","c.par_c[6]"})
      .Filter("!(bool) std::isfinite(dp)").Count(),
    (long long) *df.Count());

  TFile f(Form("%s.root", prefix), "RECREATE");
  auto prof = [&](const char *nm, const char *sel, const char *xc, const char *ti,
                  int nb, double lo, double hi) {
    auto h = d.Filter(sel).Profile1D({nm, ti, nb, lo, hi}, xc, "log10_r");
    h->Write();
  };
  const char *TT = ";log_{10}( 2R_{c} / r_{B} );#LTlog_{10}( |#Deltapos| / ULP )#GT";
  prof("tang_outer_rad", OUT_R, "log10_tang", Form("outer radii, radial normals%s", TT), 40, 0, 2);
  prof("tang_outer_perp", OUT_P, "log10_tang", Form("outer radii, track-perpendicular%s", TT), 40, 0, 2);
  prof("tang_pixel_rad", PIX_R, "log10_tang", Form("pixel radii, radial normals%s", TT), 40, 0, 2);
  prof("tang_pixel_perp", PIX_P, "log10_tang", Form("pixel radii, track-perpendicular%s", TT), 40, 0, 2);
  const char *TD = ";|#Delta#alpha| from A to B [rad];#LTlog_{10}( |#Deltapos| / ULP )#GT";
  prof("dal_rad", Form("%s || %s", OUT_R, PIX_R), "dal", Form("radial normals%s", TD), 40, 0, 3.2);
  prof("dal_perp", Form("%s || %s", OUT_P, PIX_P), "dal", Form("track-perpendicular%s", TD), 40, 0, 3.2);
  f.Close();
  P("\nwrote %s.root (6 profiles) and %s.txt\n", prefix, prefix);
  fclose(log);
}
