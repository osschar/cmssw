// What turn angle does a real track actually present between consecutive hits?
// Bare root.exe; dictionary library only.
//
//   LD_LIBRARY_PATH=. root.exe -l -b -q \
//     '../RecoTracker/MkFitCore/standalone/test/an-val-dalpha.C("val-dalpha-w1.root")'
//
// This is the weighting function for every scan made against turn angle. A
// synthetic scan chooses Delta-alpha; a detector does not.

void an_val_dalpha(const char *in_file = "val-dalpha-w1.root",
                   const char *prefix = "130-dalpha") {
  gSystem->Load("libMkFitRootDataFormats.so");
  FILE *log = fopen(Form("%s.txt", prefix), "w");
  auto P = [&](const char *fmt, ...) { va_list a; va_start(a,fmt); vprintf(fmt,a); va_end(a);
    va_start(a,fmt); vfprintf(log,fmt,a); va_end(a); };

  ROOT::RDataFrame df("steps", in_file);
  auto d = df.Define("Rc", "(double) s.pt / (0.01 * 0.299792458 * 3.8112)")
             // chord / turning-circle diameter. >= 1 cannot lie on one circle:
             // bad hit ordering, a secondary, or a looper.
             .Define("ratio", "(double) s.d_perp / (2.0 * Rc)")
             .Define("da", "(double) s.dalpha")
             .Define("ae", "std::fabs((double) s.eta)");
  auto Q=[&](ROOT::RDF::RNode n,const char*c,double p){auto t=n.Take<double>(c);
    std::vector<double>u=*t; if(u.size()<20)return -1.0;
    std::sort(u.begin(),u.end()); return u[(size_t)(0.01*p*(u.size()-1))];};

  const long n_all = *d.Count(), n_bad = *d.Filter("ratio >= 0.999").Count();
  P("### %s -- turn angle per propagation step, REAL tracks\n###\n", prefix);
  P("### definition : the helix angle turned between two consecutive hits of a\n");
  P("###              track, from the chord and the track's own R_c = pT/(0.3B)\n");
  P("### unit       : rad\n");
  P("### population : sim tracks with their rec hits, 20 events, %ld steps;\n", n_all);
  P("###              %ld (%.2f%%) discarded as geometrically impossible\n", n_bad, 100.0*n_bad/n_all);
  P("###\n");
  P("  %-10s %10s %9s %9s %9s %9s | %9s %9s\n","pT cut","n","median","p90","p99","max",">0.4 rad",">0.8 rad");
  for (double c : {0.5,0.7,1.0,2.0}) {
    auto n = d.Filter(Form("s.pt >= %g && ratio < 0.999", c));
    const long nn = *n.Count();
    P("  >= %-7.1f %10ld %9.4f %9.4f %9.4f %9.4f | %8.3f%% %8.3f%%\n", c, nn,
      Q(n,"da",50),Q(n,"da",90),Q(n,"da",99),Q(n,"da",100),
      100.0**n.Filter("da>0.4").Count()/nn, 100.0**n.Filter("da>0.8").Count()/nn);
  }
  P("\n  where the large steps are -- pT >= 0.5, clean:\n");
  P("  %-12s %10s %10s %10s %11s\n","|eta|","n","median","p99","frac > 0.4");
  std::vector<double> e={0,0.8,1.2,1.7,2.2,3.0,5.0};
  for (size_t j=0;j+1<e.size();++j){
    auto n=d.Filter(Form("s.pt>=0.5 && ratio<0.999 && ae>=%g && ae<%g",e[j],e[j+1]));
    if(*n.Count()<50) continue;
    P("  %4.2g .. %-6.2g %10lld %10.4f %10.4f %10.3f%%\n", e[j],e[j+1],(long long)*n.Count(),
      Q(n,"da",50),Q(n,"da",99),100.0**n.Filter("da>0.4").Count()/ *n.Count());
  }
  P("\n  Large steps are commonest in the BARREL, not the transition: low-pT tracks\n"
    "  curving between widely spaced barrel layers. So if the transition region is\n"
    "  where trouble lives, turn angle is not the variable that explains it.\n");

  TFile f(Form("%s.root", prefix), "RECREATE");
  for (auto [cut, nm] : {std::pair<double,const char*>{0.5,"dalpha_pt05"},
                         {0.7,"dalpha_pt07"}, {2.0,"dalpha_pt20"}}) {
    auto h = d.Filter(Form("s.pt >= %g && ratio < 0.999", cut))
              .Define("l10","da>0 ? std::log10(da) : -4.0")
              .Histo1D({nm, Form("p_{T} #geq %.1f GeV;log_{10}( #Delta#alpha / 1 rad );steps", cut),
                        90, -4, 0.5}, "l10");
    h->Write();
  }
  f.Close();
  P("\nwrote %s.root and %s.txt\n", prefix, prefix);
  fclose(log);
}
