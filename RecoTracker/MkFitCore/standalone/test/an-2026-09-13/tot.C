#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
void tot(const char *fn, float sc) {
  TFile f(fn); TTree *t=(TTree*)f.Get("bkfit"); ValFitHit *v=nullptr; t->SetBranchAddress("h",&v);
  std::map<std::pair<int,int>,std::pair<double,int>> trk; std::map<std::pair<int,int>,float> gf;
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if (v->step==0||!std::isfinite(v->chi2)) continue;
    auto k=std::make_pair(v->event,v->label); trk[k].first+=v->chi2; trk[k].second++; gf[k]=v->good_frac; }
  std::vector<double> pure, all;
  for (auto &kv:trk){ double r=kv.second.first/(2.0*kv.second.second);
    all.push_back(r); if (gf[kv.first]>=0.9999f) pure.push_back(r); }
  auto q=[](std::vector<double>&u,double p){ if(u.empty())return 0.0; return u[(size_t)(p/100.*(u.size()-1))]; };
  std::sort(pure.begin(),pure.end()); std::sort(all.begin(),all.end());
  printf("  %-9g %7zu %9.3g %9.3g %9.3g %9.3g | %9.3g %9.3g %9.3g\n", sc,
    pure.size(), q(pure,25), q(pure,50), q(pure,75), q(pure,90), q(all,50), q(all,90), all.back());
}
