#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static const char* regname(int l){
  if (l<4) return "PixB";
  if (l<10) return "TBPS tilted (4-9)";
  if (l<16) return "TOB 2S flat (10-15)";
  if (l<28) return "FPix+";
  if (l<38) return "TEC+ (28-37)";
  if (l<50) return "FPix-";
  return "TEC- (50-59)";
}
void reg(const char *fn, const char *tag){
  TFile f(fn); TTree *t=(TTree*)f.Get("bkfit"); ValFitHit *v=nullptr; t->SetBranchAddress("h",&v);
  std::map<std::string,std::vector<double>> m;
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if (v->step==0||!std::isfinite(v->chi2)) continue;
    if (v->good_frac < 0.9999f) continue;
    if (v->mc_track_id != v->sim_label) continue;      // truth-matched hits only
    m[regname(v->layer)].push_back(v->chi2); }
  auto q=[](std::vector<double>&u,double p){ if(u.empty())return 0.0; return u[(size_t)(p/100.*(u.size()-1))]; };
  printf("\n%s -- PURE tracks, truth-matched hits. Per-hit chi2, ideal median 1.386\n", tag);
  printf("  %-22s %7s %9s %9s %9s %9s\n","region","hits","p25","median","p75","p90");
  for (auto &kv:m){ auto&u=kv.second; if(u.size()<25) continue; std::sort(u.begin(),u.end());
    printf("  %-22s %7zu %9.3g %9.3g %9.3g %9.3g\n",kv.first.c_str(),u.size(),q(u,25),q(u,50),q(u,75),q(u,90)); }
}
