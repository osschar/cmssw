#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static double qq(std::vector<double>u,double p){ if(u.empty())return 0; std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
// Phase-2 layer map (CLAUDE.md): 0-3 PixB, 4-9 TOB PS, 10-15 TOB 2S,
// 16-27 FPix+, 28-37 TEC+, 38-49 FPix-, 50-59 TEC-.
static const char* grp(int l){
  if(l>=16&&l<24) return "TFPX  16-23";
  if(l>=24&&l<28) return "TEPX  24-27";
  if(l>=38&&l<46) return "TFPX- 38-45";
  if(l>=46&&l<50) return "TEPX- 46-49";
  if(l>=28&&l<38) return "TEC+  28-37";
  if(l>=50)       return "TEC-  50-59";
  return nullptr; }
void an_fpix(){
  const char*src[2]={"val-search-mat1.root","val-search-mat0.root"};
  const char*sn[2]={"ON ","OFF"};
  printf("\nMATERIAL IN THE PIXEL ENDCAP. MC-matched hits, search chi2 (2-D in the module\n");
  printf("plane; median of chi2_2 = 1.386). 'sig short' = sqrt(med/1.386) -- how much too\n");
  printf("small the combined sigma is. 100 events, T5 seeds.\n\n");
  printf("  %-12s %4s %8s | %8s %8s | %9s\n","group","mat","n","chi2 med","p90","sig short");
  std::map<std::string,std::array<std::vector<double>,2>> G;
  for(int m=0;m<2;++m){ TFile f(src[m]); TTree*t=(TTree*)f.Get("search");
    ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
      if(!h->mc_match||!h->had_kalman||h->chi2<=0) continue;
      const char*g=grp(h->layer); if(!g) continue;
      G[g][m].push_back(h->chi2); } }
  for(auto&kv:G){ for(int m=0;m<2;++m){ auto&u=kv.second[m]; if(u.size()<40) continue;
      double md=qq(u,50);
      printf("  %-12s %4s %8zu | %8.3f %8.2f | %8.2fx\n",m?"":kv.first.c_str(),sn[m],u.size(),
             md,qq(u,90),std::sqrt(md/1.38629436)); } printf("\n"); }
}
