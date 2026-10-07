#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
void prof(const char *fn) {
  TFile f(fn); TTree *t = (TTree*)f.Get("bkfit");
  ValFitHit *v = nullptr; t->SetBranchAddress("h", &v);
  std::map<int,std::vector<double>> by_step, by_lay;
  std::vector<double> pix, strp;
  for (Long64_t i=0;i<t->GetEntries();++i) { t->GetEntry(i);
    if (v->step==0 || !std::isfinite(v->chi2)) continue;
    if (v->good_frac < 0.9999f) continue;              // PURE tracks only
    by_step[std::min(v->step,9)].push_back(v->chi2);
    (v->layer < 4 || (v->layer>=16 && v->layer<28) || (v->layer>=38 && v->layer<50) ? pix : strp).push_back(v->chi2);
    by_lay[v->layer].push_back(v->chi2);
  }
  auto q=[](std::vector<double>&u,double p){ if(u.empty())return 0.0; size_t k=(size_t)(p/100.*(u.size()-1)); return u[k]; };
  printf("\nPURE tracks only. Per-hit chi2; ideal median for 2 d.o.f. = 1.386\n");
  printf("  %-8s %7s %9s %9s %9s\n","step","hits","p25","median","p75");
  for (auto &kv:by_step){ auto&u=kv.second; std::sort(u.begin(),u.end());
    printf("  %-8d %7zu %9.3g %9.3g %9.3g\n",kv.first,u.size(),q(u,25),q(u,50),q(u,75)); }
  std::sort(pix.begin(),pix.end()); std::sort(strp.begin(),strp.end());
  printf("\n  %-8s %7s %9s %9s %9s\n","type","hits","p25","median","p75");
  printf("  %-8s %7zu %9.3g %9.3g %9.3g\n","pixel",pix.size(),q(pix,25),q(pix,50),q(pix,75));
  printf("  %-8s %7zu %9.3g %9.3g %9.3g\n","strip",strp.size(),q(strp,25),q(strp,50),q(strp,75));
  printf("\n  %-8s %7s %9s   (layers with >=25 hits)\n","layer","hits","median");
  for (auto &kv:by_lay){ if(kv.second.size()<25) continue; auto&u=kv.second; std::sort(u.begin(),u.end());
    printf("  %-8d %7zu %9.3g\n",kv.first,u.size(),q(u,50)); }
}
