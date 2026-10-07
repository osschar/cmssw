#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
#include <map>
// SCATTERING or WRECKED STATE? The killed hits are ~8x out in the precise
// direction with residual_z ~ 0, so the propagation lands correctly and the
// track really was elsewhere. Two readings:
//   (a) independent HARD SCATTERS -- each layer fails on its own, so failures
//       are Bernoulli and the per-candidate count is BINOMIAL;
//   (b) a WRECKED CANDIDATE STATE -- once wrong, wrong at every later layer, so
//       failures CLUSTER and the distribution is over-dispersed.
// Group the MC-matched hits that reached the Kalman by candidate and compare
// the observed spread of the kill count against the binomial it would have if
// independent. No new run, no new field.
static double q(std::vector<double> u,double p){ if(u.empty())return 0;
  std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
void an_killclust(const char *fn="val-surfq-1.root"){
  struct Cand { int n=0, k=0; };
  std::map<long long, Cand> C;
  TFile f(fn); TTree*t=(TTree*)f.Get("search");
  ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->mc_match||!h->had_kalman||h->chi2<-900) continue;
    if(h->global_seed<0||h->event<0) continue;
    const long long key = (long long)h->event*1000000LL + h->global_seed;
    auto &c = C[key]; ++c.n; if(h->chi2>=30.f) ++c.k; }

  // overall rate, and the binomial prediction for candidates with n hits
  long N=0,K=0; for(auto&kv:C){ N+=kv.second.n; K+=kv.second.k; }
  const double p = (double)K/N;
  printf("\n===== DO THE chi2 FAILURES CLUSTER WITHIN A CANDIDATE? =====\n");
  printf("%zu candidates, %ld MC-matched hits reaching the Kalman, %ld killed\n",
         C.size(), N, K);
  printf("overall kill rate p = %.4f\n\n", p);
  printf("If failures are INDEPENDENT (hard scatters), the number killed among a\n");
  printf("candidate's n hits is Binomial(n, p) and var/mean = 1-p ~ %.2f.\n", 1-p);
  printf("A WRECKED STATE clusters them: var/mean well above 1, and far more\n");
  printf("all-killed candidates than binomial allows.\n\n");
  printf("  %-6s %9s | %8s %8s | %10s %10s | %9s %9s\n","n hits","cands",
         "mean k","var k","var/mean","binom exp","all-kill","binom");
  for(int n=2;n<=6;++n){
    std::vector<double> ks;
    for(auto&kv:C) if(kv.second.n==n) ks.push_back(kv.second.k);
    if(ks.size()<200) continue;
    double m=0; for(double v:ks) m+=v; m/=ks.size();
    double va=0; for(double v:ks) va+=(v-m)*(v-m); va/=ks.size();
    long all=0; for(double v:ks) if((int)v==n) ++all;
    printf("  %-6d %9zu | %8.4f %8.4f | %10.3f %10.3f | %8.3f%% %8.4f%%\n",
           n, ks.size(), m, va, m>0?va/m:0.0, 1-p,
           100.0*all/ks.size(), 100.0*std::pow(p,n));
  }
  printf("\n  'all-kill' vs 'binom' is the cleanest discriminator: a wrecked state\n");
  printf("  produces candidates whose every hit fails, which independence forbids.\n\n");
}
