#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
#include <map>
// The chi2 failures CLUSTER within a candidate (all-killed candidates 300-25000x
// more common than binomial), which rules out independent hard scatters. Is the
// wrong state inherited from a BAD SEED, or acquired during the search?
//
// seed_good_frac = of the seed's valid hits, the fraction from its best-matching
// sim track. A pure seed (1.0) that still fails wholesale acquired the problem;
// an impure seed was never following that track to begin with.
void an_purity(const char *fn="val-purity.root"){
  const int NG=5; const double gb[NG+1]={0.0,0.6,0.8,0.95,0.999,1.001};
  const char* gn[NG]={"<0.6","0.6-0.8","0.8-0.95","0.95-1.0","1.0 (pure)"};
  long n[NG]={}, k[NG]={};
  struct C{int n=0,k=0;float gf=-1;};
  std::map<long long,C> cand;
  TFile f(fn); TTree*t=(TTree*)f.Get("search");
  ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->mc_match||!h->had_kalman||h->chi2<-900||h->seed_good_frac<0) continue;
    int g=-1; for(int b=0;b<NG;++b) if(h->seed_good_frac>=gb[b]&&h->seed_good_frac<gb[b+1]){g=b;break;}
    if(g<0) continue;
    ++n[g]; if(h->chi2>=30.f) ++k[g];
    const long long key=(long long)h->event*1000000LL+h->global_seed;
    auto &c=cand[key]; ++c.n; if(h->chi2>=30.f) ++c.k; c.gf=h->seed_good_frac; }

  printf("\n===== IS THE WRONG STATE INHERITED FROM A BAD SEED? =====\n");
  printf("seed purity = fraction of the seed's valid hits from its best sim track.\n\n");
  printf("  %-12s %10s %10s %9s\n","seed purity","n hits","killed","kill rate");
  for(int g=0;g<NG;++g) if(n[g]>=200)
    printf("  %-12s %10ld %10ld %8.2f%%\n", gn[g], n[g], k[g], 100.0*k[g]/n[g]);

  printf("\n  and the ALL-KILLED candidates (every MC-matched hit over the cut),\n");
  printf("  which is where the clustering lives -- by seed purity:\n");
  printf("  %-12s %10s %12s %10s\n","seed purity","cands(n>=3)","all-killed","");
  long tot[NG]={}, allk[NG]={};
  for(auto&kv:cand){ if(kv.second.n<3||kv.second.gf<0) continue;
    int g=-1; for(int b=0;b<NG;++b) if(kv.second.gf>=gb[b]&&kv.second.gf<gb[b+1]){g=b;break;}
    if(g<0) continue; ++tot[g]; if(kv.second.k==kv.second.n) ++allk[g]; }
  for(int g=0;g<NG;++g) if(tot[g]>=50)
    printf("  %-12s %10ld %11ld %9.2f%%\n", gn[g], tot[g], allk[g], 100.0*allk[g]/tot[g]);
  printf("\n  If the kill rate and the all-killed fraction are FLAT in purity, the\n");
  printf("  wrong state is ACQUIRED during the search. If they rise steeply as\n");
  printf("  purity falls, it is INHERITED from the seed.\n\n");
}
