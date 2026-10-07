#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static double qq(std::vector<double>u,double p){ if(u.empty())return 0; std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
struct SR{int step=-1,layer=-1,nsim=-1; bool won=false; double c2mc=1e30,c2b=1e30; bool kal=false;};
static void one(const char*fn,const char*tag){
  TFile f(fn); TTree *t=(TTree*)f.Get("search"); ValSearchHit *h=nullptr; t->SetBranchAddress("h",&h);
  std::map<std::pair<int,int>,SR> S;                       // (event, search_id)
  std::map<std::pair<int,int>,std::pair<int,int>> srch2seed; // -> (global_seed, sim_label)
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(h->search_id<0) continue;
    auto k=std::make_pair(h->event,h->search_id);
    SR&s=S[k]; s.step=h->step; s.layer=h->layer;
    if(h->n_sim_hits_in_layer>s.nsim) s.nsim=h->n_sim_hits_in_layer;
    if(h->global_seed>=0) srch2seed[k]={h->global_seed,h->sim_label};
    if(h->had_kalman&&h->chi2>-900){ if(h->chi2<s.c2b)s.c2b=h->chi2;
      if(h->mc_match){s.kal=true; if(h->chi2<s.c2mc)s.c2mc=h->chi2;} } }
  for(auto&kv:S) kv.second.won = kv.second.kal && kv.second.c2mc<=kv.second.c2b;
  // ---- regroup by seed
  struct SD{int maxstep=-1; int navail=0,nwon=0; int firstloss=99; bool endcap=false;};
  std::map<std::pair<int,int>,SD> D;
  for(auto&kv:S){ auto it=srch2seed.find(kv.first); if(it==srch2seed.end()) continue;
    if(it->second.second<0) continue;                       // no truth for this seed
    SD&d=D[{kv.first.first,it->second.first}];
    SR&s=kv.second;
    if(s.step>d.maxstep) d.maxstep=s.step;
    if(s.layer>=16) d.endcap=true;
    if(s.nsim>0){ ++d.navail; if(s.won) ++d.nwon; else if(s.step<d.firstloss) d.firstloss=s.step; } }
  for(int ec=0;ec<2;++ec){
    std::vector<double> ms,ef,fl; long n=0,perf=0;
    for(auto&kv:D){SD&d=kv.second; if((int)d.endcap!=ec) continue; if(d.navail<1) continue;
      ++n; ms.push_back(d.maxstep); ef.push_back((double)d.nwon/d.navail);
      if(d.firstloss<99) fl.push_back(d.firstloss); else ++perf; }
    if(n<10) continue;
    printf("  %-9s %-13s %6ld | %5.1f %5.1f | %6.2f %6.2f | %6ld %5.0f%% | %5.1f\n",
      tag, ec?"endcap-touch":"barrel-only", n, qq(ms,50),qq(ms,90),
      qq(ef,50),qq(ef,10), perf,100.*perf/n, fl.empty()?-1:qq(fl,50));
  }
}
void an_seed(){
  printf("\nPER SEED, dead candidates included. A seed stops producing layer-searches when\n");
  printf("its candidate dies, so 'max step' is how deep it actually got. 'eff' = fraction\n");
  printf("of SEARCHED layers where the sim track had a hit and the true hit won.\n");
  printf("force = MC-matched hit forced to win its layer and bypass the chi2<30 cut.\n\n");
  printf("  %-9s %-13s %6s | %11s | %13s | %12s | %5s\n",
         "config","class","seeds","maxstep m/p90","eff med / p10","perfect seeds","1st loss");
  one("val-search-base.root","base");
  one("val-search-force.root","force");
  one("val-search-force16.root","force16");
}
