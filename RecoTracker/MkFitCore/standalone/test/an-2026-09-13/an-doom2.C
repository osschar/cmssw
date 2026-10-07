#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
void an_doom2(const char *fn="val-search-mat1.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("search"); ValSearchHit *h=nullptr; t->SetBranchAddress("h",&h);
  struct S{int step=-1,layer=-1,nsim=-1; bool sc=false,pre=false,kal=false,won=false;
           double c2mc=1e30,c2b=1e30;};
  std::map<std::pair<int,int>,S> SS;
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(h->search_id<0) continue;
    S&s=SS[{h->event,h->search_id}]; s.step=h->step; s.layer=h->layer;
    if(h->n_sim_hits_in_layer>s.nsim) s.nsim=h->n_sim_hits_in_layer;
    if(h->mc_match){ s.sc=true; if(h->passed_preselect) s.pre=true; }
    if(h->had_kalman&&h->chi2>-900){ if(h->chi2<s.c2b) s.c2b=h->chi2;
      if(h->mc_match){ s.kal=true; if(h->chi2<s.c2mc) s.c2mc=h->chi2; } } }
  for(auto&kv:SS) kv.second.won = kv.second.kal && kv.second.c2mc<30. && kv.second.c2mc<=kv.second.c2b;
  printf("\nPER LAYER-SEARCH FUNNEL, RESTRICTED TO SEARCHES WHERE THE SIM TRACK ACTUALLY\n");
  printf("HAS A HIT IN THAT LAYER (countSimHitsInLayer > 0). This removes the inclusive-\n");
  printf("plan ambiguity: a layer the track never crosses is no longer counted as a miss.\n");
  printf("'MC WINS' = true hit reached the Kalman, was under the cut of 30, and had the\n");
  printf("lowest chi2 of anything there -- i.e. it is actually picked up.\n\n");
  printf("  %-13s %5s %9s | %11s %11s | %11s\n",
         "region","step","searches","MC scanned","MC->Kalman","MC WINS");
  for(int ec=0;ec<2;++ec){
    long TN=0,TW=0;
    for(int st=0;st<=9;++st){
      long n=0,a=0,c=0,d=0;
      for(auto&kv:SS){S&s=kv.second; bool e=!(s.layer<16);
        if((int)e!=ec||s.step!=st||s.nsim<=0) continue;
        ++n; if(s.sc)++a; if(s.kal)++c; if(s.won)++d; }
      TN+=n; TW+=d;
      if(n<15) continue;
      printf("  %-13s %5d %9ld | %5ld %5.0f%% %5ld %5.0f%% | %5ld %5.0f%%\n",
        ec?"endcap/disks":"barrel 0-15",st,n,a,100.*a/n,c,100.*c/n,d,100.*d/n); }
    if(TN) printf("  %-13s %5s %9ld | %11s %11s | %5ld %5.0f%%\n",
      ec?"endcap/disks":"barrel 0-15","ALL",TN,"","",TW,100.*TW/TN);
    printf("\n"); }
}
