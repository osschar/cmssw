#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static double qq(std::vector<double>u,double p){ if(u.empty())return 0; std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
void an_seq(const char *fn="val-search-mat1.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("search"); ValSearchHit *h=nullptr; t->SetBranchAddress("h",&h);
  // ---- group by layer-search: the unit the best-hit decision is made over
  struct S { double c2_mc=1e30, c2_best=1e30; bool mc_present=false; int step=-1,layer=-1;
             bool endcap=false; double eta=0; };
  std::map<std::pair<int,int>,S> SS;   // (event, search_id) -- search_id is PER EVENT
  int nstep_max=0;
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->had_kalman||h->chi2<-900||h->search_id<0) continue;
    S &s=SS[{h->event,h->search_id}]; s.step=h->step; s.layer=h->layer; s.eta=h->eta;
    s.endcap = !(h->layer<16);
    if(h->chi2<s.c2_best) s.c2_best=h->chi2;
    if(h->mc_match){ s.mc_present=true; if(h->chi2<s.c2_mc) s.c2_mc=h->chi2; }
    if(h->step>nstep_max) nstep_max=h->step; }
  // ---- pickup efficiency per step, split barrel-only vs endcap-touching
  printf("\nPICKUP PER LAYER-SEARCH, grouped by TrLayerSearch (the unit the best-hit\n");
  printf("decision is made over; keyed by (event, search_id) -- search_id is PER EVENT).\n");
  printf("'MC is best' = the best MC-matched hit had the LOWEST chi2 of everything\n");
  printf("that reached the Kalman AND was under the cut of 30 -- i.e. it would actually\n");
  printf("be picked up. chi2 below is the BEST MC-matched hit per search, not the median\n");
  printf("over all of them, so it is lower than the per-hit tables by construction.\n\n");
  printf("  %-14s %5s %8s | %9s %9s | %9s %9s\n",
         "region","step","searches","MC is best","%","MC killed", "by cut %");
  for(int ec=0;ec<2;++ec){
    for(int st=0;st<=std::min(nstep_max,7);++st){
      long n=0,best=0,cut=0;
      for(auto&kv:SS){ S&s=kv.second; if(!s.mc_present) continue;
        if((int)s.endcap!=ec||s.step!=st) continue;
        ++n; if(s.c2_mc>=30.) ++cut; else if(s.c2_mc<=s.c2_best) ++best; }
      if(n<10) continue;
      printf("  %-14s %5d %8ld | %9ld %8.1f%% | %9ld %8.1f%%\n",
        ec?"endcap/disks":"barrel 0-15",st,n,best,100.*best/n,cut,100.*cut/n); }
    printf("\n"); }
  // ---- chi2 vs step
  printf("CHI2 of MC-matched hits vs step (median of chi2_2 = 1.386):\n");
  printf("  %-14s %5s %8s | %8s %8s | %9s\n","region","step","n","chi2 med","p90","sig short");
  for(int ec=0;ec<2;++ec){
    for(int st=0;st<=std::min(nstep_max,7);++st){
      std::vector<double> v;
      for(auto&kv:SS){S&s=kv.second; if(!s.mc_present)continue;
        if((int)s.endcap!=ec||s.step!=st)continue; v.push_back(s.c2_mc);}
      if(v.size()<10) continue; double m=qq(v,50);
      printf("  %-14s %5d %8zu | %8.3f %8.2f | %8.2fx\n",
        ec?"endcap/disks":"barrel 0-15",st,v.size(),m,qq(v,90),std::sqrt(m/1.38629436)); }
    printf("\n"); }
}
