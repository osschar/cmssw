#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static double qq(std::vector<double>u,double p){ if(u.empty())return 0; std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
struct A{ std::vector<double> c2,pull; long npre=0,nsc=0; };
static void one(const char*fn,const char*tag){
  TFile f(fn); TTree *t=(TTree*)f.Get("search"); ValSearchHit *h=nullptr; t->SetBranchAddress("h",&h);
  std::map<int,A> G; A all;
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->mc_match) continue;
    int g=h->layer<4?0:(h->layer<10?1:3); A&a=G[g]; ++a.nsc; ++all.nsc;
    if(h->passed_preselect){++a.npre;++all.npre;}
    if(h->had_kalman&&h->chi2>-900){ a.c2.push_back(h->chi2); all.c2.push_back(h->chi2); }
    if(h->sigma_q_trk>0&&h->hit_q_half_len>-900&&h->dq>-900){
      double st=std::sqrt(h->sigma_q_trk*h->sigma_q_trk+h->hit_q_half_len*h->hit_q_half_len/3.);
      double p=std::fabs(h->dq)/st; if(p<5){a.pull.push_back(p); all.pull.push_back(p);} } }
  const char* gn[4]={"PixB 0-3","TOB PS 4-9","","fwd disks"};
  for(auto&kv:G){A&a=kv.second; if(a.nsc<50)continue;
    printf("  %-10s %-14s %7ld %7ld | %8.3f %8.2f | %8.3f | %7.2fx\n",tag,gn[kv.first],
      a.nsc,a.npre,qq(a.c2,50),qq(a.c2,90),qq(a.pull,50),std::sqrt(qq(a.c2,50)/1.38629436)); }
  printf("  %-10s %-14s %7ld %7ld | %8.3f %8.2f | %8.3f | %7.2fx\n",tag,"ALL",
    all.nsc,all.npre,qq(all.c2,50),qq(all.c2,90),qq(all.pull,50),std::sqrt(qq(all.c2,50)/1.38629436));
  printf("\n");
}
void an_mat(){
  printf("\nDOES THE SEARCH'S MATERIAL DO ANYTHING? MC-matched hits, 100 events.\n");
  printf("chi2 is 2-D in the module plane: expectation 2, median of chi2_2 = 1.386.\n");
  printf("q pull clamped to < 5 sigma_tot; expect 0.67-0.87.\n\n");
  printf("  %-10s %-14s %7s %7s | %8s %8s | %8s | %8s\n",
    "material","region","scanned","presel","chi2 med","p90","q pull","sig short");
  one("val-search-mat1.root","ON");
  one("val-search-mat0.root","OFF");
}
