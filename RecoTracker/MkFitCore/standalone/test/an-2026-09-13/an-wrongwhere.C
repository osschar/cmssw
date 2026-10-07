#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
#include <map>
#include <vector>
// WHERE does a diverging candidate take its first wrong hit, and is chi2 able to
// tell? Split by region. NOTE the coverage: this is the INWARD (chopped-pT5)
// search only, which scans PixB and FwdPix heavily and the strips thinly, so the
// strip rows are indicative at best.
static double q(std::vector<double> u,double p){ if(u.empty())return 0;
  std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
struct Hit { int step,layer; bool mc,acc; float chi2; };
void an_wrongwhere(const char *fn="val-purity.root"){
  struct R{const char*n;int lo,hi;};
  const R RG[]={{"PixB 0-3",0,3},{"TBPS 4-9",4,9},{"TOB2S 10-15",10,15},
                {"FwdPix 16-27",16,27},{"TEC 28-37",28,37}};
  const int NR=5;
  std::vector<double> c2[NR];
  long n[NR]={}, avail[NR]={}, better[NR]={}, nohit[NR]={};
  std::map<long long,std::vector<Hit>> C;
  TFile f(fn); TTree*t=(TTree*)f.Get("search");
  ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->had_kalman||h->chi2<-900||h->global_seed<0||h->event<0) continue;
    C[(long long)h->event*1000000LL+h->global_seed].push_back(
      {h->step,h->layer,h->mc_match,h->accepted,h->chi2}); }
  for(auto&kv:C){
    auto v=kv.second; std::sort(v.begin(),v.end(),
      [](const Hit&a,const Hit&b){return a.step<b.step;});
    int added=0; for(auto&x:v) if(x.acc&&x.chi2<30.f) ++added;
    if(added<3) continue;
    for(size_t i=0;i<v.size();++i){
      if(!(v[i].acc && v[i].chi2<30.f && !v[i].mc)) continue;
      int r=-1; for(int k=0;k<NR;++k) if(v[i].layer>=RG[k].lo&&v[i].layer<=RG[k].hi){r=k;break;}
      if(r<0) break;
      ++n[r]; c2[r].push_back(v[i].chi2);
      float bmc=1e9f; bool has=false;
      for(auto&x:v) if(x.step==v[i].step&&x.layer==v[i].layer&&x.mc){has=true;bmc=std::min(bmc,x.chi2);}
      if(has){ ++avail[r]; if(v[i].chi2<bmc) ++better[r]; } else ++nohit[r];
      break; } }
  printf("\n===== WHERE THE FIRST WRONG HIT IS TAKEN =====\n");
  printf("INWARD (chopped-pT5) search only -- see the coverage caveat.\n");
  printf("A CORRECT hit has chi2 median 1.386.\n\n");
  printf("  %-14s %8s %7s | %8s %9s | %9s %10s\n","region","n","share",
         "chi2 med","chi2 p25","no MC hit","wrong won");
  long tot=0; for(int r=0;r<NR;++r) tot+=n[r];
  for(int r=0;r<NR;++r){ if(n[r]<50) continue;
    printf("  %-14s %8ld %6.1f%% | %8.2f %9.2f | %8.1f%% %9.1f%%\n",
      RG[r].n, n[r], 100.0*n[r]/tot, q(c2[r],50), q(c2[r],25),
      100.0*nohit[r]/n[r], avail[r]?100.0*better[r]/avail[r]:0.0); }
  printf("\n  'no MC hit'  = the layer had none, so a wrong hit beat a HOLE.\n");
  printf("  'wrong won'  = the true hit WAS there and the wrong one fitted better.\n\n");
}
