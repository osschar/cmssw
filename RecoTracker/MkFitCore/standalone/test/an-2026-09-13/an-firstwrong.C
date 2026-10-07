#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
#include <map>
#include <vector>
// Does a diverging candidate take its wrong hit because the wrong hit FITS WELL?
// For every impure path, find the FIRST non-MC hit it added and ask (a) what its
// chi2 was, and (b) whether an MC hit was even available in that layer.
// If the first wrong hit is typically a GOOD fit, chi2 cannot prevent this and
// the fix is not a tighter cut.
static double q(std::vector<double> u,double p){ if(u.empty())return 0;
  std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
struct Hit { int step,layer; bool mc,acc; float chi2; };
void an_firstwrong(const char *fn="val-purity.root"){
  std::map<long long,std::vector<Hit>> C;
  TFile f(fn); TTree*t=(TTree*)f.Get("search");
  ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->had_kalman||h->chi2<-900||h->global_seed<0||h->event<0) continue;
    C[(long long)h->event*1000000LL+h->global_seed].push_back(
        {h->step,h->layer,h->mc_match,h->accepted,h->chi2}); }
  std::vector<double> c2, c2_avail, c2_noavail;
  long n=0, avail=0, better=0;
  for(auto&kv:C){
    auto v=kv.second; std::sort(v.begin(),v.end(),
      [](const Hit&a,const Hit&b){return a.step<b.step;});
    int added=0; for(auto&x:v) if(x.acc&&x.chi2<30.f) ++added;
    if(added<3) continue;
    for(size_t i=0;i<v.size();++i){
      if(!(v[i].acc && v[i].chi2<30.f && !v[i].mc)) continue;   // first wrong ADDED hit
      ++n; c2.push_back(v[i].chi2);
      // was an MC hit available in that same layer, and did it fit worse?
      float best_mc = 1e9f; bool has=false;
      for(auto&x:v) if(x.step==v[i].step && x.layer==v[i].layer && x.mc){ has=true; best_mc=std::min(best_mc,x.chi2); }
      if(has){ ++avail; c2_avail.push_back(v[i].chi2); if(v[i].chi2<best_mc) ++better; }
      else c2_noavail.push_back(v[i].chi2);
      break; } }
  printf("\n===== THE FIRST WRONG HIT A DIVERGING CANDIDATE TAKES =====\n");
  printf("%ld candidates with >=3 added hits took a non-MC hit.\n\n", n);
  printf("  %-34s %8s %8s %8s %8s\n","","n","p25","median","p75");
  printf("  %-34s %8zu %8.2f %8.2f %8.2f\n","chi2 of that first wrong hit",
         c2.size(), q(c2,25), q(c2,50), q(c2,75));
  printf("  %-34s %8zu %8.2f %8.2f %8.2f\n","  ...when an MC hit WAS available",
         c2_avail.size(), q(c2_avail,25), q(c2_avail,50), q(c2_avail,75));
  printf("  %-34s %8zu %8.2f %8.2f %8.2f\n","  ...when there was NO MC hit there",
         c2_noavail.size(), q(c2_noavail,25), q(c2_noavail,50), q(c2_noavail,75));
  printf("\n  MC hit available in that layer at all : %ld of %ld = %.1f%%\n",
         avail, n, n?100.0*avail/n:0.0);
  printf("  and the wrong hit fitted BETTER than it: %ld of %ld = %.1f%%\n",
         better, avail, avail?100.0*better/avail:0.0);
  printf("\n  Expected chi2 for a CORRECT hit: median 1.386. If the wrong hit sits\n");
  printf("  near that, chi2 cannot tell them apart and a tighter cut will not help.\n\n");
}
