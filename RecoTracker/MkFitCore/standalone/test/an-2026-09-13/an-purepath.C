#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
#include <map>
#include <vector>
// PURE PATHS. A candidate's path is the hits actually ADDED: accepted (best in
// layer) AND chi2 < 30. If every added hit is MC-matched the candidate never
// took a wrong hit -- so if such a path STILL shows the ~8x chi2 failures, the
// "acquired = took a wrong hit and followed it" reading is wrong and the state
// goes bad while doing everything right.
static double q(std::vector<double> u,double p){ if(u.empty())return 0;
  std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
struct Hit { int step,layer; bool mc,acc; float chi2,rx,ry; };
void an_purepath(const char *fn="val-purity.root", int ndump=4){
  std::map<long long,std::vector<Hit>> C;
  std::map<long long,float> GF;
  TFile f(fn); TTree*t=(TTree*)f.Get("search");
  ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->had_kalman||h->chi2<-900||h->global_seed<0||h->event<0) continue;
    const long long k=(long long)h->event*1000000LL+h->global_seed;
    C[k].push_back({h->step,h->layer,h->mc_match,h->accepted,h->chi2,
                    h->residual_x*1e4f,h->residual_y*1e4f});
    if(h->seed_good_frac>=0) GF[k]=h->seed_good_frac; }

  long np=0, nim=0;                         // pure / impure paths
  long pk=0,pn=0, ik=0, in_=0;              // mc-hit kill counts on each
  std::vector<double> pres, ires;           // res_x of KILLED mc hits
  std::vector<long long> bad;               // pure paths that still all-kill
  std::vector<long long> bad2;              // IMPURE paths that all-kill -- the divergers
  for(auto&kv:C){
    int added=0, addedmc=0, mcn=0, mck=0;
    for(auto&x:kv.second){
      if(x.acc && x.chi2<30.f){ ++added; if(x.mc) ++addedmc; }
      if(x.mc){ ++mcn; if(x.chi2>=30.f) ++mck; } }
    if(added<3||mcn<2) continue;
    const bool pure = (addedmc==added);
    if(pure){ ++np; pn+=mcn; pk+=mck;
      for(auto&x:kv.second) if(x.mc&&x.chi2>=30.f) pres.push_back(std::fabs(x.rx));
      if(mck==mcn) bad.push_back(kv.first); }
    else    { ++nim; in_+=mcn; ik+=mck;
      for(auto&x:kv.second) if(x.mc&&x.chi2>=30.f) ires.push_back(std::fabs(x.rx));
      if(mck==mcn) bad2.push_back(kv.first); } }

  printf("\n===== PURE PATHS: candidates that never took a wrong hit =====\n");
  printf("path = hits actually ADDED (accepted AND chi2 < 30); >=3 added, >=2 mc hits.\n\n");
  printf("  %-26s %8s %10s %9s | %10s %10s\n","","cands","mc hits","kill rate",
         "killed res_x md","p90");
  printf("  %-26s %8ld %10ld %8.2f%% | %10.1f %10.1f\n","PURE path (all added mc)",
         np, pn, pn?100.0*pk/pn:0.0, q(pres,50), q(pres,90));
  printf("  %-26s %8ld %10ld %8.2f%% | %10.1f %10.1f\n","impure path", nim, in_,
         in_?100.0*ik/in_:0.0, q(ires,50), q(ires,90));
  printf("\n  %ld PURE paths still lost EVERY mc hit to chi2.\n", (long)bad.size());

  printf("  %ld IMPURE paths lost every mc hit -- these are the divergers.\n", (long)bad2.size());
  const std::vector<long long> &dump = bad.empty() ? bad2 : bad;
  if(!dump.empty() && ndump>0){
    printf("\n  --- hit by hit, %d of the %s ---\n", ndump,
           bad.empty() ? "IMPURE all-killed (divergers)" : "pure all-killed");
    int d=0;
    for(long long key : dump){ if(d++>=ndump) break;
      auto v=C[key]; std::sort(v.begin(),v.end(),
        [](const Hit&a,const Hit&b){return a.step<b.step;});
      printf("  event %lld seed %lld   seed purity %.2f\n",
             key/1000000LL, key%1000000LL, GF.count(key)?GF[key]:-1.0);
      printf("     %5s %6s %4s %4s %10s %10s %10s\n","step","layer","mc","add","chi2","res_x um","res_y um");
      for(auto&x:v)
        printf("     %5d %6d %4s %4s %10.2f %10.1f %10.1f\n", x.step, x.layer,
               x.mc?"yes":"-", (x.acc&&x.chi2<30.f)?"YES":"-", x.chi2, x.rx, x.ry);
    }
  }
  printf("\n");
}
