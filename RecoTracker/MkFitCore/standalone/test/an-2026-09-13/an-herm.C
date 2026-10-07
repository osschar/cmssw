#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static double qq(std::vector<double>u,double p){ if(u.empty())return 0; std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
void an_herm(const char *fn="val-search.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("search"); ValSearchHit *h=nullptr; t->SetBranchAddress("h",&h);
  struct A{ std::vector<double> d,rs,da,ad; long n=0,wc=0; };
  std::map<int,A> G; A all; long nno=0;
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->passed_preselect) continue;
    if(h->d_hermite_exact<0){ ++nno; continue; }
    int g = h->layer<4?0:(h->layer<10?1:3);
    A&a=G[g]; ++a.n; if(h->wrong_crossing) ++a.wc;
    a.d.push_back(h->d_hermite_exact); a.ad.push_back(std::fabs(h->dalpha_exact));
    if(h->rel_ds>-900) a.rs.push_back(std::fabs(h->rel_ds));
    all.d.push_back(h->d_hermite_exact); if(h->rel_ds>-900) all.rs.push_back(std::fabs(h->rel_ds));
    all.ad.push_back(std::fabs(h->dalpha_exact)); ++all.n; if(h->wrong_crossing) ++all.wc; }
  printf("\nIS THE HERMITE STEP AS GOOD AS AN EXACT SOLVE?\n");
  printf("Reference: exact uniform-B helix from the candidate's own state, solved against\n");
  printf("the module plane in DOUBLE. Uniform B is right here -- the mini-propagator uses\n");
  printf("Config::Bfield (default arg, MkFinderV2p2.cc:413), not the parametric field.\n");
  printf("Scale: float32 ULP at r 30-60 cm is 4-8e-6 cm; a pixel pitch is 25-100 um.\n");
  printf("%ld pre-selected hits had no usable reference.\n\n",nno);
  const char* gn[4]={"PixB 0-3","TOB PS 4-9","","fwd pix / disks"};
  printf("  %-16s %7s | %9s %9s %9s | %9s %9s | %8s %7s\n",
         "region","n","|dpos| med","p90","max","|ds|/s med","p90","med|da|","wrong%");
  for(auto&kv:G){A&a=kv.second; if(a.n<30)continue;
    printf("  %-16s %7ld | %9.2e %9.2e %9.2e | %9.2e %9.2e | %8.4f %7.2f\n",
      gn[kv.first],a.n,qq(a.d,50),qq(a.d,90),qq(a.d,100),
      qq(a.rs,50),qq(a.rs,90),qq(a.ad,50),100.*a.wc/a.n); }
  printf("  %-16s %7ld | %9.2e %9.2e %9.2e | %9.2e %9.2e | %8.4f %7.2f\n",
    "ALL",all.n,qq(all.d,50),qq(all.d,90),qq(all.d,100),
    qq(all.rs,50),qq(all.rs,90),qq(all.ad,50),100.*all.wc/all.n);
}
