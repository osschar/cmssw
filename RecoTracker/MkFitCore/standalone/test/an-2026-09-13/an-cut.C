#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static double qq(std::vector<double>u,double p){ if(u.empty())return 0; std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
void an_cut(const char *fn="val-search.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("search"); ValSearchHit *h=nullptr; t->SetBranchAddress("h",&h);
  struct A{ std::vector<double> c2; long n=0,g5=0,g30=0; };
  std::map<int,A> R[2];
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->had_kalman||h->chi2<-900) continue;
    int g = h->layer<4?0 : (h->layer<10?1 : 3);
    A&a=R[h->mc_match?1:0][g]; ++a.n; a.c2.push_back(h->chi2);
    if(h->chi2>5.f) ++a.g5; if(h->chi2>30.f) ++a.g30; }
  const char* gn[4]={"PixB 0-3","TOB PS 4-9","","fwd pix / disks"};
  printf("\nCHI2 OF HITS REACHING THE KALMAN, against the production cut.\n");
  printf("chi2 is 2-D in the module plane: expectation 2, median of chi2_2 = 1.386.\n");
  printf("The cut is a hardcoded bChi2 < 30 applied to the BEST hit in the layer\n");
  printf("(MkFinderV2p2.cc:938); >5 is where a held-back 'missed' candidate is added.\n");
  printf("So the '> 30' column is an UPPER BOUND on loss -- a good hit only actually\n");
  printf("dies if it was also the best one in its layer.\n\n");
  for(int m=1;m>=0;--m){
    printf("  %s:\n", m?"MC-MATCHED":"other");
    printf("    %-16s %7s | %8s %8s %8s | %8s %8s | %9s\n",
           "region","n","median","p90","p99","> 5","> 30","sigma too small by");
    for(auto&kv:R[m]){ A&a=kv.second; if(a.n<30) continue;
      double med=qq(a.c2,50);
      printf("    %-16s %7ld | %8.3f %8.1f %8.1f | %7.1f%% %7.1f%% | %9.2fx\n",
        gn[kv.first],a.n,med,qq(a.c2,90),qq(a.c2,99),
        100.*a.g5/a.n,100.*a.g30/a.n, std::sqrt(med/1.38629436)); }
    printf("\n"); }
}
