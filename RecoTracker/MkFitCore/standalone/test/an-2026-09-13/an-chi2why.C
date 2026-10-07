#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// Are the chi2-killed true hits MARGINAL (consistent with a 1.4-2x covariance
// deficit) or GROSS (so something else)?
//
// chi2 is 2-D: expectation 2, median 1.386. A sigma too small by factor f
// inflates chi2 by f^2. So:
//    f = 1.4  -> a median hit lands at  2.7
//    f = 2.0  -> 5.5
//    f = 4.6  -> 30   <- only here does the CUT start eating median hits
// If the killed ones sit far above 30, the covariance size is not what killed
// them and widening it would not recover them.
static double q(std::vector<double> u,double p){ if(u.empty())return 0;
  std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
void an_chi2why(const char *fn="val-surfq-1.root"){
  struct R{const char*n;int lo,hi;};
  const R RG[]={{"PixB 0-3",0,3},{"TBPS 4-9",4,9},{"FwdPix 16-27",16,27}};
  const int NR=3;
  std::vector<double> won[NR], kill[NR];
  TFile f(fn); TTree*t=(TTree*)f.Get("miss");
  ValSearchMiss*m=nullptr; t->SetBranchAddress("m",&m);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!m->has_sim_here||m->n_sim_in_layer<=0) continue;
    if(m->mc_chi2<-900) continue;
    int r=-1; for(int k=0;k<NR;++k) if(m->layer>=RG[k].lo&&m->layer<=RG[k].hi){r=k;break;}
    if(r<0) continue;
    if(m->verdict==7) won[r].push_back(m->mc_chi2);
    else if(m->verdict==5) kill[r].push_back(m->mc_chi2); }
  printf("\n===== ARE THE chi2-KILLED HITS MARGINAL, OR GROSS? =====\n");
  printf("chi2 is 2-D: expectation 2, median 1.386. sigma too small by f inflates\n");
  printf("chi2 by f^2, so f=1.4 -> 2.7, f=2.0 -> 5.5, f=4.6 -> 30.\n\n");
  printf("  %-14s | %7s %7s %7s | %7s %7s %7s %7s %8s\n","region",
         "n WON","med","p90","n KILL","p10","med","p90","f implied");
  for(int r=0;r<NR;++r){
    if(won[r].size()<50) continue;
    double km = kill[r].size()>=30 ? q(kill[r],50) : -1;
    printf("  %-14s | %7zu %7.2f %7.2f | %7zu %7.1f %7.1f %7.1f %8s\n",
           RG[r].n, won[r].size(), q(won[r],50), q(won[r],90), kill[r].size(),
           kill[r].size()>=30?q(kill[r],10):0.0, km, kill[r].size()>=30?q(kill[r],90):0.0,
           km>0 ? Form("%.1fx", std::sqrt(km/1.386)) : "-");
  }
  printf("\n  'f implied' = sqrt(median_killed / 1.386): the sigma error that would be\n");
  printf("  needed for those hits to be normal. Compare with the MEASURED 1.4-2.0x.\n");
  printf("  WON median vs 1.386 is the honest read of the covariance on hits we keep.\n\n");
}
