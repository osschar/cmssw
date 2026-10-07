#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static double qq(std::vector<double>u,double p){ if(u.empty())return 0; std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
void an_pix(const char *fn="val-miss-pt5c.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("miss"); ValSearchMiss *m=nullptr; t->SetBranchAddress("m",&m);
  const char* vn[8]={"0 no sim hit in layer","1 OUTSIDE the window","2 in window, not scanned",
                     "3 scanned, failed presel","4 preselected, evicted",
                     "5 Kalman, chi2 >= 30","6 Kalman, outranked","7 Kalman, WON"};
  // regions: 0 = pixel barrel 0-3, 1 = forward pixel disks 16-27 & 38-49
  long V[2][8]={}; long tot[2]={};
  std::vector<double> PH[2],QN[2],C5[2],C7[2],D3[2];
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(m->sim_label<0||m->verdict<0) continue;
    int r;
    if (m->layer<4) r=0;
    else if ((m->layer>=16&&m->layer<28)||(m->layer>=38&&m->layer<50)) r=1;
    else continue;
    ++V[r][m->verdict]; ++tot[r];
    if(m->verdict==1){ PH[r].push_back(m->sim_dphi_norm); QN[r].push_back(std::fabs(m->sim_dq_norm));
                       if(m->sim_d3d>-900) D3[r].push_back(m->sim_d3d); }
    if(m->verdict==5&&m->mc_chi2>-900) C5[r].push_back(m->mc_chi2);
    if(m->verdict==7&&m->mc_chi2>-900) C7[r].push_back(m->mc_chi2); }
  printf("\nCHOPPED pT5, PIXEL LAYERS ONLY: what happens to the sim track's own pixel hit.\n");
  printf("One row per layer-search. chi2 is 2-D in the module plane, median of chi2_2 = 1.386.\n\n");
  printf("  %-26s %10s %7s | %10s %7s\n","verdict","pix barrel","%","fwd pix disks","%");
  for(int v=0;v<8;++v)
    printf("  %-26s %10ld %6.1f%% | %10ld %6.1f%%\n",vn[v],
      V[0][v],tot[0]?100.*V[0][v]/tot[0]:0., V[1][v],tot[1]?100.*V[1][v]/tot[1]:0.);
  printf("  %-26s %10ld         | %10ld\n","TOTAL",tot[0],tot[1]);
  for(int r=0;r<2;++r){
    long has = tot[r]-V[r][0];
    if(has<10) continue;
    printf("\n  %s -- restricted to the %ld searches where the hit IS in that layer:\n",
      r?"FWD PIX DISKS":"PIXEL BARREL", has);
    for(int v=1;v<8;++v) printf("      %-26s %6.1f%%\n",vn[v],100.*V[r][v]/has);
    if(!PH[r].empty())
      printf("    outside the window by: |dphi|/phi_delta med %.2f p90 %.2f ; |dq| in q-half-widths\n"
             "      med %.2f p90 %.2f ; 3-D distance to window centre med %.3f cm p90 %.3f\n",
        qq(PH[r],50),qq(PH[r],90),qq(QN[r],50),qq(QN[r],90),qq(D3[r],50),qq(D3[r],90));
    if(!C5[r].empty()) printf("    killed by chi2>=30: chi2 med %.1f p90 %.1f (n=%zu)\n",
        qq(C5[r],50),qq(C5[r],90),C5[r].size());
    if(!C7[r].empty()) printf("    WON:                chi2 med %.2f p90 %.2f (n=%zu)\n",
        qq(C7[r],50),qq(C7[r],90),C7[r].size());
  }
}
