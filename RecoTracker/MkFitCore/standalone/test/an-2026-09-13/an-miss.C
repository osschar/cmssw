#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static double qq(std::vector<double>u,double p){ if(u.empty())return 0; std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
void an_miss(const char *fn="val-search-miss.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("miss"); ValSearchMiss *m=nullptr; t->SetBranchAddress("m",&m);
  const char* vn[6]={"0 no sim hit in layer","1 outside the window","2 in window, not scanned",
                     "3 scanned, failed presel","4 preselected, evicted","5 reached Kalman"};
  long V[2][6]={}; long tot[2]={}; long nolbl=0;
  std::vector<double> D[2], NL[2], PH[2], QQ[2];
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(m->sim_label<0){++nolbl;continue;}
    int r = m->layer<16?0:1;
    if(m->verdict<0||m->verdict>5) continue;
    ++V[r][m->verdict]; ++tot[r];
    if(m->verdict==0 && m->near_d3d>-900) { D[r].push_back(m->near_d3d); NL[r].push_back(m->near_layer); }
    if(m->verdict==1){ PH[r].push_back(m->sim_dphi_norm); QQ[r].push_back(std::fabs(m->sim_dq_norm)); } }
  printf("\nWHY THE SIM TRACK'S OWN HIT WAS NOT CONSIDERED. One row per LAYER-SEARCH,\n");
  printf("%lld searches, %ld with no sim label. Verdicts are mutually exclusive,\n",t->GetEntries(),nolbl);
  printf("ordered most-upstream first. Normalised residuals: |dphi|/phi_delta and the\n");
  printf("position in [q_min,q_max] where 0 = centre and 1 = the edge.\n\n");
  printf("  %-26s %10s %7s | %10s %7s\n","verdict","barrel n","%","endcap n","%");
  for(int v=0;v<6;++v)
    printf("  %-26s %10ld %6.1f%% | %10ld %6.1f%%\n",vn[v],
      V[0][v],tot[0]?100.*V[0][v]/tot[0]:0., V[1][v],tot[1]?100.*V[1][v]/tot[1]:0.);
  printf("  %-26s %10ld         | %10ld\n","TOTAL",tot[0],tot[1]);
  for(int r=0;r<2;++r){ if(D[r].empty()&&PH[r].empty()) continue;
    printf("\n  %s:\n", r?"ENDCAP / DISKS":"BARREL 0-15");
    if(!D[r].empty())
      printf("    verdict 0 -- searched a layer the track never hit. Nearest sim hit is\n"
             "      %.2f cm away (med), p90 %.2f; its layer med %.0f vs searched layer.\n",
             qq(D[r],50),qq(D[r],90),qq(NL[r],50));
    if(!PH[r].empty())
      printf("    verdict 1 -- hit exists but outside the window: |dphi|/phi_delta med %.2f\n"
             "      p90 %.2f ; |dq| in q-window half-widths med %.2f p90 %.2f (n=%zu)\n",
             qq(PH[r],50),qq(PH[r],90),qq(QQ[r],50),qq(QQ[r],90),PH[r].size());
  }
}
