#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
void an_beta(const char *fn="val-miss-pt5c.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("miss"); ValSearchMiss *m=nullptr; t->SetBranchAddress("m",&m);
  const double eb[]={0,0.4,0.8,1.2,1.6,2.0,2.5,9}; const int NE=7;
  const char* en[NE]={"0.0-0.4","0.4-0.8","0.8-1.2","1.2-1.6","1.6-2.0","2.0-2.5",">2.5"};
  // regions 0 pixel barrel, 1 fwd pixel disks, 2 outer tracker (strips)
  long H[3][NE]={}, W[3][NE]={}, P[3][NE]={}, C[3][NE]={}, O[3][NE]={}, X[3][NE]={};
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(m->sim_label<0||m->verdict<=0) continue;   // verdict 0 excluded by agreement
    int r; if(m->layer<4) r=0;
    else if((m->layer>=16&&m->layer<28)||(m->layer>=38&&m->layer<50)) r=1; else r=2;
    double ae=std::fabs(m->eta); int e=-1;
    for(int b=0;b<NE;++b) if(ae>=eb[b]&&ae<eb[b+1]){e=b;break;}
    if(e<0) continue;
    ++H[r][e];
    switch(m->verdict){ case 7: ++W[r][e]; break; case 3: ++P[r][e]; break;
                        case 5: ++C[r][e]; break; case 6: ++O[r][e]; break;
                        case 1: ++X[r][e]; break; }
  }
  const char* rn[3]={"PIXEL BARREL","FWD PIX DISKS","OUTER TRACKER"};
  printf("\nBEST-HIT EFFICIENCY vs |eta| -- chopped pT5, 100 events.\n");
  printf("Denominator: layer-searches where the sim track REALLY HAS a hit in that layer\n");
  printf("(verdict 0, searching a layer the track never crossed, excluded by agreement).\n");
  printf("WON = the true hit reached the Kalman, passed chi2<30, and had the lowest chi2.\n\n");
  for(int r=0;r<3;++r){
    printf("  %s:\n",rn[r]);
    printf("    %-9s %8s | %7s | %7s %7s %7s %7s\n","|eta|","n","WON %","presel","chi2>30","outrank","outside");
    for(int e=0;e<NE;++e){ if(H[r][e]<40) continue;
      printf("    %-9s %8ld | %6.1f%% | %6.1f%% %6.1f%% %6.1f%% %6.1f%%\n",en[e],H[r][e],
        100.*W[r][e]/H[r][e],100.*P[r][e]/H[r][e],100.*C[r][e]/H[r][e],
        100.*O[r][e]/H[r][e],100.*X[r][e]/H[r][e]); }
    printf("\n"); }
}
