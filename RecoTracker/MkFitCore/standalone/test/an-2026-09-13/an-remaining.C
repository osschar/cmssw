#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// WHERE DOES THE REMAINING PER-HIT INEFFICIENCY GO, now that pre-selection is
// largely fixed? One row per LAYER-SEARCH in which the sim track really HAS a
// hit -- so the inclusive layer plans do not count as failures.
//
// The verdicts are mutually exclusive and ordered by how early the hit was lost:
//   outside  : the hit exists but fell outside the opened window
//   not-scan : in the window, never visited
//   pre-sel  : scanned, failed the dq/dphi cut
//   evicted  : pre-selected, thrown out of the bounded pqueue
//   chi2     : reached the Kalman, best chi2 over the cut
//   outrank  : reached the Kalman, another hit won the layer
//   WON      : the sim track's hit was the one added
void one(const char*tag, const char*fn){
  struct R{const char*n;int lo,hi;}; 
  const R RG[]={{"PixB 0-3",0,3},{"TBPS 4-9",4,9},{"TOB2S 10-15",10,15},
                {"FwdPix 16-27",16,27},{"TEC 28-37",28,37}};
  const int NR=5, NV=7;
  long c[NR][NV]={}, tot[NR]={};
  TFile f(fn); TTree*t=(TTree*)f.Get("miss");
  if(!t){ printf("  %s: no miss tree\n", tag); return; }
  ValSearchMiss*m=nullptr; t->SetBranchAddress("m",&m);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!m->has_sim_here || m->n_sim_in_layer<=0) continue;   // truth must offer a hit
    int r=-1; for(int k=0;k<NR;++k) if(m->layer>=RG[k].lo&&m->layer<=RG[k].hi){r=k;break;}
    if(r<0) continue;
    // ValProp.cc:2180-2191 -- verdict is assigned there, do NOT re-derive it:
    //   0 no sim hit here (excluded above) | 1 outside window | 2 in window, not
    //   scanned | 3 scanned, failed pre-selection | 4 pre-selected, evicted
    //   5 chi2 (no usable chi2, or >= 30) | 6 outranked | 7 WON
    int v;
    switch(m->verdict){
      case 1: v=0; break;  case 2: v=1; break;  case 3: v=2; break;
      case 4: v=3; break;  case 5: v=4; break;  case 6: v=5; break;
      case 7: v=6; break;  default: continue; }
    ++c[r][v]; ++tot[r]; }
  printf("\n--- %s ---\n", tag);
  printf("  %-14s %8s | %6s %6s %6s %6s %6s %6s | %7s\n","region","n",
         "outsid","notscn","pre-sel","evict","chi2","outrnk","WON");
  for(int r=0;r<NR;++r){ if(tot[r]<200) continue;
    printf("  %-14s %8ld |", RG[r].n, tot[r]);
    for(int v=0;v<6;++v) printf(" %5.1f%%", 100.0*c[r][v]/tot[r]);
    printf(" | %6.1f%%\n", 100.0*c[r][6]/tot[r]); }
}
void an_remaining(){
  printf("\n===== WHERE THE REMAINING PER-HIT LOSS GOES =====\n");
  printf("Denominator: layer-searches where the sim track REALLY HAS a hit in that\n");
  printf("layer. Verdicts mutually exclusive, ordered earliest-loss first.\n");
  one("BASELINE (surface reference OFF)", "val-surfq-0.root");
  one("FIXED (surface reference ON)",     "val-surfq-1.root");
  printf("\n");
}
