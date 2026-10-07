#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// The chi2-killed true hits need a 6-10x sigma error to be normal, and we only
// have 1.4-2x. So WHAT are they? Two candidates, and ValSearchMiss can separate
// them without a new run:
//   (a) NOT REAL MATCHES -- mc_match only tests sim_lbl == hit_lbl and the layer,
//       so a hit of the right particle from a far-away part of the track passes.
//       Those show up as large |sim_dq|.
//   (b) REAL matches on a WRECKED candidate -- then |sim_dq| is ordinary and the
//       prediction itself is bad.
static double q(std::vector<double> u,double p){ if(u.empty())return 0;
  std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
void an_chi2who(const char *fn="val-surfq-1.root"){
  struct R{const char*n;int lo,hi;}; const R RG[]={{"PixB 0-3",0,3},{"TBPS 4-9",4,9},{"FwdPix 16-27",16,27}};
  const int NR=3;
  std::vector<double> dqW[NR], dqK[NR], d3W[NR], d3K[NR];
  long nK[NR]={}, nK_far[NR]={};
  TFile f(fn); TTree*t=(TTree*)f.Get("miss");
  ValSearchMiss*m=nullptr; t->SetBranchAddress("m",&m);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!m->has_sim_here||m->n_sim_in_layer<=0||m->sim_dq<-900) continue;
    int r=-1; for(int k=0;k<NR;++k) if(m->layer>=RG[k].lo&&m->layer<=RG[k].hi){r=k;break;}
    if(r<0) continue;
    const double adq=std::fabs(m->sim_dq);
    if(m->verdict==7){ dqW[r].push_back(adq); if(m->sim_d3d>-900) d3W[r].push_back(m->sim_d3d); }
    else if(m->verdict==5){ dqK[r].push_back(adq); if(m->sim_d3d>-900) d3K[r].push_back(m->sim_d3d);
      ++nK[r]; if(adq>1.0) ++nK_far[r]; } }
  printf("\n===== WHAT ARE THE chi2-KILLED HITS? =====\n");
  printf("|sim_dq| = how far the sim track's hit is from the prediction in q, cm.\n");
  printf("A REAL match is ~100 um in pixels, up to the strip length in strips.\n\n");
  printf("  %-14s | %8s %8s %8s | %8s %8s %8s | %9s\n","region",
         "WON med","WON p90","WON n","KILL med","KILL p90","KILL n","|dq|>1cm");
  for(int r=0;r<NR;++r){ if(dqW[r].size()<50||dqK[r].size()<30) continue;
    printf("  %-14s | %8.4f %8.4f %8zu | %8.3f %8.2f %8zu | %8.1f%%\n", RG[r].n,
      q(dqW[r],50), q(dqW[r],90), dqW[r].size(),
      q(dqK[r],50), q(dqK[r],90), dqK[r].size(), 100.0*nK_far[r]/nK[r]); }
  printf("\n  and the 3-D distance to the window centre [cm]:\n");
  printf("  %-14s | %8s %8s | %8s %8s\n","region","WON med","WON p90","KILL med","KILL p90");
  for(int r=0;r<NR;++r){ if(d3W[r].size()<50||d3K[r].size()<30) continue;
    printf("  %-14s | %8.4f %8.4f | %8.3f %8.2f\n", RG[r].n,
      q(d3W[r],50), q(d3W[r],90), q(d3K[r],50), q(d3K[r],90)); }
  printf("\n");
}
