#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// THE PREDICTION (slide 470 / S16c): the forward-disc loss is the WINDOW, not
// chi2. Scaling disc material should therefore move "outside the window" and
// essentially NOTHING else -- chi2 is already flat in radius. If chi2 moves too,
// the reading is wrong.
void one(const char*tag,const char*fn){
  long v[8]={}, tot=0;                       // verdicts 1..7
  long vs[3][8]={}, ts[3]={};                // by radius: <8, 8-14, >14
  const double Z[12]={25.30,32.30,41.22,52.62,67.17,84.24,110.94,139.70,
                      175.00,200.96,230.77,265.00};
  TFile f(fn); TTree*t=(TTree*)f.Get("miss");
  if(!t){printf("  %-8s no tree\n",tag);return;}
  ValSearchMiss*m=nullptr; t->SetBranchAddress("m",&m);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!m->has_sim_here||m->n_sim_in_layer<=0) continue;
    if(m->layer<16||m->layer>27||m->verdict<1||m->verdict>7) continue;
    ++v[m->verdict]; ++tot;
    const double sh=std::sinh(std::fabs(m->eta)); if(sh<1e-6) continue;
    const double r=Z[m->layer-16]/sh;
    const int b = r<8?0:(r<14?1:2);
    ++vs[b][m->verdict]; ++ts[b]; }
  if(!tot){printf("  %-8s empty\n",tag);return;}
  printf("  %-8s %8ld | %7.2f%% %7.2f%% %7.2f%% %7.2f%% | %8.2f%% || %7.2f%% %7.2f%% %7.2f%%\n",
    tag, tot, 100.0*v[1]/tot, 100.0*v[3]/tot, 100.0*v[5]/tot, 100.0*v[6]/tot,
    100.0*v[7]/tot,
    ts[0]?100.0*vs[0][1]/ts[0]:0.0, ts[1]?100.0*vs[1][1]/ts[1]:0.0,
    ts[2]?100.0*vs[2][1]/ts[2]:0.0);
}
void an_fwdmat(){
  printf("\n===== FORWARD DISC MATERIAL SCALE: WHICH CHANNEL MOVES? =====\n");
  printf("radL scaled in the forward pixel discs ONLY (|z|>22, r<26).\n");
  printf("Denominator: disc layers the sim track really crosses.\n\n");
  printf("  %-8s %8s | %8s %8s %8s %8s | %9s || %-25s\n","scale","n",
         "outside","pre-sel","chi2","outrnk","WON","outside, by r: <8  8-14  >14");
  for(const char*m : {"1.0","2.0","3.5","5.0"})
    one(m, Form("val-fwdmat-%s.root", m));
  printf("\n  PREDICTION: 'outside' falls (most at small r), chi2 stays flat.\n\n");
}
