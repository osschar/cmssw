#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static double qq(std::vector<double>u,double p){ if(u.empty())return 0; std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
void an_dqterm(const char *fn="val-miss-pt5c.root"){
  const double EX=3.0, DQF=1.2;
  TFile f(fn); TTree*t=(TTree*)f.Get("search"); ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
  const double eb[]={0,0.8,1.6,2.0,9}; const int NE=4;
  const char* en[NE]={"< 0.8","0.8-1.6","1.6-2.0","> 2.0"};
  // [0] rejected by pre-selection, [1] accepted -- compare the two
  std::vector<double> DQ[NE][2], ST[NE][2], HL[NE][2], CUT[NE][2];
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->mc_match||h->layer>=4) continue;
    if(h->sigma_q_trk<0||h->dq<-900||h->hit_q_half_len<-900) continue;
    double ae=std::fabs(h->eta); int e=-1;
    for(int b=0;b<NE;++b) if(ae>=eb[b]&&ae<eb[b+1]){e=b;break;}
    if(e<0) continue;
    int k = h->passed_preselect?1:0;
    DQ[e][k].push_back(std::fabs(h->dq));
    ST[e][k].push_back(EX*(3*h->sigma_q_trk));           // the TRACK term of the cut
    HL[e][k].push_back(EX*DQF*h->hit_q_half_len);        // the HIT term of the cut
    CUT[e][k].push_back(EX*(3*h->sigma_q_trk)+EX*DQF*h->hit_q_half_len); }
  printf("\nWHAT IS THE dq CUT MADE OF, AND HOW BIG IS THE RESIDUAL IT REJECTS?\n");
  printf("Pixel barrel, MC-matched hits. All in cm. Cut = track term + hit term, where\n");
  printf("  track term = 3 * dq_track  (dq_track is itself 3 sigma)\n");
  printf("  hit term   = 3 * 1.2 * hit_q_half_length\n");
  printf("A pixel is ~150 um long in z, so the hit term is intrinsically small.\n\n");
  printf("  %-9s %-9s %7s | %9s %9s | %9s %9s %9s\n",
         "|eta|","presel","n","|dq| med","p90","cut med","track term","hit term");
  for(int e=0;e<NE;++e) for(int k=1;k>=0;--k){ if(DQ[e][k].size()<30) continue;
    printf("  %-9s %-9s %7zu | %9.4f %9.4f | %9.4f %9.4f %9.4f\n",
      k?en[e]:"", k?"PASSED":"rejected", DQ[e][k].size(),
      qq(DQ[e][k],50), qq(DQ[e][k],90), qq(CUT[e][k],50), qq(ST[e][k],50), qq(HL[e][k],50)); }
  printf("\n  If the rejected |dq| is much larger than the PASSED |dq|, the prediction is\n");
  printf("  genuinely off for those tracks and widening the cut only lets in noise.\n");
  printf("  If they overlap, the cut is simply too narrow and widening it is the fix.\n");
}
