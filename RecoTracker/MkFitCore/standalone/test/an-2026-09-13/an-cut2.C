#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static double qq(std::vector<double>u,double p){ if(u.empty())return 0; std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
void an_cut2(const char *fn="val-search-mat1.root"){
  // The production cut, MkFinderV2p2.cc:676 --
  //   ddq   <  EXTRA_DQ*dq_track + EXTRA_DQ*DDQ_PRESEL_FAC*hit_q_half_len
  //   ddphi <  dphi_track        + DDPHI_PRESEL_FAC*HIT_PHI_HALF_EXTENT
  // dq_track / dphi_track are 3 sigma; ValSearchHit stores sigma = track/3.
  const double EXTRA_DQ=3.0, DDQ_F=1.2, DDPHI_F=2.0, HPHE=0.0123;
  TFile f(fn); TTree *t=(TTree*)f.Get("search"); ValSearchHit *h=nullptr; t->SetBranchAddress("h",&h);
  const double eb[]={0,0.8,1.2,1.6,2.0,9}; const int NE=5;
  const char* en[NE]={"< 0.8","0.8-1.2","1.2-1.6","1.6-2.0","> 2.0"};
  long F[NE]={},Dq[NE]={},Dp[NE]={},Both[NE]={};
  std::vector<double> RQ[NE], RP[NE];
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->mc_match||h->passed_preselect) continue;       // failed pre-selection
    if(h->layer>=4) continue;                             // PIXEL BARREL only
    if(h->sigma_q_trk<0||h->sigma_phi_trk<0||h->dq<-900||h->dphi<-900) continue;
    double ae=std::fabs(h->eta); int e=-1;
    for(int b=0;b<NE;++b) if(ae>=eb[b]&&ae<eb[b+1]){e=b;break;}
    if(e<0) continue;
    double rq = EXTRA_DQ*(3*h->sigma_q_trk) + EXTRA_DQ*DDQ_F*h->hit_q_half_len;
    double rp = (3*h->sigma_phi_trk) + DDPHI_F*HPHE;
    bool fq = std::fabs(h->dq) >= rq, fp = std::fabs(h->dphi) >= rp;
    ++F[e]; if(fq&&fp) ++Both[e]; else if(fq) ++Dq[e]; else if(fp) ++Dp[e];
    if(rq>0) RQ[e].push_back(std::fabs(h->dq)/rq);
    if(rp>0) RP[e].push_back(std::fabs(h->dphi)/rp); }
  printf("\nWHICH HALF OF THE PRE-SELECTION CUT REJECTS THE TRUE PIXEL-BARREL HIT?\n");
  printf("Cut (MkFinderV2p2.cc:676), EXTRA_DQ=3, DDQ_PRESEL_FAC=1.2, DDPHI_PRESEL_FAC=2,\n");
  printf("HIT_PHI_HALF_EXTENT=0.0123:\n");
  printf("   ddq   <  3*dq_track + 3*1.2*hit_q_half_len\n");
  printf("   ddphi <    dphi_track +   2*0.0123          <-- note: NO EXTRA factor\n");
  printf("Ratios below are |residual| / (its own cut), so > 1 means that half rejects.\n\n");
  printf("  %-9s %8s | %8s %8s %8s | %9s %9s | %9s %9s\n",
         "|eta|","rejected","dq only","dphi only","both","dq ratio","p90","dphi ratio","p90");
  for(int e=0;e<NE;++e){ if(F[e]<20) continue;
    printf("  %-9s %8ld | %7.1f%% %7.1f%% %7.1f%% | %9.2f %9.2f | %9.2f %9.2f\n",en[e],F[e],
      100.*Dq[e]/F[e],100.*Dp[e]/F[e],100.*Both[e]/F[e],
      qq(RQ[e],50),qq(RQ[e],90),qq(RP[e],50),qq(RP[e],90)); }
}
