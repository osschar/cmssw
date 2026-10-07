#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// The punchline plot: the cut's TRACK term and the residual it must contain,
// both vs |eta|, in cm. They diverge -- that is the whole finding.
void mk350(){
  const double EX=3.0;
  const int NB=10; const double LO=0, HI=2.5;
  std::vector<double> trk[NB], dqR[NB], dqP[NB];
  TFile f("val-miss-pt5c.root"); TTree*t=(TTree*)f.Get("search");
  ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->mc_match||h->layer>=4) continue;
    if(h->sigma_q_trk<0||h->dq<-900) continue;
    double ae=std::fabs(h->eta); if(ae>=HI) continue;
    int b=(int)(ae/(HI-LO)*NB); if(b<0||b>=NB) continue;
    trk[b].push_back(EX*(3*h->sigma_q_trk));                  // the cut's track term = 9 sigma_q
    (h->passed_preselect?dqP:dqR)[b].push_back(std::fabs(h->dq)); }
  auto med=[](std::vector<double>&u){ if(u.empty())return 0.0; std::sort(u.begin(),u.end()); return u[u.size()/2]; };
  TH1F *a=new TH1F("track_term","pixel barrel: the dq cut's track term vs the residual;|#eta|;cm",NB,LO,HI);
  TH1F *b=new TH1F("dq_rejected","",NB,LO,HI);
  TH1F *c=new TH1F("dq_passed","",NB,LO,HI);
  for(int k=0;k<NB;++k){ a->SetBinContent(k+1,med(trk[k]));
    b->SetBinContent(k+1,med(dqR[k])); c->SetBinContent(k+1,med(dqP[k])); }
  a->SetMinimum(0); a->SetMaximum(0.12);
  a->SetDirectory(nullptr); b->SetDirectory(nullptr); c->SetDirectory(nullptr);
  TFile o("/foo/matevz/mic-dev/validation-slides/data/350-q-covariance.root","RECREATE");
  a->Write(); b->Write(); c->Write(); o.Close();
  printf("350: track term %.4f -> %.4f cm ; rejected |dq| %.4f -> %.4f ; passed |dq| %.4f -> %.4f\n",
         med(trk[0]),med(trk[NB-1]),med(dqR[0]),med(dqR[NB-1]),med(dqP[0]),med(dqP[NB-1]));
}
