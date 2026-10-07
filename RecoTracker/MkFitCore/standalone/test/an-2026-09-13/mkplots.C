#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// Write the histograms the deck reads. One output .root per slide, named by the
// slide, so re-running an analysis updates that slide and nothing else.
static TH1F* mk(const char*n,const char*t,int nb,double lo,double hi){
  TH1F*h=new TH1F(n,t,nb,lo,hi); h->SetDirectory(nullptr); return h; }
static void wr(const char*fn, std::vector<TH1F*> hs){
  TFile f(fn,"RECREATE"); for(auto h:hs) h->Write(); f.Close();
  printf("  wrote %s (%zu histograms)\n",fn,hs.size()); }

// --- 200: backward-fit chi2 per hit, split by seed purity ------------------
void p200(const char*in,const char*out){
  TFile f(in); TTree*t=(TTree*)f.Get("bkfit"); ValFitHit*v=nullptr; t->SetBranchAddress("h",&v);
  const char*nm[3]={"pure","mostly","dirty"};
  std::vector<TH1F*> H;
  for(int c=0;c<3;++c) H.push_back(mk(Form("chi2_%s",nm[c]),
    "per-hit #chi^{2}, backward fit;log_{10} #chi^{2};hits",60,-2.5,3.5));
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(v->step==0||!std::isfinite(v->chi2)||v->chi2<=0||v->good_frac<0) continue;
    int c = v->good_frac>=0.9999f?0:(v->good_frac>=0.8f?1:2);
    H[c]->Fill(std::log10(v->chi2)); }
  wr(out,H);
}
// --- 240 / 250: search chi2 by region, material on/off ---------------------
static int reg(int l){ return l<4?0:(l<10?1:((l>=16&&l<28)||(l>=38&&l<50)?2:3)); }
void p240(const char*on,const char*off,const char*out){
  const char*rn[3]={"pixb","tobps","disks"};
  std::vector<TH1F*> H;
  const char*src[2]={on,off}; const char*sn[2]={"on","off"};
  for(int m=0;m<2;++m) for(int r=0;r<3;++r)
    H.push_back(mk(Form("chi2_%s_%s",rn[r],sn[m]),
      "search #chi^{2}, MC-matched hits;log_{10} #chi^{2};hits",60,-2.5,3.5));
  for(int m=0;m<2;++m){
    TFile f(src[m]); TTree*t=(TTree*)f.Get("search"); ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
      if(!h->mc_match||!h->had_kalman||h->chi2<=0) continue;
      int r=reg(h->layer); if(r>2) continue;
      H[m*3+r]->Fill(std::log10(h->chi2)); } }
  wr(out,H);
}
// --- 230: full-chain chi2 against the analytic chi2_2 ----------------------
void p230(const char*in,const char*out){
  TFile f(in); TTree*t=(TTree*)f.Get("covxport"); ValCovXport*v=nullptr; t->SetBranchAddress("v",&v);
  TH1F*r=mk("ratio","median #chi^{2} / 1.3863, all configs;ratio;configs",40,0.5,1.5);
  TH1F*l=mk("lam","whitened eigenvalues, normal incidence;#lambda;count",40,0.85,1.15);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(v->chi2_med>0) r->Fill(v->chi2_med_ratio);
    if(v->cos_inc>0.99f) for(int k=0;k<v->rank_cb&&k<6;++k) if(v->lam[k]>0) l->Fill(v->lam[k]); }
  wr(out,{r,l});
}
// --- 330 / 340: pixel-barrel pre-selection vs eta --------------------------
void p330(const char*in,const char*out){
  const double EX=3.0,DQF=1.2,DPF=2.0,HPHE=0.0123;
  TH1F*won=mk("won","pixel barrel, best-hit efficiency;|#eta|;fraction",10,0,2.5);
  TH1F*pre=mk("presel","lost to pre-selection;|#eta|;fraction",10,0,2.5);
  TH1F*den=mk("den","",10,0,2.5);
  TH1F*rq=mk("ratio_dq","rejected true hits: residual / its own cut;log_{10} ratio;hits",50,-2,1.5);
  TH1F*rp=mk("ratio_dphi","",50,-2,1.5);
  { TFile f(in); TTree*t=(TTree*)f.Get("miss"); ValSearchMiss*m=nullptr; t->SetBranchAddress("m",&m);
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
      if(m->sim_label<0||m->verdict<=0||m->layer>=4) continue;
      double ae=std::fabs(m->eta); if(ae>=2.5) continue;
      den->Fill(ae); if(m->verdict==7) won->Fill(ae); if(m->verdict==3) pre->Fill(ae); } }
  won->Divide(den); pre->Divide(den);
  { TFile f(in); TTree*t=(TTree*)f.Get("search"); ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
      if(!h->mc_match||h->passed_preselect||h->layer>=4) continue;
      if(h->sigma_q_trk<0||h->sigma_phi_trk<0||h->dq<-900||h->dphi<-900) continue;
      double cq=EX*(3*h->sigma_q_trk)+EX*DQF*h->hit_q_half_len;
      double cp=(3*h->sigma_phi_trk)+DPF*HPHE;
      if(cq>0&&h->dq!=0) rq->Fill(std::log10(std::fabs(h->dq)/cq));
      if(cp>0&&h->dphi!=0) rp->Fill(std::log10(std::fabs(h->dphi)/cp)); } }
  wr(out,{won,pre,rq,rp});
}
