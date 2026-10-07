#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// The forward pixel discs lose ~17% of the layers the track actually crosses.
// Hermeticity is NOT the explanation -- the denominator is layers where the sim
// track HAS a hit, so gaps are already excluded. Two things that could still be
// geometric or kinematic:
//   r  -- edge / ring structure. For a disc at fixed z, r = z/sinh|eta|.
//   pT -- curvature, scattering, and the window size all scale with it.
// verdicts (ValProp.cc:2180): 1 outside 2 not-scanned 3 pre-sel 4 evicted
//                             5 chi2 6 outranked 7 WON
void an_diskwhy(const char *fn="val-purity.root"){
  // disc z, mkFit layers 16..27
  const double Z[12]={25.30,32.30,41.22,52.62,67.17,84.24,110.94,139.70,
                      175.00,200.96,230.77,265.00};
  const int NP=6; const double pb[NP+1]={0,0.7,1.0,1.5,2.5,5.0,1e9};
  const char* pn[NP]={"<0.7","0.7-1","1-1.5","1.5-2.5","2.5-5",">5"};
  const int NRB=6; const double rb[NRB+1]={0,5,8,11,14,17,30};
  const char* rn[NRB]={"<5","5-8","8-11","11-14","14-17",">17"};
  long np[NP][8]={}, nr[NRB][8]={}, tp[NP]={}, tr[NRB]={};
  TFile f(fn); TTree*t=(TTree*)f.Get("miss");
  ValSearchMiss*m=nullptr; t->SetBranchAddress("m",&m);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!m->has_sim_here||m->n_sim_in_layer<=0) continue;
    if(m->layer<16||m->layer>27||m->verdict<1||m->verdict>7) continue;
    if(m->pt<0||m->eta<-900) continue;
    int p=-1; for(int b=0;b<NP;++b) if(m->pt>=pb[b]&&m->pt<pb[b+1]){p=b;break;}
    if(p>=0){ ++np[p][m->verdict]; ++tp[p]; }
    const double sh=std::sinh(std::fabs(m->eta));
    if(sh>1e-6){ const double r=Z[m->layer-16]/sh;
      int q=-1; for(int b=0;b<NRB;++b) if(r>=rb[b]&&r<rb[b+1]){q=b;break;}
      if(q>=0){ ++nr[q][m->verdict]; ++tr[q]; } } }
  auto row=[&](const char*lab,long*v,long tot){
    if(tot<200){ return; }
    printf("  %-9s %8ld | %6.1f%% %6.1f%% %6.1f%% %6.1f%% %6.1f%% | %7.1f%%\n", lab, tot,
      100.0*v[1]/tot,100.0*v[3]/tot,100.0*v[4]/tot,100.0*v[5]/tot,100.0*v[6]/tot,
      100.0*v[7]/tot); };
  printf("\n===== FORWARD PIXEL DISCS: WHY THE LAYER IS LOST =====\n");
  printf("Denominator: disc layers the sim track REALLY crosses (so hermeticity is\n");
  printf("already divided out -- a track through a gap is not in here).\n\n");
  printf("  %-9s %8s | %7s %7s %7s %7s %7s | %8s\n","pT [GeV]","n",
         "outside","pre-sel","evict","chi2","outrnk","WON");
  for(int b=0;b<NP;++b) row(pn[b],np[b],tp[b]);
  printf("\n  %-9s %8s | %7s %7s %7s %7s %7s | %8s\n","r [cm]","n",
         "outside","pre-sel","evict","chi2","outrnk","WON");
  for(int b=0;b<NRB;++b) row(rn[b],nr[b],tr[b]);
  printf("\n  r = z/sinh|eta| at the disc plane. TFPX spans r 3.2-16.1, TEPX 6.5-25.5.\n\n");
}
