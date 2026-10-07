#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// The chi2-killed true hits need a 6-10x sigma error to be normal and we only
// have 1.4-2x, so the covariance SIZE is not what kills them. Decompose their
// residual in the MODULE FRAME, which separates three different failures:
//
//   residual_z  off the module plane -- MUST be ~0. The hit lies on its own
//               module's plane by construction and the prediction is solved onto
//               that plane, so anything here is a propagation/solver failure.
//   residual_x  across strip / the PRECISE (phi) direction -- a genuinely wrong
//               prediction in phi.
//   residual_y  along strip / the COARSE direction -- in TBPS this is the tilted
//               strip projection, where the material is known over-counted ~3x.
static double q(std::vector<double> u,double p){ if(u.empty())return 0;
  std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
void an_killres(const char *fn="val-surfq-1.root"){
  struct R{const char*n;int lo,hi;};
  const R RG[]={{"PixB 0-3",0,3},{"TBPS 4-9",4,9},{"TOB2S 10-15",10,15},{"FwdPix 16-27",16,27}};
  const int NR=4;
  std::vector<double> X[NR][2], Y[NR][2], Z[NR][2];   // [region][0=kept 1=killed]
  TFile f(fn); TTree*t=(TTree*)f.Get("search");
  ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->mc_match||!h->had_kalman||h->chi2<-900) continue;
    if(h->residual_x<-900||h->residual_y<-900||h->residual_z<-900) continue;
    int r=-1; for(int k=0;k<NR;++k) if(h->layer>=RG[k].lo&&h->layer<=RG[k].hi){r=k;break;}
    if(r<0) continue;
    const int c = (h->chi2>=30.f) ? 1 : 0;
    X[r][c].push_back(std::fabs(h->residual_x)*1e4);   // um
    Y[r][c].push_back(std::fabs(h->residual_y)*1e4);
    Z[r][c].push_back(std::fabs(h->residual_z)*1e4); }
  printf("\n===== MODULE-FRAME RESIDUALS OF chi2-KILLED TRUE HITS =====\n");
  printf("All in um. x = across strip (PRECISE/phi), y = along strip (COARSE),\n");
  printf("z = OFF THE MODULE PLANE and must be ~0 -- anything there is a\n");
  printf("propagation or solver failure, not a covariance problem.\n\n");
  const char* cn[2]={"kept  ","KILLED"};
  printf("  %-13s %-7s %7s | %9s %9s | %9s %9s | %9s %9s\n","region","",
         "n","|res_x| md","p90","|res_y| md","p90","|res_z| md","p90");
  for(int r=0;r<NR;++r){
    if(X[r][0].size()<50) continue;
    for(int c=0;c<2;++c){
      if(X[r][c].size()<20) continue;
      printf("  %-13s %-7s %7zu | %9.1f %9.1f | %9.1f %9.1f | %9.2f %9.2f\n",
             c?"":RG[r].n, cn[c], X[r][c].size(),
             q(X[r][c],50), q(X[r][c],90), q(Y[r][c],50), q(Y[r][c],90),
             q(Z[r][c],50), q(Z[r][c],90)); } }
  printf("\n  Ratios KILLED/kept, which say WHICH direction failed:\n");
  printf("  %-13s | %9s %9s %9s\n","region","x (phi)","y (coarse)","z (plane)");
  for(int r=0;r<NR;++r){
    if(X[r][0].size()<50||X[r][1].size()<20) continue;
    printf("  %-13s | %9.1f %9.1f %9.1f\n", RG[r].n,
           q(X[r][1],50)/std::max(1e-9,q(X[r][0],50)),
           q(Y[r][1],50)/std::max(1e-9,q(Y[r][0],50)),
           q(Z[r][1],50)/std::max(1e-9,q(Z[r][0],50))); }
  printf("\n");
}
