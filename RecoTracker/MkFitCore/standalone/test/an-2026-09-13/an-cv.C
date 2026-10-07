#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
void an_cv(const char *fn="val-covxport.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("covxport"); ValCovXport *v=nullptr; t->SetBranchAddress("v",&v);
  printf("\nENSEMBLE vs TRANSPORTED COVARIANCE. Eigenvalues of C_B^{-1/2} S C_B^{-1/2},\n");
  printf("rank 5 in every config. Expectation 1.000; MC error at N=1e4 is sqrt(2/N)=0.014.\n");
  printf("  'on plane'       : samples as the propagator left them -- fair ONLY at cos_inc=1\n");
  printf("  'curvilinear'    : samples slid along their own momentum onto the surface the\n");
  printf("                     transported covariance actually lives on -- fair everywhere\n\n");
  printf("  %4s %5s %5s %6s %5s | %8s %8s | %8s %8s\n",
         "pt","eta","dal","cosinc","scal","onpl_min","onpl_max","curv_min","curv_max");
  int nbad=0;
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    bool bad = v->lam_min_cv<0.95 || v->lam_max_cv>1.05;
    if(bad) ++nbad;
    printf("  %4.0f %5.1f %5.2f %6.2f %5.0f | %8.4f %8.4f | %8.4f %8.4f %s\n",
      v->pt,v->eta,v->dalpha,v->cos_inc,v->cov_scale,
      v->lam_min,v->lam_max,v->lam_min_cv,v->lam_max_cv, bad?"  <--":""); }
  printf("\n  %d of %lld configs outside [0.95, 1.05] on the CURVILINEAR surface\n",nbad,t->GetEntries());
}
