#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
void an_norm(const char *fn="val-covxport.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("covxport"); ValCovXport *v=nullptr; t->SetBranchAddress("v",&v);
  printf("\nNORMAL INCIDENCE ONLY -- where the curvilinear surface IS the target plane,\n");
  printf("so the ensemble and the transported covariance describe the same thing.\n");
  printf("Eigenvalues of C_B^{-1/2} S C_B^{-1/2}, rank 5. Expect 1.000 +- 0.014 (N=1e4).\n");
  printf("cov_scale is the LINEARITY control: correct linear transport is scale-invariant.\n\n");
  printf("  %4s %5s %6s | %8s %8s  %8s %8s | %s\n",
         "pt","eta","dal","min@x1","max@x1","min@x100","max@x100","");
  double s1,s2;
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if (v->cos_inc<0.99f) continue;
    if (v->cov_scale==1.f){ s1=v->lam_min; s2=v->lam_max;
      printf("  %4.0f %5.1f %6.2f | %8.4f %8.4f ",v->pt,v->eta,v->dalpha,s1,s2); }
    else {
      bool b1 = s1<0.95||s2>1.05, b2 = v->lam_min<0.95||v->lam_max>1.05;
      printf(" %8.4f %8.4f | %s\n",v->lam_min,v->lam_max,
             b2&&!b1 ? "nonlinear at x100" : (b1?"  <-- off at x1":"")); }
  }
}
