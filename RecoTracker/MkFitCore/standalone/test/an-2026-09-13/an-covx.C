#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
void an_covx(const char *fn="val-covxport.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("covxport"); ValCovXport *v=nullptr; t->SetBranchAddress("v",&v);
  printf("\n(a) IS C_B RANK 5?  Eigenvalues of its CORRELATION matrix (dimensionless, sum = 6).\n");
  printf("    A rank-5 transport onto a surface must show one eigenvalue at the float32 floor.\n\n");
  printf("  %4s %5s %5s %6s %5s | %s\n","pt","eta","dal","cosinc","scal","eigenvalues, descending");
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if (v->cov_scale!=1.f) continue;
    printf("  %4.0f %5.1f %5.2f %6.2f %5.0f |",v->pt,v->eta,v->dalpha,v->cos_inc,v->cov_scale);
    for(int k=0;k<6;++k) printf(" %9.2e",v->cb_spec[k]);
    printf("   rank %d\n", v->rank_cb); }
  printf("\n(b) DOES THE ENSEMBLE AGREE?  Eigenvalues of C_B^{-1/2} S C_B^{-1/2} inside that rank.\n");
  printf("    Expectation 1.000 for every one; MC error at N=1e4 is sqrt(2/N) = 0.014.\n\n");
  printf("  %4s %5s %5s %6s %5s %6s | %4s %7s %7s | %8s %7s\n",
         "pt","eta","dal","cosinc","scal","nfail","rank","lam_min","lam_max","null_sig","maxpull");
  int nbad=0;
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    double mp=0; for(int k=0;k<6;++k) mp=std::max(mp,(double)std::fabs(v->mean_pull[k]));
    bool bad = (v->rank_cb!=5)||v->lam_min<0.95||v->lam_max>1.05||v->null_sig>0.05||mp>0.05;
    if(bad) ++nbad;
    printf("  %4.0f %5.1f %5.2f %6.2f %5.0f %6d | %4d %7.4f %7.4f | %8.4f %7.3f %s\n",
      v->pt,v->eta,v->dalpha,v->cos_inc,v->cov_scale,v->n_fail,
      v->rank_cb,v->lam_min,v->lam_max,v->null_sig,mp, bad?"  <--":""); }
  printf("\n  %d of %lld configs outside tolerance\n",nbad,t->GetEntries());
}
