#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
void an_chain(const char *fn="val-covxport.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("covxport"); ValCovXport *v=nullptr; t->SetBranchAddress("v",&v);
  printf("\nFULL CHAIN: transport + jacCurv2Loc + local 2-D projection, via the production\n");
  printf("kalmanComputeChi2Plane. Each sample's landing point is a hit with (0.1 um)^2 error.\n");
  printf("If the covariance the UPDATE sees is right, these are chi2 with 2 d.o.f.:\n");
  printf("  mean 2.000   median 1.3863   p90 4.6052   p99 9.2103   (all closed form)\n");
  printf("ratio = median/1.3863. ratio < 1 = over-covered (sigma too wide by 1/sqrt(ratio)).\n\n");
  printf("  %4s %5s %6s %6s %5s | %7s %7s %7s %7s | %6s %8s\n",
         "pt","eta","dal","cosinc","scal","mean","median","p90","p99","ratio","sigma x");
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if (v->chi2_med<=0) continue;
    printf("  %4.0f %5.1f %6.2f %6.2f %5.0f | %7.3f %7.3f %7.3f %7.2f | %6.3f %8.3f %s\n",
      v->pt,v->eta,v->dalpha,v->cos_inc,v->cov_scale,
      v->chi2_mean,v->chi2_med,v->chi2_p90,v->chi2_p99,v->chi2_med_ratio,
      1.0/std::sqrt(v->chi2_med_ratio),
      (v->chi2_med_ratio<0.9||v->chi2_med_ratio>1.1)?"  <--":""); }
}
