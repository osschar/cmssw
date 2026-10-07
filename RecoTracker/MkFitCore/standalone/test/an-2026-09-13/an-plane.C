#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
void an_plane(const char *fn="val-covxport.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("covxport"); ValCovXport *v=nullptr; t->SetBranchAddress("v",&v);
  printf("\nDIRECT CHECK: |n.(x - p)|, distance of the PROPAGATED state from the plane it\n");
  printf("was propagated TO. No eigen-decomposition. A successful propagation gives ~0.\n");
  printf("Scale: float32 ULP at r ~ 30-60 cm is 4e-6 cm; a sensor is 300 um = 0.03 cm.\n\n");
  printf("  %4s %5s %5s %6s %5s | %10s | %10s %10s %10s | %8s\n",
         "pt","eta","dal","cosinc","scal","ref [cm]","med [cm]","p90","max","null_sig");
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if (v->cov_scale!=1.f) continue;
    printf("  %4.0f %5.1f %5.2f %6.2f %5.0f | %10.3e | %10.3e %10.3e %10.3e | %8.3f\n",
      v->pt,v->eta,v->dalpha,v->cos_inc,v->cov_scale,
      v->d_plane_ref,v->d_plane_med,v->d_plane_p90,v->d_plane_max,v->null_sig); }
}
