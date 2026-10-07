#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
void an_inc(const char *fn="val-covxport.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("covxport"); ValCovXport *v=nullptr; t->SetBranchAddress("v",&v);
  printf("\n  %4s %5s %5s | %8s %10s | %14s\n","pt","eta","dal","cos_inc","MEASURED","sum|C_B|");
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if (v->cov_scale!=1.f) continue;
    printf("  %4.0f %5.1f %5.2f | %8.2f %10.4f | %14.7e\n",
      v->pt,v->eta,v->dalpha,v->cos_inc,v->cos_inc_meas,v->covb_sum); }
}
