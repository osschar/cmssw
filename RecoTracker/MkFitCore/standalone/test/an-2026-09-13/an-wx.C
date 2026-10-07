#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
void an_wx(const char *fn="val-search.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("search"); ValSearchHit *h=nullptr; t->SetBranchAddress("h",&h);
  long nall=0, nwx_all=0, npre=0, nwx_pre=0;
  printf("\nWRONG-CROSSING CASES: the plane root nearest the Hermite's own turn angle is NOT\n");
  printf("the root of smallest |alpha|, i.e. the cubic converged to a LATER crossing of the\n");
  printf("module's infinite plane. alpha_in/out are the layer's bounding-surface crossings\n");
  printf("solved exactly from the same state -- what sp1/sp2 should have been.\n\n");
  printf("  %3s %5s %5s %5s %5s %6s | %7s %6s | %8s %8s %8s %8s | %7s %7s %5s %s\n",
    "ev","seed","gsd","sim","lay","pt","eta","nroot","da_herm","da_exact","a_small","t_herm",
    "a_in","a_out","pre","|dpos| cm");
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(h->d_hermite_exact<0) continue;
    ++nall; if(h->wrong_crossing) ++nwx_all;
    if(h->passed_preselect){ ++npre; if(h->wrong_crossing) ++nwx_pre; }
    if(!h->wrong_crossing) continue;
    printf("  %3d %5d %5d %5d %5d %6.2f | %7.3f %6d | %8.4f %8.4f %8.4f %8.3f | %7.4f %7.4f %5s %.3e\n",
      h->event,h->seed,h->global_seed,h->sim,h->layer,h->pt,h->eta,h->n_roots,
      h->dalpha_hermite,h->dalpha_exact,h->alpha_small,h->t_hermite,
      h->alpha_in,h->alpha_out, h->passed_preselect?"yes":"no", h->d_hermite_exact); }
  printf("\n  %ld of %ld hits with a reference are wrong-crossing (%.3f%%);\n",nwx_all,nall,100.*nwx_all/std::max(1L,nall));
  printf("  %ld of %ld PRE-SELECTED ones (%.3f%%).\n",nwx_pre,npre,100.*nwx_pre/std::max(1L,npre));
}
