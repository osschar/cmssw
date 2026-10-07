#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static double qq(std::vector<double>u,double p){ if(u.empty())return 0; std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
void an_res(const char *fn="val-search.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("search"); ValSearchHit *h=nullptr; t->SetBranchAddress("h",&h);
  struct A{ std::vector<double> rx,ry,rz,c2; long n=0,acc=0; };
  std::map<int,A> L; long nfill=0, nkal=0;
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->mc_match||!h->had_kalman) continue; ++nkal;
    if(h->residual_x<-900) continue; ++nfill;
    int g = h->layer<4?0 : (h->layer<10?1 : (h->layer<16?2 : 3));
    A&a=L[g]; ++a.n; if(h->accepted) ++a.acc;
    a.rx.push_back(std::fabs(h->residual_x)); a.ry.push_back(std::fabs(h->residual_y));
    a.rz.push_back(std::fabs(h->residual_z)); a.c2.push_back(h->chi2);
  }
  printf("\nMODULE-FRAME RESIDUALS of MC-matched hits that reached the Kalman.\n");
  printf("%ld of %ld have the residual filled.\n",nfill,nkal);
  printf("x = across strip (precise / phi), y = along strip (coarse), z = off the plane (must be ~0).\n\n");
  const char* gn[4]={"PixB 0-3","TOB PS 4-9","TOB 2S 10-15","fwd pix / disks"};
  printf("  %-16s %7s | %9s %9s | %9s %9s | %9s | %8s %7s\n",
         "region","n","|res_x| med","p90","|res_y| med","p90","|res_z| med","chi2 med","acc %");
  for(auto&kv:L){ A&a=kv.second; if(a.n<30) continue;
    printf("  %-16s %7ld | %9.5f %9.5f | %9.5f %9.5f | %9.2e | %8.3f %7.1f\n",
      gn[kv.first],a.n,qq(a.rx,50),qq(a.rx,90),qq(a.ry,50),qq(a.ry,90),qq(a.rz,50),
      qq(a.c2,50),100.*a.acc/a.n); }
  printf("\n  Units are cm. Pixel pitch ~ 25-100 um = 0.0025-0.01 cm;\n");
  printf("  strip pitch ~ 90-100 um; a TOB 2S strip half-length is 2.5 cm.\n");
}
