#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static double qq(std::vector<double>u,double p){ if(u.empty())return 0; std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
struct A { std::vector<double> pull,c2,sq,sh,dq; long nk=0,na=0; };
void an_lay(const char *fn="val-search.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("search"); ValSearchHit *h=nullptr; t->SetBranchAddress("h",&h);
  std::map<int,A> L; A all, clamp;
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->mc_match) continue;
    if(!(h->sigma_q_trk>0 && h->hit_q_half_len>-900 && h->dq>-900)) continue;
    double sh = h->hit_q_half_len/std::sqrt(3.0), sq=h->sigma_q_trk;
    double st = std::sqrt(sq*sq+sh*sh), p=std::fabs(h->dq)/st;
    A &a=L[h->layer]; a.pull.push_back(p); a.sq.push_back(sq); a.sh.push_back(sh); a.dq.push_back(std::fabs(h->dq));
    if(h->had_kalman){ ++a.nk; if(h->accepted) ++a.na; if(h->chi2>-900) a.c2.push_back(h->chi2); }
    all.pull.push_back(p);
    // GEOMETRIC CLAMP: a hit further in q than a few times its own extent plus a
    // few track sigmas cannot be a measurement of this track -- it shares the
    // sim label and the layer, nothing more. 5 sigma_tot is generous.
    if (p < 5.0){ clamp.pull.push_back(p); if(h->chi2>-900&&h->had_kalman){clamp.c2.push_back(h->chi2); ++clamp.nk; if(h->accepted) ++clamp.na;} }
  }
  printf("\nMC-MATCHED HITS, inward search, by layer.\n");
  printf("Reference for |dq|/sigma_tot: a GAUSSIAN residual gives median 0.674; a residual\n");
  printf("dominated by uniform spread along the strip gives 0.866. So expect 0.67-0.87.\n");
  printf("sig_trk / sig_hit say which term dominates: sig_hit = hit_q_half_len/sqrt(3).\n\n");
  printf("  %5s %7s | %9s %9s | %8s %8s %8s | %8s %7s\n",
         "layer","n","sig_trk cm","sig_hit cm","med pull","p90","p99","chi2 med","acc %");
  for(auto &kv:L){ A&a=kv.second; if(a.pull.size()<30) continue;
    printf("  %5d %7zu | %9.4f %9.4f | %8.3f %8.2f %8.1f | %8.3f %7.1f\n",
      kv.first,a.pull.size(),qq(a.sq,50),qq(a.sh,50),
      qq(a.pull,50),qq(a.pull,90),qq(a.pull,99),
      qq(a.c2,50), a.nk? 100.*a.na/a.nk : -1.0); }
  printf("\n  %-34s n=%zu  med %.3f  p90 %.2f\n","ALL, no clamp",all.pull.size(),qq(all.pull,50),qq(all.pull,90));
  printf("  %-34s n=%zu  med %.3f  p90 %.2f   chi2 med %.3f   accepted %.1f%%\n",
    "ALL, |dq| < 5 sigma_tot",clamp.pull.size(),qq(clamp.pull,50),qq(clamp.pull,90),
    qq(clamp.c2,50), clamp.nk? 100.*clamp.na/clamp.nk : -1.0);
}
