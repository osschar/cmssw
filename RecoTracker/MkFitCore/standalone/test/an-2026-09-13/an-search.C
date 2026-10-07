#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
static double q(std::vector<double>&u,double p){ if(u.empty())return 0; std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
static double qs(std::vector<double>&u,double p){ if(u.empty())return 0; return u[(size_t)(p/100.*(u.size()-1))]; }
void an_search(const char *fn="val-search.root"){
  TFile f(fn); TTree *t=(TTree*)f.Get("search"); ValSearchHit *h=nullptr; t->SetBranchAddress("h",&h);
  long n[2]={0,0}, npre[2]={0,0}, npq[2]={0,0}, nkal[2]={0,0}, nacc[2]={0,0}, nout=0;
  std::vector<double> c2[2], pq[2], pphi[2], dqa[2];
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if (h->is_outward) { ++nout; continue; }
    int m = h->mc_match?1:0;
    ++n[m];
    if (h->passed_preselect) ++npre[m];
    if (h->passed_pqueue)    ++npq[m];
    if (h->had_kalman){ ++nkal[m]; if(h->accepted) ++nacc[m];
      if (h->chi2>-900) c2[m].push_back(h->chi2); }
    if (h->sigma_q_trk>0 && h->hit_q_half_len>-900 && h->dq>-900){
      double st = std::sqrt(h->sigma_q_trk*h->sigma_q_trk
                          + h->hit_q_half_len*h->hit_q_half_len/3.0);
      pq[m].push_back(std::fabs(h->dq)/st); dqa[m].push_back(std::fabs(h->dq)); }
    if (h->sigma_phi_trk>0 && h->dphi>-900) pphi[m].push_back(std::fabs(h->dphi)/h->sigma_phi_trk);
  }
  printf("\nINWARD (backward) search, T5 after backward fit. %lld records, %ld outward (dropped).\n",
         t->GetEntries(), nout);
  printf("\nFUNNEL, split on MC truth of the hit:\n");
  printf("  %-24s %10s %10s   %s\n","stage","MC-matched","other","survival of MC-matched");
  printf("  %-24s %10ld %10ld\n","scanned",n[1],n[0]);
  printf("  %-24s %10ld %10ld   %6.1f%%\n","passed pre-selection",npre[1],npre[0],100.*npre[1]/std::max(1L,n[1]));
  printf("  %-24s %10ld %10ld   %6.1f%%\n","survived pqueue",npq[1],npq[0],100.*npq[1]/std::max(1L,npre[1]));
  printf("  %-24s %10ld %10ld   %6.1f%%\n","reached Kalman",nkal[1],nkal[0],100.*nkal[1]/std::max(1L,npq[1]));
  printf("  %-24s %10ld %10ld   %6.1f%%  <== accepted / reached Kalman\n","ACCEPTED",nacc[1],nacc[0],100.*nacc[1]/std::max(1L,nkal[1]));
  const char *nm[2]={"other","MC-matched"};
  printf("\nPULLS, |residual| / sigma. Expectation: half-normal, median 0.674, p90 1.645.\n");
  printf("  q uses sigma_tot = sqrt(sigma_q_trk^2 + (hit_q_half_len)^2/3), i.e. track AND hit extent.\n\n");
  printf("  %-12s %8s %8s %8s %8s | %8s %8s %8s\n","sample","n","|dq|/sig","  p90","  p99","|dphi|/s","  p90","  p99");
  for(int m=1;m>=0;--m)
    printf("  %-12s %8zu %8.3g %8.3g %8.3g | %8.3g %8.3g %8.3g\n",nm[m],pq[m].size(),
      q(pq[m],50),qs(pq[m],90),qs(pq[m],99), q(pphi[m],50),qs(pphi[m],90),qs(pphi[m],99));
  printf("\nKALMAN chi2 (2-D in the module plane; expectation 2, median of chi2_2 = 1.386):\n");
  printf("  %-12s %8s %8s %8s %8s %8s\n","sample","n","median","p90","p99","max");
  for(int m=1;m>=0;--m)
    printf("  %-12s %8zu %8.3g %8.3g %8.3g %8.3g\n",nm[m],c2[m].size(),
      q(c2[m],50),qs(c2[m],90),qs(c2[m],99), c2[m].empty()?0:c2[m].back());
  printf("\n|dq| itself [cm], MC-matched: median %.3g  p90 %.3g  p99 %.3g  max %.3g\n",
    q(dqa[1],50),qs(dqa[1],90),qs(dqa[1],99), dqa[1].empty()?0:dqa[1].back());
}
