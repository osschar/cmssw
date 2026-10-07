#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// Does the surface reference explain the "5x already at handover" leg of the
// 5.30 x 4.21 = 22.3 decomposition?
//
// That decomposition used sigma_q_trk -- the WINDOW covariance -- for "quoted",
// so the missing surface term is in both legs, but not symmetrically:
//   * the 4.21x accumulated IN the search is a ratio of ratios on the SAME
//     track, and theta barely changes along a track, so cosh^2(eta) largely
//     cancels -> that leg should survive;
//   * the 5.30x at HANDOVER is absolute, so it carries the full factor -> it
//     should shrink, and if it shrinks to ~1 the backward fit is exonerated.
// step = TrCandState::step, how many layers into THIS search. step<=1 is the
// handover end, step>=4 is deep.
static double qtl(std::vector<double> u,double p){
  if(u.empty())return 0; std::sort(u.begin(),u.end());
  return u[(size_t)(p/100.*(u.size()-1))]; }

void an_surfq_step(const char*f0="val-surfq-0.root", const char*f1="val-surfq-1.root"){
  const double eb[]={0,0.8,1.6,2.6}; const int NE=3;
  const char* en[NE]={"< 0.8","0.8-1.6","1.6-2.5"};
  const int NS=3; const char* sn[NS]={"step<=1 (handover)","step 2-3","step>=4 (deep)"};
  std::vector<double> dq[2][NE][NS], sg[2][NE][NS];
  const char* fn[2]={f0,f1};
  for(int m=0;m<2;++m){
    TFile f(fn[m]); TTree*t=(TTree*)f.Get("search");
    ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
      if(h->layer<0||h->layer>=4||!h->mc_match) continue;     // pixel barrel
      if(h->sigma_q_trk<=0||h->dq_s<-900||h->eta<-900||h->step<0) continue;
      double ae=std::fabs(h->eta); int e=-1;
      for(int b=0;b<NE;++b) if(ae>=eb[b]&&ae<eb[b+1]){e=b;break;}
      if(e<0) continue;
      int st = h->step<=1 ? 0 : (h->step<=3 ? 1 : 2);
      dq[m][e][st].push_back(h->dq_s); sg[m][e][st].push_back(h->sigma_q_trk); } }

  printf("\n===== ratio (core/quoted) vs DEPTH IN THE SEARCH, pixel barrel =====\n");
  printf("If the 'already too tight at handover' leg is really the missing surface\n");
  printf("term, the step<=1 column should collapse when the correction is on.\n\n");
  for(int e=0;e<NE;++e){
    printf("  |eta| %s\n", en[e]);
    printf("    %-22s %8s | %9s %9s | %7s\n","depth","n","ratio OFF","ratio ON","factor");
    for(int st=0;st<NS;++st){
      if(dq[0][e][st].size()<40||dq[1][e][st].size()<40){
        printf("    %-22s %8zu | %9s %9s |\n", sn[st], dq[0][e][st].size(),"-","-"); continue; }
      double c0=(qtl(dq[0][e][st],75)-qtl(dq[0][e][st],25))/1.349, q0=qtl(sg[0][e][st],50);
      double c1=(qtl(dq[1][e][st],75)-qtl(dq[1][e][st],25))/1.349, q1=qtl(sg[1][e][st],50);
      printf("    %-22s %8zu | %9.2f %9.2f | %7.2f\n",
             sn[st], dq[0][e][st].size(), c0/q0, c1/q1, (c0/q0)/(c1/q1)); }
    // the "accumulated in the search" leg: deep / handover, each mode separately
    if(dq[0][e][0].size()>=40&&dq[0][e][2].size()>=40&&dq[1][e][0].size()>=40&&dq[1][e][2].size()>=40){
      auto R=[&](int m,int st){ return ((qtl(dq[m][e][st],75)-qtl(dq[m][e][st],25))/1.349)/qtl(sg[m][e][st],50); };
      printf("    -> accumulated deep/handover:  OFF %.2f   ON %.2f\n", R(0,2)/R(0,0), R(1,2)/R(1,0)); }
    printf("\n"); }
  printf("  Reading: the 'factor' column is what the surface reference removed at that\n");
  printf("  depth. The accumulated line is the leg that should SURVIVE, because theta\n");
  printf("  is nearly constant along a track so cosh^2(eta) cancels between depths.\n\n");
}
