#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// Per-region A/B of the surface-referenced q window.
//
// The correction runs in MkBins, before any hit is known, so it can only use the
// LAYER surface -- I take the barrel normal as radial. That is exact for the
// pixel barrel (IT tilt is measured zero) but WRONG for tilted TBPS, where the
// modules face the IP and a track hits them near-normally while its incidence to
// the radial direction is large. There the cylinder projection should OVER-widen.
// This macro is the check: watch for ratio dropping well below 1 in layers 4-9,
// and for hits scanned blowing up there.
static double qtl(std::vector<double> u, double p){
  if(u.empty()) return 0; std::sort(u.begin(),u.end());
  return u[(size_t)(p/100.*(u.size()-1))]; }

// TBPS P and S must be SEPARATE: both are tilted, but the P macro-pixel's q
// extent is ~0.04 cm against the S strip's ~0.6 cm, so the hit term swamps the
// track term in S and does NOT in P. Lumping 4-9 averages the one informative
// sub-layer away.  parity: even = P (post-10781bd48eb), odd = S.
struct R { const char *name; int lo, hi, parity; };   // parity <0 = both
static const R REG[] = {
  {"PixB 0-3",          0,  3, -1},
  {"TBPS-P 4/6/8 tilt", 4,  9,  0},
  {"TBPS-S 5/7/9 tilt", 4,  9,  1},
  {"TOB 2S 10-15",     10, 15, -1},
  {"fwd pix 16-27",    16, 27, -1},
  {"TEC 28-37",        28, 37, -1},
};
static const int NR = sizeof(REG)/sizeof(R);

void an_surfq_reg(const char *f0="val-surfq-0.root", const char *f1="val-surfq-1.root"){
  std::vector<double> dq[2][NR], sg[2][NR];
  long nsc[2][NR]={}, ntr[2][NR]={}, npass[2][NR]={};
  const char *fn[2]={f0,f1};
  for(int m=0;m<2;++m){
    TFile f(fn[m]); TTree*t=(TTree*)f.Get("search");
    ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
      int r=-1;
      for(int k=0;k<NR;++k)
        if(h->layer>=REG[k].lo && h->layer<=REG[k].hi &&
           (REG[k].parity<0 || (h->layer&1)==REG[k].parity)) {r=k;break;}
      if(r<0) continue;
      ++nsc[m][r];
      if(!h->mc_match||h->sigma_q_trk<=0||h->dq_s<-900) continue;
      dq[m][r].push_back(h->dq_s); sg[m][r].push_back(h->sigma_q_trk);
      ++ntr[m][r]; if(h->passed_preselect) ++npass[m][r]; } }

  printf("\n=========== SURFACE-REFERENCED q WINDOW, BY REGION ===========\n");
  printf("ratio = (IQR/1.349 of signed dq) / (median sigma_q_trk). 1.0 = correctly\n");
  printf("sized. Well BELOW 1 after the change = over-widened, which is what the\n");
  printf("radial-normal approximation is expected to do in TILTED TBPS.\n\n");
  printf("  %-20s %8s | %7s %7s | %8s %8s | %10s %10s\n",
         "region","n true","ratio 0","ratio 1","presel 0","presel 1","scanned 0","scanned 1");
  for(int r=0;r<NR;++r){
    if(dq[0][r].size()<40||dq[1][r].size()<40){
      printf("  %-20s %8zu | %7s %7s | %8s %8s | %10ld %10ld\n",
             REG[r].name, dq[0][r].size(), "-","-","-","-", nsc[0][r], nsc[1][r]);
      continue; }
    double c0=(qtl(dq[0][r],75)-qtl(dq[0][r],25))/1.349, q0=qtl(sg[0][r],50);
    double c1=(qtl(dq[1][r],75)-qtl(dq[1][r],25))/1.349, q1=qtl(sg[1][r],50);
    printf("  %-20s %8zu | %7.2f %7.2f | %7.1f%% %7.1f%% | %10ld %10ld  %+5.1f%%\n",
           REG[r].name, dq[0][r].size(), c0/q0, c1/q1,
           100.*npass[0][r]/ntr[0][r], 100.*npass[1][r]/ntr[1][r],
           nsc[0][r], nsc[1][r], 100.*(nsc[1][r]-nsc[0][r])/std::max(1L,nsc[0][r])); }
  long t0=0,t1=0; for(int r=0;r<NR;++r){t0+=nsc[0][r];t1+=nsc[1][r];}
  printf("\n  TOTAL scanned %ld -> %ld  (%+.1f%%)\n", t0,t1, 100.*(t1-t0)/t0);
  printf("\n  NOTE: in strip layers the cut is carried by the HIT extent (a 2.5 cm\n");
  printf("  strip), not by the track term, so pre-selection there is insensitive to\n");
  printf("  this change by construction. The number to read is 'scanned'.\n\n");
}
