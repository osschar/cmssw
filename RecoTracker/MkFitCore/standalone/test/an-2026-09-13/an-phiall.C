#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// PHI covariance across the whole detector -- localising the eta-rising
// component. phi is the only coordinate measurable everywhere (cluster spans are
// 2.2-3.3 in the precise direction in every region; q is 100 % single-cell in
// every strip region and therefore unmeasurable there).
//
// Reports for every region x eta cell:
//   short  = (measured core / sigma_phi_trk) / expected, with the HIT term
//            divided out: expected = sqrt(sig_trk^2 + sig_hit^2)/sig_trk
//   trkf   = sig_trk^2 / (sig_trk^2 + sig_hit^2) -- the TRACK's share of the
//            residual. A cell with trkf near 0 says NOTHING about the
//            covariance however confident 'short' looks; printed so an
//            uninformative cell is visibly uninformative.
static double qtl(std::vector<double> u,double p){
  if(u.empty())return 0; std::sort(u.begin(),u.end());
  return u[(size_t)(p/100.*(u.size()-1))]; }
static double core(std::vector<double>&u){
  return u.size()<50?-1:(qtl(u,75)-qtl(u,25))/1.349; }

struct Reg { const char *name; int lo, hi, parity; };
static const Reg R[] = {
  {"PixB 0-3",       0,  3, -1}, {"TBPS-P 4/6/8",  4,  9,  0},
  {"TBPS-S 5/7/9",   4,  9,  1}, {"TOB2S 10-15",  10, 15, -1},
  {"FwdPix 16-27",  16, 27, -1}, {"TEC 28-37",    28, 37, -1},
};
static const int NR = sizeof(R)/sizeof(Reg);
static const double eb[] = {0, 0.8, 1.6, 2.6};
static const int NE = 3;
static const char *en[NE] = {"<0.8", "0.8-1.6", "1.6-2.5"};
static bool g_clamp = false;

static void one(const char *tag, const char *fn){
  TFile f(fn); TTree*t=(TTree*)f.Get("search");
  if(!t){ printf("  %s: no tree in %s\n", tag, fn); return; }
  ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
  std::vector<double> dp[NR][NE], sp[NR][NE], hp[NR][NE];
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->mc_match||h->sigma_phi_trk<=0||h->dphi_s<-900||h->eta<-900) continue;
    if(h->hit_phi_sigma<0) continue;
    // GEOMETRIC CLAMP on the truth match, and it is ORTHOGONAL to what is being
    // measured: mc_match only tests sim_lbl == hit_lbl and the layer, so a hit
    // of the right particle but metres away in z passes. Cutting on |dq| kills
    // those without touching phi, so the phi measurement stays unbiased.
    // Generous by design: 3 cm, or 5x the hit's own q extent for long strips.
    if(g_clamp){
      if(h->dq_s<-900||h->hit_q_half_len<-900) continue;
      const double lim = std::max(3.0, 5.0*h->hit_q_half_len);
      if(std::fabs(h->dq_s) > lim) continue;
    }
    int r=-1;
    for(int k=0;k<NR;++k)
      if(h->layer>=R[k].lo&&h->layer<=R[k].hi&&(R[k].parity<0||(h->layer&1)==R[k].parity)){r=k;break;}
    if(r<0) continue;
    double ae=std::fabs(h->eta); int e=-1;
    for(int b=0;b<NE;++b) if(ae>=eb[b]&&ae<eb[b+1]){e=b;break;}
    if(e<0) continue;
    dp[r][e].push_back(h->dphi_s); sp[r][e].push_back(h->sigma_phi_trk);
    hp[r][e].push_back(h->hit_phi_sigma); }

  printf("\n--- %s ---\n", tag);
  printf("  %-14s", "region");
  for(int e=0;e<NE;++e) printf(" |%18s", en[e]);
  printf("\n  %-14s","");
  for(int e=0;e<NE;++e) printf(" |%7s%7s%5s","n","short","trkf");
  printf("\n");
  for(int r=0;r<NR;++r){
    printf("  %-14s", R[r].name);
    for(int e=0;e<NE;++e){
      double c=core(dp[r][e]);
      if(c<0){ printf(" |%7s%7s%5s","-","-","-"); continue; }
      double st=qtl(sp[r][e],50), sh=qtl(hp[r][e],50);
      double exp=std::sqrt(st*st+sh*sh)/st;
      double trkf=st*st/(st*st+sh*sh);
      printf(" |%7zu%7.2f%5.2f", dp[r][e].size(), (c/st)/exp, trkf);
    }
    printf("\n");
  }
}

void an_phiall(bool clamp = false){
  g_clamp = clamp;
  printf("\n%s\n", clamp ? "########## WITH the geometric |dq| clamp ##########"
                           : "########## NO clamp (raw mc_match) ##########");
  printf("\n===== PHI COVARIANCE ACROSS THE DETECTOR =====\n");
  printf("short = how many times too small the phi covariance is, hit error\n");
  printf("divided out. 1.0 = correct. trkf = track's share of the residual\n");
  printf("variance; a cell with small trkf is UNINFORMATIVE regardless of 'short'.\n");
  one("REAL OUTWARD  (ProcessEventStd, real cmssw seeds)", "val-phi-out.root");
  one("CHOPPED pT5 INWARD (ProcessEventHlt)",              "val-phi-in.root");
  one("SIM-SEEDED OUTWARD (invented seed covariance)",     "val-phi-sim1.root");
  printf("\n");
}
