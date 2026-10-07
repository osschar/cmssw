#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// The ratio core/sigma_trk is BIASED HIGH: the measured residual carries the hit
// error too, var = sig_trk^2 + sig_hit^2. Divide it out and ask what is really
// left, separately for q and phi -- phi being the coordinate the surface
// reference never touched, so it is the clean baseline.
static double qtl(std::vector<double> u,double p){
  if(u.empty())return 0; std::sort(u.begin(),u.end());
  return u[(size_t)(p/100.*(u.size()-1))]; }
static double core(std::vector<double>&u){
  return u.size()<40?-1:(qtl(u,75)-qtl(u,25))/1.349; }

void an_hitterm(const char *fn="val-surfq-1.root"){
  const double eb[]={0,0.8,1.6,2.0,9}; const int NE=4;
  const char* en[NE]={"< 0.8","0.8-1.6","1.6-2.0","> 2.0"};
  std::vector<double> dq[NE],sq[NE],hq[NE],dp[NE],sp[NE],hp[NE];
  TFile f(fn); TTree*t=(TTree*)f.Get("search");
  ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(h->layer<0||h->layer>=4||!h->mc_match) continue;
    if(h->sigma_q_trk<=0||h->dq_s<-900||h->eta<-900||h->hit_q_half_len<-900) continue;
    double ae=std::fabs(h->eta); int e=-1;
    for(int b=0;b<NE;++b) if(ae>=eb[b]&&ae<eb[b+1]){e=b;break;}
    if(e<0) continue;
    dq[e].push_back(h->dq_s); sq[e].push_back(h->sigma_q_trk);
    hq[e].push_back(h->hit_q_half_len/3.0);   // PIXELS: HitStructures.cc uses hl_fac = 3 (sqrt(3) is for strips)
    if(h->dphi_s>-900&&h->sigma_phi_trk>0&&h->hit_phi_sigma>=0){
      dp[e].push_back(h->dphi_s); sp[e].push_back(h->sigma_phi_trk);
      hp[e].push_back(h->hit_phi_sigma); } }

  printf("\n===== WHAT IS LEFT AFTER DIVIDING OUT THE HIT ERROR =====\n");
  printf("Pixel barrel, MC-matched, surface reference ON.\n");
  printf("expected = sqrt(sig_trk^2 + sig_hit^2)/sig_trk -- what the ratio would be\n");
  printf("if the track covariance were exactly right. 'short by' = measured/expected.\n\n");
  printf("  %-9s %7s | %8s %8s %7s %7s %8s | %8s %8s %7s %7s %8s\n",
         "|eta|","n","sig_trk","sig_hit","ratio","expect","q SHORT",
         "sigP_trk","sigP_hit","ratio","expect","phi SHORT");
  for(int e=0;e<NE;++e){
    double cq=core(dq[e]); if(cq<0) continue;
    double st=qtl(sq[e],50), sh=qtl(hq[e],50);
    double exp_q=std::sqrt(st*st+sh*sh)/st;
    double cp=core(dp[e]), sp50=qtl(sp[e],50);
    // hit phi sigma: use the pixel's precise pitch, ~25um/sqrt(12), at this radius.
    // taken from the data instead: no per-hit phi extent is stored, so quote the
    // phi columns WITHOUT a hit term and say so.
    double shp = qtl(hp[e],50);
    double exp_p = sp50 > 0 ? std::sqrt(sp50*sp50 + shp*shp)/sp50 : 1.0;
    printf("  %-9s %7zu | %8.1f %8.1f %7.2f %7.3f %8.2f | %8.3f %8.3f %7.2f %7.3f %8.2f\n",
           en[e], dq[e].size(), st*1e4, sh*1e4, cq/st, exp_q, (cq/st)/exp_q,
           sp50*1e3, shp*1e3, cp/sp50, exp_p, (cp/sp50)/exp_p);
  }
  printf("\n  q sigmas in um, phi sigmas in mrad. BOTH coordinates are now hit-corrected:\n");
  printf("  hit_phi_sigma is projected from the hit covariance in val_search_event,\n");
  printf("  since LayerOfHits carries no per-hit phi extent (it uses a flat\n");
  printf("  HIT_PHI_HALF_EXTENT = 0.0123 rad).\n\n");
}
