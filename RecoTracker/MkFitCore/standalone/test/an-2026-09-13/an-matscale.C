#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// Does scaling radL close the isotropic covariance deficit left after the
// surface reference? ratio = measured core / quoted sigma; 1.0 = correct.
// Reported for BOTH coordinates: multiple scattering is isotropic, so if this is
// the mechanism, dq and dphi must improve TOGETHER.
static double qtl(std::vector<double> u,double p){
  if(u.empty())return 0; std::sort(u.begin(),u.end());
  return u[(size_t)(p/100.*(u.size()-1))]; }
static double core(std::vector<double>&u){
  return u.size()<40?-1:(qtl(u,75)-qtl(u,25))/1.349; }

void an_matscale(){
  const char* sc[]={"1.0","10.0","100.0","1000.0"}; const int NM=4;
  const double eb[]={0,0.8,1.6,2.6}; const int NE=3;
  const char* en[NE]={"< 0.8","0.8-1.6","1.6-2.5"};
  printf("\n===== radL SCALE vs the residual covariance deficit =====\n");
  printf("Pixel barrel, MC-matched, surface reference ON. ratio = core/quoted.\n");
  printf("If under-counted scattering is the mechanism, BOTH columns walk to 1.\n\n");
  printf("  %-6s %-9s %8s | %9s %9s %7s | %9s %9s %7s\n",
         "scale","|eta|","n","dq CORE","quoted","ratio","dphi CORE","quoted","ratio");
  for(int m=0;m<NM;++m){
    TFile f(Form("val-mat-%s.root",sc[m]));
    TTree*t=(TTree*)f.Get("search"); if(!t){printf("  %-6s -- no tree\n",sc[m]);continue;}
    ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
    std::vector<double> dq[NE],sq[NE],dp[NE],sp[NE];
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
      if(h->layer<0||h->layer>=4||!h->mc_match) continue;
      if(h->sigma_q_trk<=0||h->dq_s<-900||h->eta<-900) continue;
      double ae=std::fabs(h->eta); int e=-1;
      for(int b=0;b<NE;++b) if(ae>=eb[b]&&ae<eb[b+1]){e=b;break;}
      if(e<0) continue;
      dq[e].push_back(h->dq_s); sq[e].push_back(h->sigma_q_trk);
      if(h->dphi_s>-900&&h->sigma_phi_trk>0){ dp[e].push_back(h->dphi_s); sp[e].push_back(h->sigma_phi_trk);} }
    for(int e=0;e<NE;++e){
      double cq=core(dq[e]); if(cq<0) continue;
      double qq=qtl(sq[e],50), cp=core(dp[e]), pq=qtl(sp[e],50);
      printf("  %-6s %-9s %8zu | %9.1f %9.1f %7.2f | %9.4f %9.4f %7.2f\n",
             m?"":sc[m], en[e], dq[e].size(), cq*1e4, qq*1e4, cq/qq,
             cp*1e3, pq*1e3, cp>0?cp/pq:0.0);
      if(e==0) printf("");
    }
    printf("  %-6s\n", sc[m]); }
  printf("\n  (q in um, phi in mrad)\n");
  printf("  CORE is the MEASURED residual width. If it SHRINKS with material, the\n");
  printf("  larger gain is letting hits correct unmodelled systematics -- i.e. the\n");
  printf("  deficit closes from both ends, not just by inflating sigma.\n\n");
}
