#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// Is the phi deficit a property of the REGION or of the DEPTH in the search?
//
// Last night's localisation -- pixels short, strips clean at the same eta -- is
// confounded: the inward search reaches the strips EARLY in a candidate's life
// and the pixel barrel LATE. Same trap as the material-vs-path-length question.
// `step` (TrCandState::step) is already in the trace, so this needs no new run.
static double qtl(std::vector<double> u,double p){
  if(u.empty())return 0; std::sort(u.begin(),u.end());
  return u[(size_t)(p/100.*(u.size()-1))]; }
static double core(std::vector<double>&u){
  return u.size()<50?-1:(qtl(u,75)-qtl(u,25))/1.349; }

void an_phidepth(const char *fn="val-phi-in.root"){
  const int NS=4; const int sb[NS+1]={0,2,4,6,99};
  const char* sn[NS]={"step 0-1","step 2-3","step 4-5","step >=6"};
  struct Reg{const char*n;int lo,hi;};
  const Reg R[]={{"PixB 0-3",0,3},{"TBPS 4-9",4,9},{"FwdPix 16-27",16,27}};
  const int NR=3;
  // Second question: is the ETA dependence itself just DEPTH in disguise?
  // Bin by eta at FIXED step, all regions pooled.
  const int NEE=4; const double ee[NEE+1]={0,0.8,1.6,2.0,9};
  const char* een[NEE]={"<0.8","0.8-1.6","1.6-2.0",">2.0"};
  std::vector<double> ed[NEE][NS], es[NEE][NS], eh[NEE][NS];
  std::vector<double> dp[NR][NS], sp[NR][NS], hp[NR][NS];
  TFile f(fn); TTree*t=(TTree*)f.Get("search");
  ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->mc_match||h->sigma_phi_trk<=0||h->dphi_s<-900||h->hit_phi_sigma<0) continue;
    if(h->step<0||h->eta<-900) continue;
    double ae=std::fabs(h->eta); if(ae<1.6||ae>2.5) continue;
    if(h->dq_s<-900||h->hit_q_half_len<-900) continue;
    if(std::fabs(h->dq_s) > std::max(3.0,5.0*h->hit_q_half_len)) continue;  // clamp
    int r=-1; for(int k=0;k<NR;++k) if(h->layer>=R[k].lo&&h->layer<=R[k].hi){r=k;break;}
    if(r<0) continue;
    int st=-1; for(int b=0;b<NS;++b) if(h->step>=sb[b]&&h->step<sb[b+1]){st=b;break;}
    if(st<0) continue;
    dp[r][st].push_back(h->dphi_s); sp[r][st].push_back(h->sigma_phi_trk);
    hp[r][st].push_back(h->hit_phi_sigma); }
  // second pass for the eta-at-fixed-step table (no eta restriction)
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->mc_match||h->sigma_phi_trk<=0||h->dphi_s<-900||h->hit_phi_sigma<0) continue;
    if(h->step<0||h->eta<-900||h->dq_s<-900||h->hit_q_half_len<-900) continue;
    if(std::fabs(h->dq_s) > std::max(3.0,5.0*h->hit_q_half_len)) continue;
    double ae=std::fabs(h->eta); int e=-1;
    for(int b=0;b<NEE;++b) if(ae>=ee[b]&&ae<ee[b+1]){e=b;break;}
    int st=-1; for(int b=0;b<NS;++b) if(h->step>=sb[b]&&h->step<sb[b+1]){st=b;break;}
    if(e<0||st<0) continue;
    ed[e][st].push_back(h->dphi_s); es[e][st].push_back(h->sigma_phi_trk);
    eh[e][st].push_back(h->hit_phi_sigma); }

  printf("\n===== IS IT THE REGION OR THE DEPTH? =====\n");
  printf("phi short by, |eta| 1.6-2.5, clamped truth match. If the deficit follows\n");
  printf("the REGION it stays put down each column; if it follows the DEPTH it rises\n");
  printf("along each row regardless of region.\n\n");
  printf("  %-14s", "region");
  for(int b=0;b<NS;++b) printf(" |%14s", sn[b]);
  printf("\n  %-14s","");
  for(int b=0;b<NS;++b) printf(" |%8s%6s","n","short");
  printf("\n");
  for(int r=0;r<NR;++r){
    printf("  %-14s", R[r].n);
    for(int b=0;b<NS;++b){
      double c=core(dp[r][b]);
      if(c<0){ printf(" |%8s%6s","-","-"); continue; }
      double st=qtl(sp[r][b],50), sh=qtl(hp[r][b],50);
      double exp=std::sqrt(st*st+sh*sh)/st;
      printf(" |%8zu%6.2f", dp[r][b].size(), (c/st)/exp);
    }
    printf("\n");
  }
  printf("\n  --- and is the ETA dependence just DEPTH in disguise? ---\n");
  printf("  all regions pooled. If eta matters in its own right, each COLUMN rises.\n\n");
  printf("  %-10s", "|eta|");
  for(int b=0;b<NS;++b) printf(" |%14s", sn[b]);
  printf("\n  %-10s","");
  for(int b=0;b<NS;++b) printf(" |%8s%6s","n","short");
  printf("\n");
  for(int e=0;e<NEE;++e){
    printf("  %-10s", een[e]);
    for(int b=0;b<NS;++b){
      double c=core(ed[e][b]);
      if(c<0){ printf(" |%8s%6s","-","-"); continue; }
      double st=qtl(es[e][b],50), sh=qtl(eh[e][b],50);
      double exp=std::sqrt(st*st+sh*sh)/st;
      printf(" |%8zu%6.2f", ed[e][b].size(), (c/st)/exp);
    }
    printf("\n");
  }
  printf("\n");
}
