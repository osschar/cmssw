#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
// THE GEOMETRY SIDE'S FALSIFIABLE PREDICTION (geo-stuff-9d, ANSWERS section 12a):
// out-of-envelope material rises monotonically with |z| across the barrel, 21%
// -> 52%. A track at fixed eta moves to larger |z| as it moves outward, so the
// fraction of material we OMIT grows layer by layer -- a compounding mechanism.
// Prediction: the PER-STEP degradation rate should itself worsen with |eta|,
// because a higher-eta track reaches the high-omission region at an earlier
// layer index. If the rate is FLAT in eta, the services are not the driver.
//
// Fit log(short) vs step per eta bin; the slope IS the per-step rate.
static double qtl(std::vector<double> u,double p){
  if(u.empty())return 0; std::sort(u.begin(),u.end());
  return u[(size_t)(p/100.*(u.size()-1))]; }
static double core(std::vector<double>&u){
  return u.size()<50?-1:(qtl(u,75)-qtl(u,25))/1.349; }

void an_phirate(const char *fn="val-phi-in.root", int mode=0){
  // mode 0 = split by |eta|; mode 1 = split by REGION, which is the direct
  // test of "the rate does not care how much material the step crossed":
  // per-layer budgets differ 2.4x between PixB (0.0185) and TBPS (0.045).
  const int NST=8;                       // step 0..7+
  const int NEE=4;
  const double ee[NEE+1]={0,0.8,1.6,2.0,9};
  const char* een_eta[NEE]={"<0.8","0.8-1.6","1.6-2.0",">2.0"};
  const char* een_reg[NEE]={"PixB 0-3","TBPS 4-9","TOB2S 10-15","FwdPix 16-27"};
  const int rlo[NEE]={0,4,10,16}, rhi[NEE]={3,9,15,27};
  const char** een = mode ? een_reg : een_eta;
  std::vector<double> d[NEE][NST], sg[NEE][NST], hh[NEE][NST];
  TFile f(fn); TTree*t=(TTree*)f.Get("search");
  ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->mc_match||h->sigma_phi_trk<=0||h->dphi_s<-900||h->hit_phi_sigma<0) continue;
    if(h->step<0||h->eta<-900||h->dq_s<-900||h->hit_q_half_len<-900) continue;
    if(std::fabs(h->dq_s) > std::max(3.0,5.0*h->hit_q_half_len)) continue;
    int e=-1;
    if(mode){ for(int b=0;b<NEE;++b) if(h->layer>=rlo[b]&&h->layer<=rhi[b]){e=b;break;} }
    else { double ae=std::fabs(h->eta);
           for(int b=0;b<NEE;++b) if(ae>=ee[b]&&ae<ee[b+1]){e=b;break;} }
    if(e<0) continue;
    int st = std::min(h->step, NST-1);
    d[e][st].push_back(h->dphi_s); sg[e][st].push_back(h->sigma_phi_trk);
    hh[e][st].push_back(h->hit_phi_sigma); }

  printf("\n===== PER-STEP DEGRADATION RATE vs %s =====\n", mode?"REGION":"|eta|");
  printf("phi short by, per step. Slope from an unweighted least-squares fit of\n");
  printf("log(short) vs step over the populated steps; rate = exp(slope) is the\n");
  printf("per-layer multiplier.\n\n");
  printf("  %-13s", mode?"region":"|eta|");
  for(int k=0;k<NST;++k) printf("%7d", k);
  printf("   | %8s %8s\n", "rate/step", "n steps");
  for(int e=0;e<NEE;++e){
    printf("  %-13s", een[e]);
    std::vector<double> xs, ys;
    for(int k=0;k<NST;++k){
      double c=core(d[e][k]);
      if(c<0){ printf("%7s","-"); continue; }
      double st=qtl(sg[e][k],50), sh=qtl(hh[e][k],50);
      double ex=std::sqrt(st*st+sh*sh)/st;
      double v=(c/st)/ex;
      printf("%7.2f", v);
      if(v>0){ xs.push_back(k); ys.push_back(std::log(v)); }
    }
    if(xs.size()>=3){
      double sx=0,sy=0,sxx=0,sxy=0; int n=xs.size();
      for(int i=0;i<n;++i){ sx+=xs[i]; sy+=ys[i]; sxx+=xs[i]*xs[i]; sxy+=xs[i]*ys[i]; }
      double slope=(n*sxy-sx*sy)/(n*sxx-sx*sx);
      printf("   | %8.3f %8d", std::exp(slope), n);
    } else printf("   | %8s %8zu", "-", xs.size());
    printf("\n");
  }
  printf("\n  PREDICTION: rate/step RISES with |eta| if the high-|z| services drive it.\n");
  printf("  FLAT means they do not, and something else is compounding.\n\n");
}
