#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
#include <map>
#include <vector>
// 33-56% of divergences take a wrong hit in a layer where the sim track has NO
// hit. Is that GEOMETRIC -- the track cannot reach that layer at its eta?
// Reachability of a barrel layer is eta_max = asinh(z_half / R), from
// matlayers.txt (r_centroid, z_half): IT1 3.06/19.96, IT2 6.18/20.12,
// IT3 10.52/20.12, IT4 14.71/20.12, OT1 23.57/120.76, OT2 36.32/120.90,
// OT3 51.45/120.87. mkFit splits each OT layer in two, so 4,5->OT1 etc.
static double q(std::vector<double> u,double p){ if(u.empty())return 0;
  std::sort(u.begin(),u.end()); return u[(size_t)(p/100.*(u.size()-1))]; }
struct Hit { int step,layer; bool mc,acc; float chi2,eta; };
void an_nohit(const char *fn="val-purity.root"){
  // eta_max per mkFit barrel layer
  double emax[16];
  const double Rc[7]={3.06,6.18,10.52,14.71,23.57,36.32,51.45};
  const double Zh[7]={19.96,20.12,20.12,20.12,120.76,120.90,120.87};
  for(int L=0;L<4;++L) emax[L]=asinh(Zh[L]/Rc[L]);
  for(int L=4;L<10;++L){ int o=4+(L-4)/2; emax[L]=asinh(Zh[o]/Rc[o]); }
  for(int L=10;L<16;++L) emax[L]=asinh(117.68/(L<12?68.8:(L<14?86.07:108.36)));

  std::map<long long,std::vector<Hit>> C;
  TFile f(fn); TTree*t=(TTree*)f.Get("search");
  ValSearchHit*h=nullptr; t->SetBranchAddress("h",&h);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!h->had_kalman||h->chi2<-900||h->global_seed<0||h->eta<-900) continue;
    C[(long long)h->event*1000000LL+h->global_seed].push_back(
      {h->step,h->layer,h->mc_match,h->accepted,h->chi2,h->eta}); }

  long n=0, beyond=0, inside=0; std::vector<double> mE, mR;
  long nb[4]={},bb[4]={};                        // per pixel layer
  for(auto&kv:C){
    auto v=kv.second; std::sort(v.begin(),v.end(),
      [](const Hit&a,const Hit&b){return a.step<b.step;});
    int added=0; for(auto&x:v) if(x.acc&&x.chi2<30.f) ++added;
    if(added<3) continue;
    for(size_t i=0;i<v.size();++i){
      if(!(v[i].acc && v[i].chi2<30.f && !v[i].mc)) continue;
      bool has=false;
      for(auto&x:v) if(x.step==v[i].step&&x.layer==v[i].layer&&x.mc) has=true;
      if(!has && v[i].layer>=0 && v[i].layer<16){
        ++n; const double ae=std::fabs(v[i].eta), em=emax[v[i].layer];
        mE.push_back(ae); mR.push_back(ae/em);
        if(ae>em) ++beyond; else ++inside;
        if(v[i].layer<4){ ++nb[v[i].layer]; if(ae>em) ++bb[v[i].layer]; } }
      break; } }
  printf("\n===== IS \"NO MC HIT IN THAT LAYER\" GEOMETRIC? =====\n");
  printf("Barrel layers only. eta_max = asinh(z_half/R) is where a track from the\n");
  printf("origin stops crossing that layer at all.\n\n");
  printf("  no-MC-hit divergences in a barrel layer : %ld\n", n);
  printf("  of those, |eta| BEYOND that layer's reach: %ld = %.1f%%\n",
         beyond, n?100.0*beyond/n:0.0);
  printf("  |eta|/eta_max  p25 %.2f  median %.2f  p75 %.2f  p90 %.2f\n",
         q(mR,25), q(mR,50), q(mR,75), q(mR,90));
  printf("  |eta|          p25 %.2f  median %.2f  p75 %.2f\n", q(mE,25), q(mE,50), q(mE,75));
  printf("\n  per pixel-barrel layer:\n");
  for(int L=0;L<4;++L) if(nb[L]>=20)
    printf("    layer %d (eta_max %.2f) : %5ld cases, %5.1f%% beyond reach\n",
           L, emax[L], nb[L], 100.0*bb[L]/nb[L]);
  printf("\n  >50%% beyond reach = the search is being sent into layers the track\n");
  printf("  cannot cross, i.e. the layer plan (S14), not the hit selection.\n\n");
}
