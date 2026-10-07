#include "RecoTracker/MkFitCore/standalone/DataFormats/ValStructs.h"
#include <map>
#include <set>
// IN-LAYER COVERAGE (not 'hermeticity' -- that is whole-detector solid-angle
// coverage, a calorimetry/MET property, and this is not it). Geometry-free. A sim track with disc hits at layers 16,17,19
// demonstrably CROSSED 18 and left nothing there. Counting those interpolation
// gaps measures "crossed but no hit" without needing to propagate anything.
//
// It is what the 82.9 % per-layer efficiency deliberately EXCLUDES: that
// denominator is layers where the sim track HAS a hit, so it measures the
// algorithm and says nothing about the detector.
//
// WHAT IT LUMPS TOGETHER and cannot separate: inter-module gaps in r and phi,
// dead modules, and modules crossed that did not register a hit. Call it
// active-area acceptance x module efficiency, not hermeticity.
//
// CAVEAT that bounds this from above: `has_sim_here` comes from
// countSimHitsInLayer (sim->rec, NOT gated by the bestTkIdx arbitration), so a
// discarded truth LABEL does not create a false gap here. Dead modules, module
// edges and genuine inter-module gaps all do.
void an_layercov(const char *fn="val-purity.root"){
  // (event, sim_label) -> layers searched, and which had a sim hit
  std::map<long long, std::pair<std::set<int>,std::set<int>>> T;  // searched, withhit
  TFile f(fn); TTree*t=(TTree*)f.Get("miss");
  ValSearchMiss*m=nullptr; t->SetBranchAddress("m",&m);
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(m->sim_label<0||m->layer<0) continue;
    auto &e=T[(long long)m->event*10000000LL+m->sim_label];
    e.first.insert(m->layer);
    if(m->has_sim_here && m->n_sim_in_layer>0) e.second.insert(m->layer); }

  struct Acc { long span=0, gap=0, trk=0; };
  Acc disk, pixb;
  for(auto&kv:T){
    const auto &srch=kv.second.first, &has=kv.second.second;
    // forward pixel discs 16-27, and pixel barrel 0-3, treated separately
    for(int pass=0;pass<2;++pass){
      const int lo = pass? 0:16, hi = pass? 3:27;
      int a=99,b=-1;
      for(int L : has) if(L>=lo&&L<=hi){ a=std::min(a,L); b=std::max(b,L); }
      if(b<0||b<=a) continue;
      long sp=0, gp=0;
      for(int L=a; L<=b; ++L){
        if(!srch.count(L)) continue;       // only layers the search visited
        ++sp; if(!has.count(L)) ++gp; }
      Acc &A = pass? pixb : disk;
      A.span+=sp; A.gap+=gp; ++A.trk; } }

  printf("\n===== IN-LAYER COVERAGE, from interpolation gaps =====\n");
  printf("Between a track's innermost and outermost hit in a region, how many of\n");
  printf("the layers it demonstrably crossed carry NO hit of that track?\n\n");
  printf("  %-16s %8s %10s %10s %9s\n","region","tracks","layers","empty","gap rate");
  printf("  %-16s %8ld %10ld %10ld %8.2f%%\n","fwd pix 16-27", disk.trk, disk.span,
         disk.gap, disk.span?100.0*disk.gap/disk.span:0.0);
  printf("  %-16s %8ld %10ld %10ld %8.2f%%\n","PixB 0-3", pixb.trk, pixb.span,
         pixb.gap, pixb.span?100.0*pixb.gap/pixb.span:0.0);
  printf("\n  This is the piece the 82.9%% per-layer number divides out. Combined:\n");
  if(disk.span) printf("    discs: 82.9%% of hit-bearing layers x (1 - %.4f) = %.1f%% of layers CROSSED\n",
         (double)disk.gap/disk.span, 82.9*(1.0-(double)disk.gap/disk.span));
  printf("\n");
}
