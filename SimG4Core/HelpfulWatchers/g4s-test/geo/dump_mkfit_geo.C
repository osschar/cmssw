// dump_mkfit_geo.C(bin, out): mkFit's TrackerInfo from a geometry binary into a ROOT file: tree "M" (one entry per
// module: layer, is_barrel, is_stereo, detid, pos, zdir, xdir, radl, bbxi), tree "L" (per layer: id, barrel, stereo,
// rin, rout, zmin, zmax) and the (|z|, r) material grid as TH2F "radl" and "bbxi".
// Compiles mkFit's TrackerInfo.cc in, from a CMSSW src tree with MkFitCore (format-5 reader, #190 and later), with
// any ROOT (here the system one on black):
//   root -l -b -q -e 'gSystem->AddIncludePath("-I/foo/matevz/mic-dev/cmssw_20_v2p2/src -DMKFIT_STANDALONE");' \
//        -e '.x dump_mkfit_geo.C+("/foo/matevz/mic-dev/CMS-phase2.bin-v8-Run4D121", "mkfit_geo_D121.root")'
// D121's dump serves D128: the tracker is T35 in both.
#include "RecoTracker/MkFitCore/src/TrackerInfo.cc"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"

void dump_mkfit_geo(const char* bin, const char* out) {
  mkfit::TrackerInfo ti;
  ti.read_bin_file(bin);
  TFile f(out, "RECREATE");
  TTree m("M", "mkFit modules");
  int lay, brl, ste; unsigned int detid; float px, py, pz, nx, ny, nz, xx, xy, xz, radl, bbxi;
  m.Branch("layer", &lay); m.Branch("barrel", &brl); m.Branch("stereo", &ste); m.Branch("detid", &detid);
  m.Branch("px", &px); m.Branch("py", &py); m.Branch("pz", &pz);
  m.Branch("nx", &nx); m.Branch("ny", &ny); m.Branch("nz", &nz);
  m.Branch("xx", &xx); m.Branch("xy", &xy); m.Branch("xz", &xz);
  m.Branch("radl", &radl); m.Branch("bbxi", &bbxi);
  TTree L("L", "mkFit layers");
  float rin, rout, zmin, zmax;
  L.Branch("layer", &lay); L.Branch("barrel", &brl); L.Branch("stereo", &ste);
  L.Branch("rin", &rin); L.Branch("rout", &rout); L.Branch("zmin", &zmin); L.Branch("zmax", &zmax);
  for (int l = 0; l < ti.n_layers(); ++l) {
    const auto& li = ti.layer(l);
    lay = l; brl = li.is_barrel(); ste = li.is_stereo();
    rin = li.rin(); rout = li.rout(); zmin = li.zmin(); zmax = li.zmax();
    L.Fill();
    for (int i = 0; i < li.n_modules(); ++i) {
      const auto& mi = li.module_info(i);
      detid = mi.detid;
      px = mi.pos[0]; py = mi.pos[1]; pz = mi.pos[2];
      nx = mi.zdir[0]; ny = mi.zdir[1]; nz = mi.zdir[2];
      xx = mi.xdir[0]; xy = mi.xdir[1]; xz = mi.xdir[2];
      radl = mi.radl; bbxi = mi.bbxi;
      m.Fill();
    }
  }
  const int nzb = ti.mat_nbins_z(), nrb = ti.mat_nbins_r();
  TH2F hr("radl", "mkFit material grid radl;|z| [cm];r [cm]", nzb, 0, ti.mat_range_z(), nrb, 0, ti.mat_range_r());
  TH2F hx("bbxi", "mkFit material grid bbxi;|z| [cm];r [cm]", nzb, 0, ti.mat_range_z(), nrb, 0, ti.mat_range_r());
  for (int i = 0; i < nzb; ++i)
    for (int j = 0; j < nrb; ++j) {
      hr.SetBinContent(i + 1, j + 1, ti.material_radl(i, j));
      hx.SetBinContent(i + 1, j + 1, ti.material_bbxi(i, j));
    }
  f.Write();
  printf("%d layers, %lld modules, grid %d x %d over |z| < %.0f, r < %.0f\n", ti.n_layers(), m.GetEntries(), nzb, nrb,
         ti.mat_range_z(), ti.mat_range_r());
}
