// v2tgeo.C(snitch, tgeo, offsets, nmax): map each entry of a G4Snitch file's volume table (tree V, path "name:copy/..."
// from the Geant4 world down) to a node of the TGeo export of the same geometry (Fireworks' TGeoMgrFromDdd), level by
// level: the daughter whose volume name, without the DDD namespace, and copy number match. Geant4's copy number is
// DDD's plus DDG4Builder's CopyNoOffset + CopyNoTag (SpecPars of the muon numbering); the table "offsets" (from
// copyno_offsets.py: "regex<TAB>offset" per line) gives them, and they are subtracted. The export has no CMSE level
// (an identity container): a level the TGeo file lacks is skipped when the next level is found under the same node.
// Prints the counts and the first failures; with ROOT's node ids (TGeoNodeCache::BuildIdArray()) a mapped path also
// has a physical node id, printed for the first few.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#include "TGeoCache.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
static std::string bare(const std::string &n) { auto p = n.find(':'); return p == std::string::npos ? n : n.substr(p + 1); }
void v2tgeo(const char *sfn, const char *gfn, const char *ofn, long nmax = -1) {
  std::vector<std::pair<std::regex, int>> offsets;
  {
    std::ifstream in(ofn);
    std::string re;
    int off;
    while (in >> re >> off) offsets.emplace_back(std::regex(re), off);
  }
  std::map<std::string, int> off_cache;
  auto offset_of = [&](const std::string &name) {
    auto it = off_cache.find(name);
    if (it != off_cache.end()) return it->second;
    int off = 0;
    for (auto &o : offsets) if (std::regex_match(name, o.first)) { off = o.second; break; }
    off_cache[name] = off;
    return off;
  };
  TFile gf(gfn);
  TGeoManager *gm = (TGeoManager *)gf.Get("cmsGeo");
  gm->GetCache()->BuildIdArray();
  TGeoNode *top = gm->GetTopNode();
  TFile sf(sfn);
  TTree *v = (TTree *)sf.Get("V");
  std::string *path = nullptr;
  double vx = 0, vy = 0, vz = 0;
  v->SetBranchAddress("path", &path);
  v->SetBranchAddress("x", &vx);
  v->SetBranchAddress("y", &vy);
  v->SetBranchAddress("z", &vz);
  double max_d = 0;
  long n_far = 0;
  long n = v->GetEntries(), ok = 0, with_offset = 0;
  if (nmax >= 0 && nmax < n) n = nmax;
  std::map<std::string, long> skipped, failed_at;
  int shown = 0;
  for (long e = 0; e < n; ++e) {
    v->GetEntry(e);
    std::vector<std::pair<std::string, int>> lv;
    std::stringstream ss(*path);
    std::string tok;
    while (std::getline(ss, tok, '/')) { auto c = tok.rfind(':'); lv.emplace_back(tok.substr(0, c), std::stoi(tok.substr(c + 1))); }
    for (auto &x : lv) { int o = offset_of(x.first); if (o) { x.second -= o; ++with_offset; } }
    auto find = [&](TGeoNode *node, size_t l) -> TGeoNode * {
      for (int i = 0; i < node->GetNdaughters(); ++i) {
        TGeoNode *d = node->GetDaughter(i);
        if (d->GetNumber() == lv[l].second && bare(d->GetVolume()->GetName()) == lv[l].first) return d;
      }
      return nullptr;
    };
    TGeoNode *node = top;
    std::string tpath = std::string("/") + top->GetName();
    bool good = true;
    for (size_t l = 1; l < lv.size(); ++l) {  // level 0 is the Geant4 world
      TGeoNode *next = find(node, l);
      if (!next && l + 1 < lv.size() && find(node, l + 1)) { ++skipped[lv[l].first]; continue; }
      if (!next) {
        if (failed_at[lv[l].first]++ < 2)
          printf("FAIL %s\n  at level %zu (%s, DDD copy %d) under %s\n", path->c_str(), l, lv[l].first.c_str(), lv[l].second, node->GetName());
        good = false;
        break;
      }
      node = next;
      tpath += std::string("/") + node->GetName();
    }
    if (good) {
      ++ok;
      // The volume's origin as Geant4 saw it against the mapped node's global position.
      gm->cd(tpath.c_str());
      const double *t = gm->GetCurrentMatrix()->GetTranslation();
      double d = std::max(std::max(std::abs(t[0] - vx), std::abs(t[1] - vy)), std::abs(t[2] - vz));
      if (d > max_d) max_d = d;
      if (d > 1e-3 && n_far++ < 3) printf("FAR %s: TGeo (%g, %g, %g), Geant4 (%g, %g, %g)\n", path->c_str(), t[0], t[1], t[2], vx, vy, vz);
      if (shown < 3 && lv.size() > 6) {
        gm->cd(tpath.c_str());
        printf("  %s\n   -> %s, node id %d\n", path->c_str(), tpath.c_str(), gm->GetCurrentNodeId());
        ++shown;
      }
    }
  }
  printf("%ld of %ld volumes mapped to a TGeo node (%ld path levels had a copy-number offset)\n", ok, n, with_offset);
  printf("  origin of the mapped node vs Geant4's: largest difference %.3g cm, %ld beyond 10 um\n", max_d, n_far);
  for (auto &x : skipped) printf("  level absent in TGeo, skipped: %s (%ld paths)\n", x.first.c_str(), x.second);
  for (auto &x : failed_at) printf("  failed at %s (%ld)\n", x.first.c_str(), x.second);
}
