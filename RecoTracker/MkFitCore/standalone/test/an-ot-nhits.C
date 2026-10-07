// Hits per OT layer on a track: how many overlaps the search picks up.
//
// For every track and every outer-tracker mkFit layer it has at least one hit
// in, count the hits it took in that layer. The +z and -z endcaps are folded
// together: mkFit layers 50-59 are counted as 28-37. Series:
//
//   sim        the sim tracks that v2p2 found, their own hits with a matching
//              mcTrackID (truth-bound, so a neighbour's hit does not count).
//              What the detector offered, counted in hits.
//   sim_mod    the same, counted in distinct MODULES per layer. Two hits of
//              one module (a split cluster) count once, because the in-layer
//              search vetoes a second hit in the same module. This is the
//              ceiling the search can reach.
//   v2p2       v2p2's final tracks (Event::candidateTracks_), all hits, with
//              the in-layer search (the production configuration).
//   v2p2_true  the same tracks, only hits whose mcTrackID is the track's label.
//   besthit    v2p2 with the in-layer search off: one best hit per layer.
//   besthit_true  the same, truth-matched hits only.
//   prod       the production tracks in the sample (Event::cmsswTracks_; needs
//              --read-cmssw-tracks and a sample written with them), all hits.
//   prod_true  the same, truth-matched hits only.
//
// A reco track's label is the majority mcTrackID over its valid hits, and the
// track is used only if that majority is at least half of them (the quality-val
// association rule, without the seed-hit exclusion). The seeds are pixel-only,
// so no OT hit is a seed hit.
//
// Output: one TH2F per series, OT layer (x) vs hits in that layer (y).
//
//   ot_nh_reset(); per event: ot_nh_ev(s.event(), mode); ot_nh_write("out.root");
//   mode 0: in-layer search on (fills sim, v2p2, prod); mode 1: off (besthit).

#include "RecoTracker/MkFitCMS/standalone/Shell.h"

#include "TFile.h"
#include "TH2F.h"

#include <algorithm>
#include <map>
#include <vector>

namespace {
  constexpr int kNSer = 8;
  constexpr int kMaxN = 8;  // counts above this go into the last bin
  const char *ot_ser_name[kNSer] = {"sim", "v2p2", "v2p2_true", "besthit", "besthit_true", "prod", "prod_true", "sim_mod"};

  // 22 OT layers: 4-15 (TBPS 4-9, TB2S 10-15), 28-37 (TEDD, both endcaps).
  int ot_slot(int lay) {
    if (lay >= 50 && lay <= 59) lay -= 22;
    if (lay >= 4 && lay <= 15) return lay - 4;
    if (lay >= 28 && lay <= 37) return 12 + (lay - 28);
    return -1;
  }
  int ot_layer_of_slot(int s) { return s < 12 ? s + 4 : s - 12 + 28; }

  long ot_cnt[kNSer][22][kMaxN + 1];
  long ot_n_tracks[kNSer];

  int hit_mc_track(const mkfit::Event *ev, int lay, int idx) {
    const int mc_hit = ev->layerHits_[lay][idx].mcHitID();
    if (mc_hit < 0 || mc_hit >= (int)ev->simHitsInfo_.size())
      return -1;
    return ev->simHitsInfo_[mc_hit].mcTrackID();
  }

  // Majority label of a reco track, or -1 if the majority is under half.
  int reco_label(const mkfit::Event *ev, const mkfit::Track &t) {
    std::map<int, int> votes;
    int n_valid = 0;
    for (int i = 0; i < t.nTotalHits(); ++i) {
      const int idx = t.getHitIdx(i), lay = t.getHitLyr(i);
      if (idx < 0 || lay < 0)
        continue;
      ++n_valid;
      const int mc = hit_mc_track(ev, lay, idx);
      if (mc >= 0)
        ++votes[mc];
    }
    int best = -1, best_n = 0;
    for (auto &[l, n] : votes)
      if (n > best_n) { best = l; best_n = n; }
    return (n_valid > 0 && 2 * best_n >= n_valid) ? best : -1;
  }

  // Fill one track: per OT layer, hits (all) and hits matching `label` (true).
  void fill_track(const mkfit::Event *ev, const mkfit::Track &t, int label, int ser_all, int ser_true,
                  int ser_mod = -1) {
    int n_all[22] = {0}, n_true[22] = {0};
    std::vector<std::vector<unsigned int>> mods(22);
    for (int i = 0; i < t.nTotalHits(); ++i) {
      const int idx = t.getHitIdx(i), lay = t.getHitLyr(i);
      if (idx < 0 || lay < 0)
        continue;
      const int s = ot_slot(lay);
      if (s < 0)
        continue;
      ++n_all[s];
      if (label >= 0 && hit_mc_track(ev, lay, idx) == label) {
        ++n_true[s];
        const unsigned int m = ev->layerHits_[lay][idx].detIDinLayer();
        if (std::find(mods[s].begin(), mods[s].end(), m) == mods[s].end())
          mods[s].push_back(m);
      }
    }
    for (int s = 0; s < 22; ++s) {
      if (ser_all >= 0 && n_all[s] > 0)
        ++ot_cnt[ser_all][s][std::min(n_all[s], kMaxN)];
      if (ser_true >= 0 && n_true[s] > 0)
        ++ot_cnt[ser_true][s][std::min(n_true[s], kMaxN)];
      if (ser_mod >= 0 && !mods[s].empty())
        ++ot_cnt[ser_mod][s][std::min((int)mods[s].size(), kMaxN)];
    }
    if (ser_mod >= 0) ++ot_n_tracks[ser_mod];
    if (ser_all >= 0) ++ot_n_tracks[ser_all];
    if (ser_true >= 0) ++ot_n_tracks[ser_true];
  }
}  // namespace

void ot_nh_reset() {
  for (int k = 0; k < kNSer; ++k) {
    ot_n_tracks[k] = 0;
    for (int s = 0; s < 22; ++s)
      for (int n = 0; n <= kMaxN; ++n)
        ot_cnt[k][s][n] = 0;
  }
}

void ot_nh_ev(const mkfit::Event *ev, int mode = 0) {
  if (!ev) return;
  if (mode == 1) {
    for (const auto &t : ev->candidateTracks_) {
      const int label = reco_label(ev, t);
      if (label >= 0) fill_track(ev, t, label, 3, 4);
    }
    return;
  }
  // v2p2, and the sim tracks it found (each sim track once per event).
  std::map<int, bool> found_sim;
  for (const auto &t : ev->candidateTracks_) {
    const int label = reco_label(ev, t);
    if (label < 0) continue;
    fill_track(ev, t, label, 1, 2);
    found_sim[label] = true;
  }
  for (auto &[label, _] : found_sim) {
    if (label >= (int)ev->simTracks_.size()) continue;
    fill_track(ev, ev->simTracks_[label], label, -1, 0, 7);
  }
  // production
  for (const auto &t : ev->cmsswTracks_) {
    const int label = reco_label(ev, t);
    if (label < 0) continue;
    fill_track(ev, t, label, 5, 6);
  }
}

void ot_nh_write(const char *fname) {
  TFile f(fname, "RECREATE");
  for (int k = 0; k < kNSer; ++k) {
    TH2F h(Form("h_%s", ot_ser_name[k]),
           Form("%s;OT layer;hits in the layer", ot_ser_name[k]),
           22, 0, 22, kMaxN, 0.5, kMaxN + 0.5);
    for (int s = 0; s < 22; ++s) {
      h.GetXaxis()->SetBinLabel(s + 1, Form("%d", ot_layer_of_slot(s)));
      for (int n = 1; n <= kMaxN; ++n)
        h.SetBinContent(s + 1, n, ot_cnt[k][s][n]);
    }
    h.Write();
  }
  f.Close();
  // Console summary: per layer group, fraction of layer crossings with >= 2 hits.
  const char *grp[4] = {"TBPS-P (4,6,8)", "TBPS-S (5,7,9)", "TB2S (10-15)", "TEDD (28-37, both z)"};
  auto in_grp = [](int g, int s) {
    const int l = ot_layer_of_slot(s);
    if (g == 0) return l >= 4 && l <= 9 && l % 2 == 0;
    if (g == 1) return l >= 4 && l <= 9 && l % 2 == 1;
    if (g == 2) return l >= 10 && l <= 15;
    return l >= 28 && l <= 37;
  };
  printf("\not_nh: tracks used: sim %ld, v2p2 %ld, besthit %ld, prod %ld\n", ot_n_tracks[0], ot_n_tracks[1],
         ot_n_tracks[3], ot_n_tracks[5]);
  printf("%-22s", "crossings with >=2 hits");
  for (int k = 0; k < kNSer; ++k) printf(" %11s", ot_ser_name[k]);
  printf("\n");
  for (int g = 0; g < 4; ++g) {
    printf("%-22s", grp[g]);
    for (int k = 0; k < kNSer; ++k) {
      long tot = 0, two = 0;
      for (int s = 0; s < 22; ++s) {
        if (!in_grp(g, s)) continue;
        for (int n = 1; n <= kMaxN; ++n) { tot += ot_cnt[k][s][n]; if (n >= 2) two += ot_cnt[k][s][n]; }
      }
      printf(" %10.2f%%", tot ? 100.0 * two / tot : 0.0);
    }
    printf("\n");
  }
}
