#ifndef RecoTracker_MkFitCore_standalone_DataFormats_RntConversions_h
#define RecoTracker_MkFitCore_standalone_DataFormats_RntConversions_h

#include "RecoTracker/MkFitCore/standalone/DataFormats/RntStructs.h"

#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/src/Matrix.h"
#include "RecoTracker/MkFitCore/src/MiniPropagators.h"

namespace mkfit {
  namespace miprops = mkfit::mini_propagators;

  inline EVec3 state2pos(const miprops::State &s) { return {s.x, s.y, s.z}; }
  inline EVec3 state2mom(const miprops::State &s) { return {s.px, s.py, s.pz}; }
  inline EBiVec3 state2bivec3(const miprops::State &s) { return {state2pos(s), state2mom(s)}; }

  inline EVec3 statep2pos(const miprops::StatePlex &s, int i) { return {s.x[i], s.y[i], s.z[i]}; }
  inline EVec3 statep2mom(const miprops::StatePlex &s, int i) { return {s.px[i], s.py[i], s.pz[i]}; }
  inline EBiVec3 statep2bivec3(const miprops::StatePlex &s, int i) { return {statep2pos(s, i), statep2mom(s, i)}; }
  inline PropInfo statep2propinfo(const miprops::StatePlex &s, int i) {
    return {statep2bivec3(s, i), s.dalpha[i], s.fail_flag[i]};
  }

  inline EVec3 hit2pos(const Hit &h) { return {h.x(), h.y(), h.z()}; }
  inline EVec3 track2pos(const TrackBase &s) { return {s.x(), s.y(), s.z()}; }
  inline EVec3 track2mom(const TrackBase &s) { return {s.px(), s.py(), s.pz()}; }
  inline EBiVec3 track2bivec3(const TrackBase &s) { return {track2pos(s), track2mom(s)}; }

  inline SimSeedInfo evsi2ssinfo(const Event *ev, int seed_idx) {
    SimSeedInfo ssi;
    Event::SimInfoFromHits sifh = ev->simInfoForCurrentSeed(seed_idx);
    if (sifh.is_set()) {
      ssi.s_sim = track2bivec3(ev->simTracks_[sifh.label]);
      ssi.sim_lbl = sifh.label;
      ssi.n_hits = sifh.n_hits;
      ssi.n_match = sifh.n_match;
      ssi.has_sim = true;
    }
    auto seed = ev->currentSeed(seed_idx);
    ssi.s_seed = track2bivec3(seed);
    ssi.seed_lbl = seed.label();
    ssi.seed_idx = seed_idx;
    return ssi;
  }
}  // namespace mkfit

#endif
