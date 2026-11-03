#ifndef RecoTracker_MkFitCore_interface_V2p2Params_h
#define RecoTracker_MkFitCore_interface_V2p2Params_h

namespace mkfit {

  // Coefficients of the layer-step score (src/V2p2Score.h), one set per search
  // direction.
  struct V2p2ScoreParams {
    float hit_bonus = 3.0f;      // per hit taken; a hole wins above chi2 11, see the design notes
    float overlap_bonus = 0.0f;  // per hit beyond the first in a layer
    float chi2_weight = 1.0f;    // per unit of chi2
    float miss_penalty = 8.0f;   // a real hole
    float edge_penalty = 0.0f;   // crossing at the layer boundary
    float gap_penalty = 0.0f;    // crossing in an inactive region
    float stop_penalty = 8.0f;   // out of hole budget
    float hit_eff = 0.90f;       // per-layer hit efficiency, likelihood mode only
    // Per layer group, likelihood mode only; a negative value takes hit_eff. Groups, see
    // v2p2_hit_eff_group(): 0 pixel barrel, 1 pixel endcap, 2 strip barrel below layer 10
    // (phase-2 TBPS), 3 strip barrel from layer 10 (TB2S), 4 strip endcap (TEDD).
    static constexpr int n_hit_eff_grp = 5;
    float hit_eff_grp[n_hit_eff_grp] = {-1.0f, -1.0f, -1.0f, -1.0f, -1.0f};
  };

  // The MkFinderV2p2 tunings that can be set per iteration. In CMSSW they are
  // IterationConfig::m_v2p2, set by MkFitIterationConfigESProducer from its python
  // PSet "v2p2" and not saved to or loaded from JSON; of the score sets the PSet
  // carries hit_eff and hit_eff_grp, the fields the log-likelihood score reads. In
  // the standalone the mkFit options and the shell setters change them in
  // Config::ItrInfo. The values below are the defaults of both.
  struct V2p2Params {
    // Outward: a phase-2 tuning (CMSSW MTV, initialStep with the offline seeds, cleaner and
    // refit, 2026-10): eps 0.10, TBPS 0.14, TB2S 0.20. It will move with the seeds and the
    // final fit.
    V2p2ScoreParams score_outward{.hit_eff = 0.10f, .hit_eff_grp = {-1.0f, -1.0f, 0.14f, 0.20f, -1.0f}};
    V2p2ScoreParams score_inward{};

    // Kalman hit acceptance, chi2_acc < hit_chi2_cut, see Config::V2p2::Policy. v2p2 does
    // not use IterationParams::chi2Cut_min.
    float hit_chi2_cut = 9.21f;

    int max_presel_hits = 6;  // reduction cap per sub-layer
    int max_sec_depth = 4;    // most hits one path may take in a layer

    static constexpr int max_presel_hits_limit = 12;  // inline storage
    static constexpr int max_sec_depth_limit = 8;     // array size

    const V2p2ScoreParams& score(bool is_outward) const { return is_outward ? score_outward : score_inward; }
  };

}  // namespace mkfit

#endif
