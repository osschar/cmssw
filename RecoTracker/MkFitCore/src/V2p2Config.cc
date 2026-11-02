// Standalone definitions of the knobs declared in V2p2Config.h. In the CMSSW
// build they are constexpr and this file defines nothing.
//
// Must be the first include, so the header is expanded with definitions.
#define MKFIT_V2P2_CONFIG_DEFINE
#include "RecoTracker/MkFitCore/src/V2p2Config.h"

#if defined(MKFIT_STANDALONE)

namespace mkfit::Config::V2p2 {

  void set_extra_dq(float f) {
    Window::dq_trk_fac = f;
    Window::dq_hit_fac = f * 1.2f;  // the old EXTRA_DQ ratio
  }

  V2p2Params v2p2_params_from_knobs() {
    V2p2Params p;
    p.score_outward = Score::fwd;
    p.score_inward = Score::bkw;
    p.hit_chi2_cut = Policy::hit_chi2_cut;
    p.max_presel_hits = InLayer::max_presel_hits;
    p.max_sec_depth = InLayer::max_sec_depth;
    return p;
  }

}  // namespace mkfit::Config::V2p2

#endif
