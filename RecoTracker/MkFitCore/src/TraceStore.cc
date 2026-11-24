// The trace sink of TraceStore.h and apply_trace_truth(). Empty without MKFIT_TRACE.

#ifdef MKFIT_TRACE

#include "RecoTracker/MkFitCore/interface/TraceStore.h"

#include "RecoTracker/MkFitCore/interface/cms_common_macros.h"

namespace mkfit {

  namespace {
    // Set before the first event and cleared after the last (TraceStore.h).
    CMS_SA_ALLOW TraceSink g_trace_sink;
  }  // namespace

  void set_trace_sink(TraceSink sink) { g_trace_sink = std::move(sink); }

  const TraceSink &trace_sink() { return g_trace_sink; }

  void apply_trace_truth(TraceStore &s, const TraceTruthSource &t) {
    const int event = t.tr_event_id();
    for (TrCandMeta &m : s.trCandMetas_) {
      m.event = event;
      if (m.cand >= 0)
        m.global_seed = t.tr_seed_label(m.seed);
    }
    for (TrHitMatch &h : s.trHitMatches_) {
      const int seed = s.trCandMetas_[s.trCandStates_[h.state_id].meta_id].seed;
      h.mc_match = t.tr_mc_match(seed, t.tr_mc_hit_id(h.layer, h.hit));
    }
    for (TrBkFitUpdate &u : s.trBkFitUpdates_)
      u.mc_track_id = t.tr_mc_track_id(t.tr_mc_hit_id(u.layer, u.hit));
  }

}  // namespace mkfit

#endif
