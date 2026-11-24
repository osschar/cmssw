// The trace sink of TraceStore.h. Empty without MKFIT_TRACE.

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

}  // namespace mkfit

#endif
