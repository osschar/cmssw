#include "RecoTracker/MkFitCore/interface/MkBuilderWrapper.h"
#include "RecoTracker/MkFitCore/interface/MkBuilder.h"
#ifdef MKFIT_TRACE
#include "RecoTracker/MkFitCore/interface/TraceStore.h"
#endif

namespace mkfit {
  MkBuilderWrapper::MkBuilderWrapper(bool silent) : builder_(MkBuilder::make_builder(silent)) {
#ifdef MKFIT_TRACE
    trace_ = std::make_unique<TraceStore>();
#endif
  }
  MkBuilderWrapper::~MkBuilderWrapper() = default;

  void MkBuilderWrapper::populate() { MkBuilder::populate(); }
  void MkBuilderWrapper::clear() { MkBuilder::clear(); }

}  // namespace mkfit
