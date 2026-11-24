#ifndef RecoTracker_MkFitCore_interface_MkBuilderWrapper_h
#define RecoTracker_MkFitCore_interface_MkBuilderWrapper_h

#include <memory>

namespace mkfit {
  class MkBuilder;
#ifdef MKFIT_TRACE
  class TraceStore;
#endif

  /**
   * The purpose of this class is to hide the header of MkBuilder.h
   * from CMSSW. The headers included by MkBuilder.h contain uses of
   * the build-time configuration macros, that should remain as
   * internal details of MkFit package.
   */
  class MkBuilderWrapper {
  public:
    MkBuilderWrapper(bool silent);
    ~MkBuilderWrapper();

    MkBuilder& get() { return *builder_; }
#ifdef MKFIT_TRACE
    // The trace store of this builder's stream (MKFIT_TRACE builds only).
    TraceStore& trace() { return *trace_; }
#endif

    static void populate();
    static void clear();

  private:
    std::unique_ptr<MkBuilder> builder_;
#ifdef MKFIT_TRACE
    std::unique_ptr<TraceStore> trace_;
#endif
  };
}  // namespace mkfit

#endif
