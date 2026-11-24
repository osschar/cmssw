#ifndef MkFitFWFF_Looper_FWFFLoopCounter_h
#define MkFitFWFF_Looper_FWFFLoopCounter_h

#include <atomic>

namespace fwff {
  // The loop number set by FWFFMkFitLooper at the start of each loop, read by the modules of the same plugin
  // library (MkFitOutputDigest) to label their output. It stays -1 in a job without the looper.
  inline std::atomic<int>& currentLoop() {
    static std::atomic<int> s_loop{-1};
    return s_loop;
  }
}  // namespace fwff

#endif
