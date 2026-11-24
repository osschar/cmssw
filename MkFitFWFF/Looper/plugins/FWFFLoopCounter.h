#ifndef MkFitFWFF_Looper_FWFFLoopCounter_h
#define MkFitFWFF_Looper_FWFFLoopCounter_h

#include <atomic>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace fwff {
  // The loop number set by FWFFMkFitLooper at the start of each loop, read by the modules of the same plugin
  // library (MkFitOutputDigest) to label their output. It stays -1 in a job without the looper.
  inline std::atomic<int>& currentLoop() {
    static std::atomic<int> s_loop{-1};
    return s_loop;
  }

  // The MKFITDIGEST lines of the current loop. MkFitOutputDigest adds its line when the looper is present
  // (currentLoop() >= 0); the looper takes them at the end of the loop for its status file.
  struct DigestLines {
    std::mutex mutex;
    std::vector<std::string> lines;
  };
  inline DigestLines& digestLines() {
    static DigestLines s_lines;
    return s_lines;
  }
  inline void addDigestLine(std::string line) {
    auto& d = digestLines();
    std::lock_guard<std::mutex> lock(d.mutex);
    d.lines.push_back(std::move(line));
  }
  inline std::vector<std::string> takeDigestLines() {
    auto& d = digestLines();
    std::lock_guard<std::mutex> lock(d.mutex);
    std::vector<std::string> out;
    out.swap(d.lines);
    return out;
  }
}  // namespace fwff

#endif
