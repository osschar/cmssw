#ifndef MkFitFWFF_Trace_TraceFileDigest_h
#define MkFitFWFF_Trace_TraceFileDigest_h

// Reads a trace file written by TraceRNTupleWriter back, record by record
// through the records' dictionaries, and computes for each traced event the
// digest of TraceDigest.h, as the writer service computes it from the store in
// memory. Equal digest lines from the file and from memory show that the file
// holds the store exactly. Also counts a few things for sanity checks.

#include "MkFitFWFF/Trace/interface/TraceDigest.h"

#include <array>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace mkfit_fwff {

  struct TraceFileEvent {
    // From the Events RNTuple.
    std::uint32_t run = 0, lumi = 0;
    std::uint64_t event = 0;
    std::string module;
    int pass = 0;
    int n_seeds_in = 0, n_seeds = 0, n_tracks = 0;
    std::array<int, 7> n_rec{};  // record counts as written into Events
    std::vector<int> out_labels;

    // From the record RNTuples.
    std::array<trace_digest::VecDigest, 7> digest{};
    int n_exported = 0;       // TrCandMeta with cand >= 0
    int n_labeled = 0;        // TrCandMeta with global_seed >= 0
    int n_out_matched = 0;    // out_labels equal to the global_seed of an exported TrCandMeta
    int n_cand_distinct = 0;  // distinct cand values of the exported TrCandMeta
    int n_stage[3]{};         // TrCandStage by stage: 0 FwdSearch, 1 BkwFit, 2 BkwSearch
    int n_stage_final[3]{};   // of those, with final_state_id >= 0
    int n_on_final_path = 0;  // TrCandState with on_final_path
  };

  // Throws std::runtime_error on a missing RNTuple or a record whose key is not in Events.
  std::vector<TraceFileEvent> readTraceFile(const std::string& fileName);

  // One line per record vector, the same from the file and from memory:
  //   MKFITTRACE digest <run>:<lumi>:<event> <module> pass=<p> <vector> <size> <hash>
  std::string digestLine(std::uint32_t run,
                         std::uint32_t lumi,
                         std::uint64_t event,
                         std::string_view module,
                         int pass,
                         const trace_digest::VecDigest& d);

}  // namespace mkfit_fwff

#endif
