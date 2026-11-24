#ifndef MkFitFWFF_Trace_TraceRNTupleWriter_h
#define MkFitFWFF_Trace_TraceRNTupleWriter_h

// Writes the mkFit V2p2 trace (RecoTracker/MkFitCore/interface/TraceStore.h)
// into one ROOT file, as plain RNTuples:
//
// - one RNTuple per record type, named as the record (TrCandMeta, TrCandStage,
//   TrCandState, TrLayerSearch, TrHitMatch, TrKalmanUpdate, TrBkFitUpdate), one
//   entry per record, in store order. Fields: run, lumi, event, module, pass,
//   and the record itself as the class field "rec". A record's id is its index
//   in its event's vector, so the ids link records within one
//   (run, lumi, event, module, pass).
// - the RNTuple "Events", one entry per traced event: the same five key fields,
//   n_seeds_in, n_seeds (after cleaning), n_tracks (mkFit's output), the number
//   of records of each type (n_TrCandMeta ...) and out_labels, the labels of the
//   output tracks in order. An event without records still has its entry.
//
// "pass" counts how often the same (run, lumi, event, module) was written in
// this job, from 0: a looper re-running an event gives passes 0, 1, ...
//
// ROOT's implicit multi-threading is off for every writer. write() is thread
// safe (one mutex); the file is written from whichever stream calls it.

#include <memory>
#include <mutex>
#include <string>

namespace mkfit {
  class TraceStore;
  struct TraceEventInfo;
}  // namespace mkfit

namespace mkfit_fwff {

  class TraceRNTupleWriter {
  public:
    static constexpr const char* kEventsName = "Events";

    // Creates (recreates) the file. Throws std::runtime_error if it cannot.
    explicit TraceRNTupleWriter(const std::string& fileName);
    ~TraceRNTupleWriter();

    TraceRNTupleWriter(const TraceRNTupleWriter&) = delete;
    TraceRNTupleWriter& operator=(const TraceRNTupleWriter&) = delete;

    // Appends one event's trace. Returns its pass. Thread safe.
    int write(const mkfit::TraceStore& store, const mkfit::TraceEventInfo& info);

    // Commits every RNTuple and closes the file; write() afterwards throws.
    // Called by the destructor if not called before.
    void close();

    long long nEventsWritten() const { return m_n_events; }

  private:
    struct Impl;
    std::mutex m_mutex;
    std::unique_ptr<Impl> m_impl;
    long long m_n_events = 0;
  };

}  // namespace mkfit_fwff

#endif
