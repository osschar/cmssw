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
// Writing. With Options::threads > 0, write() takes the store's seven record
// vectors by swapping them with empty ones (the vectors are mutable members),
// queues the event and returns. The store is left empty, with the capacity of
// an earlier event's vectors. Writer threads, each owning some of the RNTuples,
// fill the RNTuples in the order of the write() calls. At most
// Options::maxEventsInFlight events are held; write() waits for one to finish
// beyond that. The records are bound to the RNTuple entries in place, not
// copied. With threads == 0, write() fills every RNTuple itself and leaves the
// store as it is. The file is the same either way.
//
// ROOT's implicit multi-threading is off for every RNTuple. write() is thread
// safe. Writing pages into the file is serialised by one mutex. An exception in
// a writer thread is rethrown by the next write() and by close().

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <string_view>

namespace mkfit {
  class TraceStore;
  struct TraceEventInfo;
}  // namespace mkfit

namespace mkfit_fwff {

  class TraceRNTupleWriter {
  public:
    static constexpr const char* kEventsName = "Events";

    struct Options {
      // ROOT's compression setting, 100 * algorithm + level (Compression.h); 501 is zstd at level 1.
      int compression = 501;
      // The size, in bytes, up to which an uncompressed page grows.
      std::size_t maxUnzippedPageSize = 4 * 1024 * 1024;
      // Writer threads; 0 writes in the caller of write().
      int threads = 4;
      // Events taken from stores and not yet written, at most (threads > 0).
      int maxEventsInFlight = 2;
    };

    // One event, once its records are in every RNTuple.
    struct EventDone {
      unsigned int run = 0, lumi = 0;
      unsigned long long event = 0;
      std::string_view module;
      int pass = 0;
      std::size_t n_records = 0;
      double take_ms = 0;   // time spent in write()
      double write_ms = 0;  // time spent filling this event into the RNTuples, summed over them
      double done_ms = 0;   // from the call of write() until the event is in every RNTuple
    };
    // Called once per event with the event's records, from a writer thread when threads > 0 (one call at a time),
    // else from write(). An exception it throws is treated as a writing error.
    using EventDoneFn = std::function<void(const mkfit::TraceStore& records, const EventDone& done)>;

    // Creates (recreates) the file and starts the writer threads. Throws std::runtime_error if it cannot.
    TraceRNTupleWriter(const std::string& fileName, const Options& options, EventDoneFn onDone = nullptr);
    explicit TraceRNTupleWriter(const std::string& fileName) : TraceRNTupleWriter(fileName, Options()) {}
    ~TraceRNTupleWriter();

    TraceRNTupleWriter(const TraceRNTupleWriter&) = delete;
    TraceRNTupleWriter& operator=(const TraceRNTupleWriter&) = delete;

    // Appends one event's trace. Returns its pass. Thread safe. With threads > 0 it takes the store's records
    // (above); info is used only during the call.
    int write(const mkfit::TraceStore& store, const mkfit::TraceEventInfo& info);

    // Writes the queued events, commits every RNTuple and closes the file; write() afterwards throws. Called by the
    // destructor if not called before.
    void close();

    long long nEventsWritten() const;

  private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
  };

}  // namespace mkfit_fwff

#endif
