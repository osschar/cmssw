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
// - the RNTuples "Seeds" and "Tracks", one entry per seed after mkFit's
//   cleaning (TrCandMeta::seed indexes them) and per output track of mkFit,
//   written with the records. Fields: the five key fields, idx (the position in
//   its event's vector), label, state (mkfit::TrackState), chi2, score, status
//   (Track::Status as its 32 bits), n_found, hit_lyr and hit_idx (the hits on
//   the track, in order, as mkFit's HitOnTrack layer and index).
// - the truth of an event, written by writeTruth() when a truth module gives
//   it (TraceTruthTables below), keyed like the records:
//   "TruthHits", one entry per layer: layer, n_hits (1 + the largest index of a
//   hit in the layer), hit (the indices of the hits that have a sim track) and
//   sim (their sim track, an index into the event's SimTracks);
//   "SimTracks", one entry per sim track, with the fields of Seeds;
//   "TruthSeeds", one entry per input seed: seed (its index in the input seed
//   collection, which is the mkFit seed's label), sim and sim_frac (the best
//   sim track of the seed-to-TrackingParticle association and its quality),
//   best_sim and best_frac (the sim track with most of the seed's clusters, at
//   least 3, and the fraction of the clusters it has).
//
// "pass" counts how often the same (run, lumi, event, module) was written in
// this job, from 0: a looper re-running an event gives passes 0, 1, ... The
// truth tables count their passes on their own, in the order of the
// writeTruth() calls.
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

#include "RecoTracker/MkFitCore/interface/Track.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace mkfit {
  class TraceStore;
  struct TraceEventInfo;
}  // namespace mkfit

namespace mkfit_fwff {

  // The truth of one traced event, for TraceRNTupleWriter::writeTruth(); the RNTuples are described above.
  struct TraceTruthTables {
    unsigned int run = 0, lumi = 0;
    unsigned long long event = 0;
    std::string module;  // label of the MkFitProducer the truth is for

    struct LayerHits {
      int layer = -1;
      int n_hits = 0;        // 1 + the largest index of a hit in the layer
      std::vector<int> hit;  // indices of the hits that have a sim track
      std::vector<int> sim;  // their sim track, an index into simTracks
    };
    std::vector<LayerHits> layers;

    // The sim tracks; a track's label is its index.
    std::vector<mkfit::Track> simTracks;

    struct SeedSim {
      int seed = -1;  // index in the input seed collection, the mkFit seed's label
      int sim = -1;   // best sim track of the seed-to-TrackingParticle association, -1 none
      float sim_frac = 0;
      int best_sim = -1;  // the sim track with most of the seed's clusters (at least 3), -1 none
      float best_frac = 0;
    };
    std::vector<SeedSim> seeds;
  };

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

    // Appends the truth tables of one event, in the caller. Returns their pass. Thread safe.
    int writeTruth(const TraceTruthTables& truth);

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
