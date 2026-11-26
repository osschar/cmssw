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
//   These three are the truth source "tp" (TraceTruthTables::source empty or
//   "tp"; cmsRun: MkFitTraceTruth, TrackingParticles and the cluster charge).
// - the truth of another source (TraceTruthTables::source, a name of letters
//   and digits other than "tp"; cmsRun: MkFitTraceTruthNtuple writes
//   "ntuple"), keyed like the records, with passes of its own per source:
//   "TruthHits_<source>" and "SimTracks_<source>" with the fields of TruthHits
//   and SimTracks; "TruthSeeds_<source>" with the fields of TruthSeeds and
//   hit_lyr, hit_idx (the input seed's hits, as mkFit's layer and the
//   cluster's index, so that a seed can be matched by its hits); and
//   "SimHitStates_<source>", one entry per layer: layer, n_hits, and per hit
//   that has a sim hit chosen by the source, in index order: hit, sh_idx (the
//   source's index of the sim hit), x, y, z, px, py, pz (the sim hit's global
//   position and momentum at entry: the standalone's Event::simHitStates_),
//   ts_q and ts_par (charge and CCS parameters of the state at the sim hit:
//   the standalone's Event::simTrackStates_, with no covariance), linked (1
//   when a digi-sim link on the hit's cluster names the sim hit under the
//   decoding of the links' sim hit index that the event's links agree with, 0
//   when not, -1 when the job did not check). A source's RNTuples are created
//   with its first event, so a file without the source has none of them.
//   A source whose sim tracks are not TrackingParticles (cmsRun:
//   MkFitTraceTruthSim writes "sim", whose sim tracks are SimTracks) can also
//   write "SimTrackMap_<source>", one entry per event: per sim track, in the
//   order of SimTracks_<source>, g4_id (SimTrack::trackId()), bx and ev (its
//   EncodedEventId: bunch crossing, and event number in it; 0 and 0 for the
//   signal), pdg (SimTrack::type()), proc (the creating process, the
//   processType() of its SimVertex; -1 when the vertex is not found),
//   parent_g4 (that SimVertex's parentIndex(), the parent's g4_id; -1 none),
//   gen (SimTrack::genpartIndex()), n_own_simhits (its own sim hits in the
//   tracker), and one field tp_<name> per other source named by the writing
//   module (the index of the TrackingParticle that lists the SimTrack among
//   its g4Tracks, in that source's sim tracks; -1 none). It is created with
//   the source's first event when that event's tables carry the map.
//
// - mkFit's inputs of an event, written by writeInputs() when a module gives
//   them (TraceInputTables below; cmsRun: MkFitTraceInputs), keyed like the
//   records. They need no simulation.
//   "Hits", one entry per layer: layer, n_hits (1 + the largest index of a hit
//   in the layer), and per hit of the layer, in index order: idx (the index the
//   records use, which is the cluster's index in its cluster collection), x, y,
//   z and err (the hit as mkFit has it: Hit::position() and Hit::errArray(),
//   the symmetric 3x3 covariance as xx, xy, yy, xz, yz, zz), detid (the CMSSW
//   DetId), det_in_layer (Hit::detIDinLayer()), charge_pcm
//   (Hit::chargePerCM()), span_rows and span_cols (Hit::spanRows(),
//   spanCols()).
//   "Inputs", one entry per event: the beam spot of mkFit's EventOfHits
//   (bs_x, bs_y, bs_z, bs_sigma_z, bs_width_x, bs_width_y, bs_dxdz, bs_dydz),
//   n_hits_pixel and n_hits_outer (the sizes of the two hit vectors, i.e. the
//   index spaces of the pixel and the outer-tracker layers), the dead bins of
//   the EventOfHits (n_phi_bins, and per dead bin dead_layer, dead_q_bin,
//   dead_phi_bin), and the identity of mkFit's TrackerInfo: geom_md5 (the md5
//   of TrackerInfo::write_bin_file() with geom_version as its stamp) and
//   geom_version.
//   "Geometry", one entry per distinct geom_md5 in the file, written with the
//   first event that has it: md5, version, and bin, the bytes of that file.
//
// - the origin of each hit's charge, written by writeHitOrigin() when a module
//   gives it (TraceHitOriginTables below; cmsRun: MkFitTraceHitOrigin), keyed
//   like the records. It needs the simulation's sim hits, sim tracks and
//   digi-sim links, not the TrackingParticles.
//   "TruthHitOrigin", one entry per layer: layer, n_hits (1 + the largest
//   index of a hit in the layer), and per hit with simulated charge (at least
//   one digi-sim link on its cluster), in index order: hit (its index, as the
//   records use it), sim (the sim track, an index into SimTracks, of the
//   signal SimTrack with most of the hit's charge; -1 when none or when the
//   job had no TrackingParticles), f_own, f_sec, f_pu and f_unm (fractions of
//   the hit's charge: ionisation by a SimTrack's own particle; by a secondary
//   that is not a SimTrack but carries a SimTrack's id or the id of no
//   SimTrack; by pileup, whose sim hits the job does not have; signal charge
//   whose sim hit was not found; the rest is charge without a link, noise),
//   cls (the largest of own, secondary and unknown = pileup + unmatched:
//   1, 2, 3; 0 when all three are 0), sec_pdg and sec_proc (the particle type
//   and the creating process of the secondary sim hit with most charge in the
//   hit, 0 and 0 when none). These fields describe the signal's charge only,
//   whether or not the job had the pileup's sim hits.
//   Then the pileup's part, appended after them: f_pu_own, f_pu_sec and
//   f_pu_unm (the pileup's charge, f_pu, split by the same rule with the
//   pileup event's own SimTracks: own, secondary, and pileup charge whose sim
//   hit was not found or not given), cls_all (the largest of own, secondary
//   and unknown over all the charge: f_own + f_pu_own, f_sec + f_pu_sec,
//   f_unm + f_pu_unm; values as cls), sim_all (as sim, over the signal and the
//   pileup SimTracks with classified charge) and sec_pdg_all, sec_proc_all (as
//   sec_pdg, sec_proc, over the signal and the pileup sim hits). Without the
//   pileup's sim hits f_pu_unm equals f_pu, f_pu_own and f_pu_sec are 0, and
//   cls_all, sim_all, sec_pdg_all, sec_proc_all equal cls, sim, sec_pdg,
//   sec_proc.
//
// "pass" counts how often the same (run, lumi, event, module) was written in
// this job, from 0: a looper re-running an event gives passes 0, 1, ... The
// truth tables count their passes on their own, in the order of the
// writeTruth() calls and per truth source, and so do the input tables (writeInputs()) and the hit
// origin (writeHitOrigin()).
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

#include <array>
#include <cstddef>
#include <cstdint>
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
    // The truth source: empty or "tp" writes TruthHits, SimTracks and TruthSeeds; another name (letters and digits)
    // writes TruthHits_<source>, SimTracks_<source>, TruthSeeds_<source> and SimHitStates_<source>.
    std::string source;

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
      // The input seed's hits (mkFit's layer, the cluster's index); written for a named source only.
      std::vector<std::int32_t> hit_lyr, hit_idx;
    };
    std::vector<SeedSim> seeds;

    // The sim hit a named source chose for each hit (SimHitStates_<source>); not written for "tp".
    struct LayerSimHitStates {
      int layer = -1;
      int n_hits = 0;                    // 1 + the largest index of a hit in the layer
      std::vector<std::int32_t> hit;     // indices of the hits with a chosen sim hit
      std::vector<std::int32_t> sh_idx;  // the source's index of that sim hit
      std::vector<float> x, y, z, px, py, pz;
      std::vector<std::int8_t> ts_q;
      std::vector<std::array<float, 6>> ts_par;
      std::vector<std::int8_t> linked;
    };
    std::vector<LayerSimHitStates> simHitStates;

    // Each sim track's SimTrack and its TrackingParticles in other sources (SimTrackMap_<source>), for a named source
    // whose sim tracks are SimTracks; written when present is set in the source's first event, and then required in
    // every event of the source, with the same tp_names.
    struct SimTrackMap {
      bool present = false;
      std::vector<std::uint32_t> g4_id;
      std::vector<std::int32_t> bx, ev, pdg, proc, parent_g4, gen, n_own_simhits;
      std::vector<std::string> tp_names;                // the other sources, one field tp_<name> each
      std::vector<std::vector<std::int32_t>> tp_index;  // per name, per sim track: its TrackingParticle there, -1 none
    };
    SimTrackMap simTrackMap;
  };

  // mkFit's inputs of one traced event, for TraceRNTupleWriter::writeInputs(); the RNTuples Hits, Inputs and Geometry
  // are described above.
  struct TraceInputTables {
    unsigned int run = 0, lumi = 0;
    unsigned long long event = 0;
    std::string module;  // label of the MkFitProducer the inputs are for

    struct LayerHits {
      int layer = -1;
      int n_hits = 0;  // 1 + the largest index of a hit in the layer
      std::vector<std::int32_t> idx;
      std::vector<float> x, y, z;
      std::vector<std::array<float, 6>> err;  // xx, xy, yy, xz, yz, zz
      std::vector<std::uint32_t> detid;
      std::vector<std::uint16_t> det_in_layer, charge_pcm;
      std::vector<std::uint8_t> span_rows, span_cols;
    };
    std::vector<LayerHits> layers;

    // mkfit::BeamSpot of the EventOfHits.
    float bs_x = 0, bs_y = 0, bs_z = 0, bs_sigma_z = 0, bs_width_x = 0, bs_width_y = 0, bs_dxdz = 0, bs_dydz = 0;
    std::int32_t n_hits_pixel = 0, n_hits_outer = 0;
    // The dead bins of the EventOfHits' layers, LayerOfHits::isBinDead(phi bin, q bin).
    std::int32_t n_phi_bins = 0;
    std::vector<std::int32_t> dead_layer, dead_q_bin, dead_phi_bin;

    std::string geom_md5, geom_version;
    // The bytes of the geometry file; written into Geometry when its md5 is not there yet. May be null.
    const std::vector<std::uint8_t>* geom_bin = nullptr;
  };

  // The origin of each hit's charge in one traced event, for TraceRNTupleWriter::writeHitOrigin(); the RNTuple
  // TruthHitOrigin is described above.
  struct TraceHitOriginTables {
    unsigned int run = 0, lumi = 0;
    unsigned long long event = 0;
    std::string module;  // label of the MkFitProducer the table is for

    enum Class : std::uint8_t { kNone = 0, kOwn = 1, kSecondary = 2, kUnknown = 3 };

    struct LayerHits {
      int layer = -1;
      int n_hits = 0;                 // 1 + the largest index of a hit in the layer
      std::vector<std::int32_t> hit;  // indices of the hits with simulated charge
      std::vector<std::int32_t> sim;  // the sim track of the signal SimTrack with most charge, -1 none
      std::vector<float> f_own, f_sec, f_pu, f_unm;
      std::vector<std::uint8_t> cls;
      std::vector<std::int32_t> sec_pdg, sec_proc;
      // The pileup's part (see TruthHitOrigin above).
      std::vector<float> f_pu_own, f_pu_sec, f_pu_unm;
      std::vector<std::uint8_t> cls_all;
      std::vector<std::int32_t> sim_all, sec_pdg_all, sec_proc_all;
    };
    std::vector<LayerHits> layers;
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

    // Appends the truth tables of one event, in the caller, into the RNTuples of truth.source (created with the
    // source's first event). Returns their pass. Thread safe. Throws std::runtime_error for a source name that is not
    // letters and digits.
    int writeTruth(const TraceTruthTables& truth);

    // Appends the input tables of one event, in the caller, and the geometry when its md5 is new to the file.
    // Returns their pass; geomWritten tells whether the geometry was written. Thread safe.
    int writeInputs(const TraceInputTables& inputs, bool* geomWritten = nullptr);

    // Appends the hit-origin table of one event, in the caller. Returns its pass. Thread safe.
    int writeHitOrigin(const TraceHitOriginTables& origin);

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
