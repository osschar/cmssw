#ifndef RecoTracker_MkFitCore_interface_TraceStore_h
#define RecoTracker_MkFitCore_interface_TraceStore_h

//==============================================================================
// The trace of MkFinderV2p2 (MKFIT_TRACE): its records and the store the hooks
// write them into.
//
// MkFitCore includes this header only under MKFIT_TRACE. The standalone's ROOT
// data formats (mkFit-external, standalone/DataFormats/RntStructs.h) include it
// in every ROOT build, as they did when the records were defined there.
//
// Positions and momenta in the records are ROOT's REveVector.
//==============================================================================

#include "RecoTracker/MkFitCore/interface/TrackState.h"

#include "ROOT/REveVector.hxx"

#include <cassert>
#include <functional>
#include <string_view>
#include <utility>
#include <vector>

//==============================================================
// Basic kine structs
//==============================================================

using EVec3 = ROOT::Experimental::REveVector;

struct EBiVec3 {
  EVec3 pos, mom;
};

struct PropInfo : public EBiVec3 {
  float dalpha;  // helix angle during propagation
  int fail_flag;
};

//==============================================================
// Trace and RDF stuff
//==============================================================

struct TrCandMeta {
  int id;
  int event;
  int seed = -1;  // index of seed in Event::trSeeds_
  // The seed's label. It equals the seed's index in Event::seedTracks_ only after
  // Event::relabelSeedTracksSequentially().
  int global_seed = -1;
  int sim = -1;                  // index of sim track in event simTracks_
  int cand = -1;                 // candidateTracks_ index (valid for final stage)
  int stage_ids[3]{-1, -1, -1};  // indices of stages
  // The forward search's final candidate: its worst end-of-layer selection rank
  // along its ancestry, that layer, and the first layer where the rank exceeded
  // Config::V2p2::Diag::rank_over (TrackCand::m_rank_*).
  int fwd_rank_max = 0;
  int fwd_rank_max_layer = -1;
  int fwd_rank_over_layer = -1;
  int fwd_stop_reason = 0;               // its TrackCand::StopReason_e, 0 = not stopped
  unsigned long long fwd_hits_hash = 0;  // hits_hash() (src/TraceConversions.h) of its hits
};

struct TrCandStage {
  int id = -1;
  int meta_id = -1;          // link to global meta
  int parent_stage_id = -1;  // previous stage (for multi-stage linkage)

  int stage = 0;          // 0=FwdSearch, 1=BkwFit, 2=BkwSearch, as in SteeringParams::IterationType_e
  int iteration_idx = 0;  // which iteration

  int root_state_id = -1;   // first state in this stage
  int final_state_id = -1;  // best final state in this stage
};

struct TrCandState {
  int id = -1;
  int parent_id = -1;
  int meta_id = -1;
  int stage_id = -1;

  int layer;
  int step;  // not sure -- for now just depth number, incresed for every new state

  // missing: n_hits, n_missed_hits

  EBiVec3 kine;
  mkfit::TrackState state;  // TrackState at this point

  int search_id = -1;          // set in post processing (might need multiple searches for thick layers)
  bool has_children = false;   // true if it has children
  bool on_final_path = false;  // true if it is on the final selected candidate path (final candidate and its ancestors)
};

// TrLayerSearch -- one hit-search window, for one candidate state, on one layer.
//
// Filled by MkFinderV2p2::process_layer_batch(), stage 1 ("initial propagation"):
// the mini-propagator is run to the two layer bounding surfaces, the full
// propagation supplies the covariance, and the two together give the phi/q
// window and the bin ranges that the stage-2 hit loop then scans.
//
// This is the v2p2 counterpart of the old RntDumper BinSearch plus the two
// PropInfos in CandInfo (ps_min / ps_max) -- keep them comparable.
struct TrLayerSearch {
  int id = -1;
  int state_id = -1;
  int layer = -1;
  int layer_sec = -1;  // second sub-layer of a double layer, -1 when not double
  bool is_barrel = false;
  bool is_outward = false;  // false for the inward search into the pixels

  // The within-sensitive-region verdict (TrackerInfo.h WithinSensitiveRegion_e:
  // -1 undef, 0 inside, 1 edge, 2 outside), from MkFinderV2p2::determine_wsr()
  // on the sp1/sp2 crossings fuzzed by 5 sigma of dq_track. RECORD-ONLY here --
  // the search itself skips a WSR_Outside candidate entirely when
  // Config::V2p2::Policy::use_wsr is on, but the record is written either way, so any
  // analysis over layer searches MUST filter on this or it counts searches that
  // never scanned a hit. in_gap says the absence is explained by a disc's
  // small-r hole rather than by the track missing the layer.
  signed char wsr = -1;
  bool wsr_in_gap = false;

  // Stage 1a -- mini-propagator (PA_Exact) onto the layer bounding surfaces,
  // in propagation order: entry is crossed first, exit second. NOT sorted by
  // q or phi. dalpha is the helix angle turned, fail_flag the propagation status.
  // The Hermite cubic used for the per-hit propagation is built from these two,
  // so its interpolation error scales as (dalpha_exit - dalpha_entry)^4.
  PropInfo prop_entry;  // MkBins::m_sp1
  PropInfo prop_exit;   // MkBins::m_sp2

  // Stage 1b -- search window from MkBins::determine_bin_windows().
  // phi_delta, dphi_track and dq_track are HALF-widths; q_min/q_max a full range.
  float phi_center = -999.99f, phi_delta = -999.99f;
  float q_min = -999.99f, q_max = -999.99f;
  float q_center() const { return 0.5f * (q_min + q_max); }
  float dphi_track = -999.99f, dq_track = -999.99f;  // 3 sigma track errors

  // Covariance the window was built from (MkBinTrackCovExtract), i.e. err(0,0),
  // err(0,1), err(1,1), err(2,2) after the full propagation to prop_exit.
  // NOTE: the phi jacobian is evaluated at min(r_entry, r_exit) -- the smallest
  // radius the track crosses in the layer, where sigma_phi is largest, so one
  // window covers the whole layer. It used to be taken at m_isp (the entry edge)
  // while the covariance came from m_sp2 (the exit), which scaled dphi_track by
  // rout/rin: 1.33 at pixel layer 0. The ENDCAP dq jacobian is still at m_isp,
  // but |grad r| = 1 there so there is no scale factor. Recorded
  // here so the size of that mismatch can actually be measured.
  float cov_xx = -999.99f, cov_xy = -999.99f, cov_yy = -999.99f, cov_zz = -999.99f;

  // Stage 1c -- bin ranges. p2/q2 are exclusive; p wraps around, q does not.
  unsigned short p1 = 0, p2 = 0, q1 = 0, q2 = 0;                  // primary layer
  unsigned short p1_sec = 0, p2_sec = 0, q1_sec = 0, q2_sec = 0;  // secondary sub-layer

  // Stage 2 tallies over the hit loop. n_hits_pqueue is capped at
  // V2p2Params::max_presel_hits per sub-layer and equals the number of
  // TrHitMatch with passed_pqueue.
  int n_hits_scanned = 0;  // hits visited inside the bin ranges
  int n_hits_masked = 0;   // of those, rejected by the iteration hit mask
  int n_hits_precut = 0;   // of the rest, rejected by the line pre-cut
  int n_hits_presel = 0;   // of those, passed the dq/dphi pre-selection
  int n_hits_pqueue = 0;   // of those, survived the priority queue -> Kalman

  // Only computed with MKFIT_TRACE_PROP_COMPARE, otherwise left at zero. The
  // member itself is unconditional: unlike the per-hit compare state this is one
  // vector per search, not per scanned hit, and keeping it out of the #ifdef
  // means one less layout the ROOT dictionary has to agree with.
  //
  // How far the Hermite cubic strays from the helix in the middle of the layer,
  // where it is furthest from the two endpoints it was built to interpolate:
  // evaluate the cubic at t = 0.5, then propagate the mini-propagator exactly to
  // the bounding surface (r or z) that point lands on, and difference the two.
  // Both are on the same surface, so this is a pure in-surface deviation.
  EVec3 hermite_mid_dev{0, 0, 0};
};

struct TrHitMatch {
  int id;
  int state_id;
  int search_id = -1;  // TrLayerSearch this hit was scanned for
  int layer;           // NOTE: TrCandState.layer we are pointing to is PREVIUOUS layer
  int hit;
  bool mc_match;  // is hit mc-matching

  // hit info -- can get from event->layer-of-hits ....
  // EVec3 hit_pos;
  // float hit_phi;
  // float hit_q;
  // float hit_qbar;
  // int   hit_lbl;

  // pre-selection quantities.
  // kine_on_plane is the state the dphi / dq / residual_* below were computed
  // from -- currently the mini-propagator PA_Line step onto the module plane.
  EBiVec3 kine_on_plane{EVec3(), EVec3()};
  // Hermite parameter at the module-plane crossing. Inside [0,1] the crossing is
  // between the layer bounding surfaces; outside it the cubic is extrapolating
  // and both kine_on_plane and the window it was pre-selected against are suspect.
  float t_hermite = -999.99;
  // Distance from the Hermite point at t_hermite to the module plane, in cm.
  // Hermite3DOnPlane::solve() takes a SINGLE Newton step, so this is not
  // machine-zero by construction -- it is the measure of whether one step
  // sufficed for this hit.
  float d_plane_h3 = -999.99;
  float dphi = -999.99;
  float dq = -999.99;
  // Half-length of the hit in the q direction, i.e. half the strip length in the
  // barrel. Needed to make dq interpretable: for a long strip the residual is
  // dominated by this, not by the track error, so dq / sigma_track is NOT a
  // covariance pull there. It is also the second term of the dq cut.
  float hit_q_half_len = -999.99;
  float score = -999.99;  // usually dphi, could add dq checks
  bool passed_preselect = false;

  float residual_x = -999.99;  // distance from hit to track in precise / phi direction, for accepted hits
  float residual_y = -999.99;  // distance from hit to track in coarse direction, for accepted hits
  float residual_z = -999.99;  // distance from detector plane, should be about 0

  // Which SOURCE this hit came from. Pre-selection runs per sub-layer -- one
  // bounded pqueue each -- so that a shower in one sensor cannot evict the
  // other's hits, and so that each sub-layer's pre-selection quality can be
  // measured on its own. Do NOT infer this from `layer` once the split layers
  // are merged into fat layers: both sensors are then in ONE layer and the layer
  // number stops distinguishing them. This flag survives that.
  bool is_sec_layer = false;

  // after pre-selection via priority-queue
  //
  // TWO ranks, because there are two populations and they answer different
  // questions. Both are by score (ddphi), 1 = best, -1 = did not pass the pqueue.
  //   sub_rank  -- within this hit's OWN sub-layer. Measures that sub-layer's
  //                pre-selection quality, and is what the single `rank` field
  //                meant back when only one sub-layer was ever scanned.
  //   full_rank -- across the merged layer, both sub-layers together. Says
  //                whether the layer's best hit sits in P or in S.
  // NEITHER is the position in the step-ordered list the combinatorial search
  // walks: that is ordered by path length, not by score, and is a third thing.
  int sub_rank = -1;  // (would be nice to have for rejected hits too, std::vec instead of pqueue for MKFIT_TRACE?)
  int full_rank = -1;
  bool passed_pqueue = false;

  int kalman_id = -1;

  // Optional -- with MKFIT_TRACE_PROP_COMPARE.
  // The OTHER of the two cheap propagations onto the module plane: the code runs
  // both a mini-propagator PA_Line step and a Hermite-cubic plane solve for every
  // scanned hit, and uses one of them (see kine_on_plane) for the cuts. This is
  // the one it did not use, so the two can be differenced for ALL scanned hits,
  // including rejected ones.
  // For the pre-selected hits the exact reference is already in the trace and
  // needs no define: TrKalmanUpdate::propagated_state, via kalman_id.
  // Set the macro in Makefile.config -- it must be in CPPFLAGS so that the code
  // and the ROOT dictionary agree on the layout.
#ifdef MKFIT_TRACE_PROP_COMPARE
  EBiVec3 kine_on_plane_cmp{EVec3(), EVec3()};
#endif
};

// One hit of the BACKWARD FIT. Deliberately NOT TrHitMatch + TrKalmanUpdate:
// those are split because in the SEARCH one state scans many hits and only a
// few reach a Kalman update, and the split carries that fan-out. The fit has no
// fan-out -- one state, one hit, one update, always -- so a single record at
// 1:1 with hits collapses them and costs no join.
//
// The fit is a CHAIN, not a tree: TrCandStates linked by parent_id with one of
// these hanging between each consecutive pair. Its stage is TrCandStage with
// stage == 1 (BkwFit), which is also the discriminator telling a consumer to
// read this collection rather than trHitMatches_ / trKalmanUpdates_.
struct TrBkFitUpdate {
  int id;
  int state_id_in = -1;   // TrCandState before this update
  int state_id_out = -1;  // TrCandState after it; the fit rejects nothing, so
                          // this is always set, and the updated state is
                          // reachable there -- hence no updated_state member.
  int layer = -1;
  int hit = -1;  // index in the layer
  // simHitsInfo_[hit.mcHitID()].mcTrackID(), i.e. WHICH PARTICLE made this hit.
  // Stored as the id and NOT as a "matches the track" bool: the comparison is
  // against the track's SIM label, which is not known at fill time -- the only
  // label MkFinder has is the seed's, and after
  // relabelSeedTracksSequentially() that is the seed's INDEX, a different
  // namespace. Resolve it in analysis via TrCandMeta::global_seed.
  int mc_track_id = -2;  // -2 = not resolved
  int step = -1;         // hit index along the fit, 0 = outermost (first)

  // Propagation status of the propagation to the plane. The fit SOLVES for the
  // plane with getS (propagateHelixToPlaneMPlex).
  // NOTE it is currently ALWAYS 0: Leonardo's getS carries no bracket and so
  // detects no failure. Non-zero values only appear if a solver with failure
  // detection is put back (the Hermite bracket did, and that is how the fit's
  // missing state gate was found). Do not read 0 as "the fit did not fail".
  int fail = 0;

  // step == 0 is a propagation to the plane the state is ALREADY on -- the seed
  // sits on the outermost hit -- so it is zero-length and its chi2 re-uses the
  // hit that produced the state. EXCLUDE step 0 from any chi2 measurement.

  float chi2 = -999.99f;
  float chi2_cum = -999.99f;  // running total BEFORE this hit

  // Propagated state on the module plane, and the residual to the hit decomposed
  // in the MODULE frame: x = plDir (across strip, precise), y = plNrm x plDir
  // (along strip, coarse), z = plNrm (normal, must be ~0).
  EBiVec3 kine_on_plane{EVec3(), EVec3()};
  float residual_x = -999.99f;
  float residual_y = -999.99f;
  float residual_z = -999.99f;

  // Optional -- with MKFIT_TRACE_KALMAN_DEBUG. Only the PRE-update state needs
  // a member: the post-update one is trCandStates_[state_id_out].state.
#ifdef MKFIT_TRACE_KALMAN_DEBUG
  mkfit::TrackState propagated_state{};
#endif
};

struct TrKalmanUpdate {
  int id;
  int hit_match_id;
  int state_id_in;        // state before update; could get it from hitmatch.state_id
  int state_id_out = -1;  // state after update, -1 if not accepted

  float chi2 = -999.99f;
  float chi2_trk = -999.99f;
  // The acceptance chi2 (V2p2Params::chi2_trk_fac on the track covariance), which the
  // cut V2p2Params::hit_chi2_cut is applied to. v2p2 only.
  float chi2_acc = -999.99f;
  // This hit advanced the state. In the in-layer combinatorial search, that is:
  // the node is on a TrackCand that survived the end-of-layer selection, so a
  // hit that passed the chi2_acc cut can still be not accepted.
  bool accepted = false;

  // Propagated track parameters are same as in TrHitMatch (as we pass it to Kalman update),
  // The local dx, dy, dz distance from point to hit is the same as for hit match.

  // Optional -- with MKFIT_TRACE_KALMAN_DEBUG.
  // For an accepted update the post-update state is also reachable as
  // trCandStates_[state_id_out].state; this member is what gives it for the
  // REJECTED ones too, which is why it is worth the space when studying chi2
  // and error-matrix behaviour. Set the macro in Makefile.config -- it must be
  // in CPPFLAGS so that the code and the ROOT dictionary agree on the layout.
#ifdef MKFIT_TRACE_KALMAN_DEBUG
  mkfit::TrackState propagated_state{};  // pre-update, propagated onto the hit plane
  mkfit::TrackState updated_state{};     // post-update
#endif
};

namespace mkfit {

  //==============================================================================
  // TraceTruthSource -- what the hooks need beyond the records: the event
  // number, the seeds' labels and the sim truth of seeds and hits. The four
  // tr_* calls the hooks make are here; tr_mc_hit_id() is not asked by any hook
  // (they read the hit's mcHitID themselves) and serves apply_trace_truth(). The
  // defaults answer "unknown".
  //==============================================================================

  class TraceTruthSource {
  public:
    virtual ~TraceTruthSource() = default;

    // The event number written into TrCandMeta::event.
    virtual int tr_event_id() const { return -1; }
    // The label of seed seed_idx of the current search, for TrCandMeta::global_seed.
    virtual int tr_seed_label(int seed_idx) const { return -1; }
    // Whether the hit with sim index mc_hit_id was made by the sim track of seed
    // seed_idx, for TrHitMatch::mc_match.
    virtual bool tr_mc_match(int seed_idx, int mc_hit_id) const { return false; }
    // The sim track that made the hit with sim index mc_hit_id, -2 when not
    // resolved, for TrBkFitUpdate::mc_track_id.
    virtual int tr_mc_track_id(int mc_hit_id) const { return -2; }
    // The sim index (Hit::mcHitID()) of hit hit_idx of layer layer, the index
    // the records keep (TrHitMatch::hit, TrBkFitUpdate::hit); -1 when not known.
    virtual int tr_mc_hit_id(int layer, int hit_idx) const { return -1; }
  };

  //==============================================================================
  // TraceStore -- the seven record vectors of one event and the calls the hooks
  // make to append and reach records. Each append sets the record's id to its
  // index. Not thread safe within an event: one store per event in flight. The
  // vectors are mutable and the calls const, so that a hook holding a const
  // pointer can append.
  //
  // The standalone Event derives from TraceStore. The hooks ask the store the
  // questions of TraceTruthSource. A plain TraceStore answers "unknown", and the
  // standalone Event answers from its own sim information.
  //==============================================================================

  class TraceStore : public TraceTruthSource {
  public:
    TraceStore() = default;
    TraceStore(const TraceStore &) = default;
    TraceStore(TraceStore &&) = default;
    TraceStore &operator=(const TraceStore &) = default;
    TraceStore &operator=(TraceStore &&) = default;
    ~TraceStore() override = default;

    // ---- The records.

    mutable std::vector<TrCandMeta> trCandMetas_;
    mutable std::vector<TrCandStage> trCandStages_;
    mutable std::vector<TrCandState> trCandStates_;
    mutable std::vector<TrLayerSearch> trLayerSearches_;
    mutable std::vector<TrHitMatch> trHitMatches_;
    mutable std::vector<TrKalmanUpdate> trKalmanUpdates_;
    mutable std::vector<TrBkFitUpdate> trBkFitUpdates_;

    void clear_trace() {
      trCandMetas_.clear();
      trCandStages_.clear();
      trCandStates_.clear();
      trLayerSearches_.clear();
      trHitMatches_.clear();
      trKalmanUpdates_.clear();
      trBkFitUpdates_.clear();
    }

    TrCandMeta &tr_candmeta(int i) const { return trCandMetas_[i]; }
    TrCandStage &tr_candstage(int i) const { return trCandStages_[i]; }
    TrCandState &tr_candstate(int i) const { return trCandStates_[i]; }
    TrLayerSearch &tr_layersearch(int i) const { return trLayerSearches_[i]; }
    TrHitMatch &tr_hitmatch(int i) const { return trHitMatches_[i]; }
    TrKalmanUpdate &tr_kalmanupdate(int i) const { return trKalmanUpdates_[i]; }
    TrBkFitUpdate &tr_bkfitupdate(int i) const { return trBkFitUpdates_[i]; }

    TrCandMeta &trace_candmeta(TrCandMeta &&cm) const {
      int s = trCandMetas_.size();
      auto &t = trCandMetas_.emplace_back(cm);
      t.id = s;
      return t;
    }
    TrCandStage &trace_candstage(TrCandStage &&cs) const {
      int s = trCandStages_.size();
      auto &t = trCandStages_.emplace_back(cs);
      t.id = s;
      return t;
    }
    TrCandState &trace_candstate(TrCandState &&cs) const {
      int s = trCandStates_.size();
      auto &t = trCandStates_.emplace_back(cs);
      t.id = s;
      return t;
    }
    TrLayerSearch &trace_layersearch(TrLayerSearch &&ls) const {
      int s = trLayerSearches_.size();
      auto &t = trLayerSearches_.emplace_back(ls);
      t.id = s;
      return t;
    }
    TrHitMatch &trace_hitmatch(TrHitMatch &&hm) const {
      int s = trHitMatches_.size();
      auto &t = trHitMatches_.emplace_back(hm);
      t.id = s;
      return t;
    }
    TrKalmanUpdate &trace_kalmanupdate(TrKalmanUpdate &&ku) const {
      int s = trKalmanUpdates_.size();
      auto &t = trKalmanUpdates_.emplace_back(ku);
      t.id = s;
      return t;
    }
    TrBkFitUpdate &trace_bkfitupdate(TrBkFitUpdate &&bu) const {
      int s = trBkFitUpdates_.size();
      auto &t = trBkFitUpdates_.emplace_back(bu);
      t.id = s;
      return t;
    }

    int trace_new_cand_meta(int event, int seed_index) const {
      auto &cm = trace_candmeta({-1, event, seed_index});
      return cm.id;
    }
    // always do stage and initial state together
    // int trace_new_cand_stage(int meta_id, int parent_stage_id, int stage) {
    //   auto &cstg = trace_candstage({ -1, meta_id, parent_stage_id, stage });
    //   return cstg.id;
    // }
    std::pair<int, int> trace_new_cand_stage_and_state(
        int meta_id, int parent_stage_id, int stage, int layer, const EBiVec3 &kine, const TrackState &state) const {
      assert(stage >= 0 && stage <= 2 && "stage expected to be between 0 and 2");
      auto &cstage = trace_candstage({-1, meta_id, parent_stage_id, stage});
      auto &cstate = trace_candstate({-1, -1, meta_id, cstage.id, layer, 0, kine, state});
      cstage.root_state_id = cstate.id;
      return {cstage.id, cstate.id};
    }
    int trace_new_cand_state(int parent_state_id, int layer, const EBiVec3 &kine, const TrackState &state) const {
      auto &pcs = trCandStates_[parent_state_id];
      pcs.has_children = true;
      auto &cs = trace_candstate({-1, parent_state_id, pcs.meta_id, pcs.stage_id, layer, pcs.step + 1, kine, state});
      return cs.id;
    }
    int trace_new_kalman_update(int hit_match_id, int state_id_in, float chi2, float chi2_trk) const {
      auto &ku = trace_kalmanupdate({-1, hit_match_id, state_id_in, -1, chi2, chi2_trk});
      return ku.id;
    }
  };

  //==============================================================================
  // apply_trace_truth -- sets the record fields the hooks fill from the store's
  // TraceTruthSource answers, from truth instead, through the same four calls
  // and from each record's own content:
  //   TrCandMeta::event        = tr_event_id(), every meta;
  //   TrCandMeta::global_seed  = tr_seed_label(meta.seed), every meta with
  //                              cand >= 0: the export hook asks for the label
  //                              exactly when it sets cand
  //                              (MkBuilder::export_best_comb_cands()); the other
  //                              metas keep theirs;
  //   TrHitMatch::mc_match     = tr_mc_match(seed, tr_mc_hit_id(layer, hit)), with
  //                              seed the meta's seed of the match's state;
  //   TrBkFitUpdate::mc_track_id = tr_mc_track_id(tr_mc_hit_id(layer, hit)).
  // A trace recorded through a store without truth and then given the same
  // truth here equals the trace recorded with that truth. This holds as long as
  // the export hook is the last writer of global_seed for an exported meta (no
  // backward fit traced after the export in the same event).
  // Defined in src/TraceStore.cc (MKFIT_TRACE builds only).
  //==============================================================================

  void apply_trace_truth(TraceStore &store, const TraceTruthSource &truth);

  class Track;

  //==============================================================================
  // The trace sink -- how a store leaves the event. The CMSSW producer calls the
  // registered sink at the end of each event with its stream's store and what
  // it knows about the event; a writer outside MkFitCore registers it. MkFitCore
  // itself does no I/O. Defined in src/TraceStore.cc (MKFIT_TRACE builds only).
  //
  // Not thread safe: set the sink before the first event and clear it after the
  // last one. The sink itself is called from every stream and must lock.
  //==============================================================================

  struct TraceEventInfo {
    unsigned int run = 0;
    unsigned int lumi = 0;
    unsigned long long event = 0;
    std::string_view module;                         // label of the producing module
    int n_seeds_in = 0;                              // seeds given to mkFit
    const std::vector<Track> *seeds = nullptr;       // after cleaning: TrCandMeta::seed indexes these
    const std::vector<Track> *out_tracks = nullptr;  // mkFit's output, after the duplicate cleaner
  };

  using TraceSink = std::function<void(const TraceStore &, const TraceEventInfo &)>;

  void set_trace_sink(TraceSink sink);
  const TraceSink &trace_sink();

}  // namespace mkfit

#endif
