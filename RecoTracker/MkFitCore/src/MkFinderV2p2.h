#ifndef RecoTracker_MkFitCore_src_MkFinderV2p2_h
#define RecoTracker_MkFitCore_src_MkFinderV2p2_h

#include "RecoTracker/MkFitCore/interface/SteeringParams.h"
#include "RecoTracker/MkFitCore/interface/HitStructures.h"
#include "RecoTracker/MkFitCore/interface/TrackStructures.h"
#include "RecoTracker/MkFitCore/interface/MkJob.h"

#include "MkRZLimits.h"
#include "MkFinderV2p2Structures.h"

#include "MkBins.h"

#include <atomic>
#include <functional>
#include <list>
#include <utility>
#include <vector>

#if defined(MKFIT_STANDALONE)
#include "RecoTracker/MkFitCore/standalone/V2p2Diag.h"
#else
#define V2P2_COUNT(field) ((void)0)
#define V2P2_COUNT_ADD(field, n) ((void)0)
#endif

namespace mkfit {

  // Switches and parameters are in V2p2Config.h, namespace Config::V2p2.

  // Per-layer policy counters, standalone only: see standalone/V2p2Diag.h. The
  // finder increments them through V2P2_COUNT() and V2P2_COUNT_ADD(), which expand
  // to nothing in the CMSSW build.

  class FindingFoos;
  class IterationParams;
  class IterationLayerConfig;
  class SteeringParams;
  struct LayerControl;
  class LayerInfo;
  class Event;

  class MkJob;

  class MkFinderV2p2 {
    friend class MkBuilder;

    //-------------------------------------------------------------------------
    class BatchManager {
    //-------------------------------------------------------------------------
      friend class MkFinderV2p2;

      EventOfCombCandidates *mp_eoccs = nullptr;
      int m_begin;
      int m_end;

      int m_n_dormant;
      int m_n_finding;
      // int m_n_to_finalize;
      int m_n_finished;

      // Not needed, have the list of CCandReps
      // // Cursors into EOCCS for activation
      // int m_Cc_pos;  // current CombCand index
      // int m_pTc_pos; // current primary TrackCand index

      // something for the derived ones --- teritary, to handle in-layer combinatorials
      // or maybe we'll need a sub-manager for those? Or will the finder do that
    public:
      CombCandidate& ccand(int i) const { return (*mp_eoccs)[i]; }
      CombCandidate* ccand_ptr(int i) const { return &(*mp_eoccs)[i]; }

      void setup(EventOfCombCandidates &eoccs, int seed_begin, int seed_end) {
        mp_eoccs = &eoccs;
        m_begin = seed_begin;
        m_end = seed_end;
        m_n_dormant = 0;
        m_n_finding = 0;
        m_n_finished = 0;
        reset_for_new_layer();
      }
      void release() {
        mp_eoccs = nullptr;
      }

      void reset_for_new_layer() {
        // Prepare for next layer / extraction.
        // We have list of active CCandReps now
        // m_Cc_pos = m_begin;
        // m_pTc_pos = 0;
      }

      int n_total() const { return m_end - m_begin; }
      int n_dormant() const { return m_n_dormant; }
      int n_finding() const { return m_n_finding; }
      int n_finished() const { return m_n_finished; }
      bool are_all_ccands_finished() const { return n_finished() == n_total(); }
      bool has_dormant_ccands() const { return n_dormant() > 0; }

      // Iteration over all CombCandidates for top-level administrative tasks.
      class iterator {
        CombCandidate *m_ccand;
      public:
        iterator(CombCandidate *bm) : m_ccand(bm) {}
        CombCandidate& operator*() { return *m_ccand; }
        iterator& operator++() {
          ++m_ccand;
          return *this;
        }
        bool operator!=(const iterator &i) const { return m_ccand != i.m_ccand; }
      };

      iterator begin() const { return iterator(ccand_ptr(m_begin)); }
      iterator end() const { return iterator(ccand_ptr(m_end - 1) + 1); }
    //-------------------------------------------------------------------------
    }; // end class BatchManager
    //-------------------------------------------------------------------------

  public:
    MkFinderV2p2() = default;

    //----------------------------------------------------------------------------

    void setup(const MkJob *job, EventOfCombCandidates &eoccs, int seed_begin, int seed_end,
               SteeringParams::iterator &sp_it, const Event *ev);
              //  int region, // sp_it->region()
              //  const PropagationConfig &pc, // trk_info.prop_config() and m_job->m_trk_info
              //  const IterationConfig &ic, // m_job->m_iter_config
              //  const IterationParams &ip, // m_job->params_cur();
              //  const IterationLayerConfig &ilc, // m_job->m_iter_config.m_layer_configs[curr_layer]
              //  const SteeringParams &sp, // m_job->steering_params(region)
              //  const std::vector<bool> *ihm, // m_job->get_mask_for_layer(curr_layer)
              //  bool infwd);               // m_job->m_in_fwd
    void release();

    int awaken_candidates();

    void begin_layer();

    bool any_Ccreps_to_begin() const { return m_active_ccreps_pos != m_active_ccreps.end(); }
    void begin_next_Ccrep_in_layer();

    int  cand_queue_size() const { return (int) m_cand_queue.size() - m_cand_queue_head; }
    bool enough_work_for_batch() const { return cand_queue_size() >= NN; }
    bool any_work_for_batch() const { return cand_queue_size() > 0; }
    void process_layer_batch();

    void end_layer();

    void process_layer();

    // Mostly for debug printouts in MkBuilder steering code.
    const BatchManager& batch_mgr() const { return m_batch_mgr; }

    //----------------------------------------------------------------------------

  private:
    //----------------------------------------------------------------------------
    // Per-pass state of process_layer_batch(). LayerBatch is NN candidates wide
    // and indexed by i; HitBatch is NN (candidate, hit) pairs wide and indexed by
    // h, with prim_idcs[h] naming the candidate.
    struct LayerBatch {
      int N_proc = 0;
      PrimTCandRep *ptc[NN];             // the candidates in this pass
      MkBins B { 0 };                    // isp + the two bounding-surface crossings
      MkBinTrackCovExtract TCE;          // position block of the window covariance
      MkBinLimits BL_p, BL_s;            // binnor ranges, primary / secondary layer
      mini_propagators::Hermite3D H;     // cubic through sp1, sp2 -- the trajectory model
      // Local hit density per candidate, ln(hits / cm^2), over the bins walked;
      // for the likelihood score.
      int   n_scanned[NN] = {0};
      float rho_dphi[NN] = {0.0f};  // Score::rho_region 2: the cut region's padding, max over sub-layers
      float rho_dq[NN] = {0.0f};
      float log_rho[NN] = {0.0f};
#ifdef MKFIT_TRACE
      int tr_layersearch_ids[NN];
#endif
    };

    struct HitBatch {
      int fill_pos = 0;
      MPlexQI  prim_idcs;                // -> LayerBatch::ptc
      MPlexQUI hit_idcs;                 // index within the LayerOfHits
      MPlexQUI hit_orig_idcs;            // index in the event hit vector
      mini_propagators::InitialStatePlex is_plex;
#ifdef MKFIT_TRACE_PROP_COMPARE
      mini_propagators::StatePlex h_plex;  // PA_Line cross-check against the Hermite
#endif
    };

    void prop_to_layer_edges(LayerBatch &b);
    void determine_search_windows(LayerBatch &b);
    void determine_wsr(LayerBatch &b);
    void select_hits(LayerBatch &b);
    void preselect_hit_batch(LayerBatch &b, HitBatch &hb, const LayerOfHits &L, int N_proc_hits,
                             bool is_sec_layer);
    // The track's q error referenced to the hit's module plane. See MkFinderV2p2.cc.
    static float surface_referenced_dq(float dq_track_fallback,
                                       const MkBinTrackCovExtract &TCE, int pi,
                                       const mini_propagators::StatePlex &h3_state, int h,
                                       const MPlex3V &module_norm, bool is_barrel);
    // Candidate-stopping cuts applied at pull-in, where only the candidate's own
    // state is needed. Returns the reason, or SR_NotStopped to keep going.
    TrackCand::StopReason_e stop_cuts_at_pickup(const TrackCand &tc) const;
    // Which fake HoT a candidate that took no hit in this layer gets: the hole
    // limits, the WSR override, and the gap case, in V1's order.
    int fake_hit_index(const TrackCand &tc, const WSR_Result &wsr) const;

    void prepare_kalman_workload(LayerBatch &b);
    void kalman_update(LayerBatch &b);
    void process_kalman_results(LayerBatch &b);

    // The in-layer combinatorial search, replacing the two phases above when
    // Config::V2p2::InLayer::comb is on. expand_in_layer() grows the SecTCandRep
    // tree breadth-first by depth; select_and_materialise() runs the end-of-layer
    // selection and registers the survivors into the CombCandidate.
    void expand_in_layer(LayerBatch &b);
    void select_and_materialise(CCandRep &ccrep);
    void drop_sister_dominated_holes(const CombCandidate &ccand);
#if defined(MKFIT_STANDALONE)
    void count_dominated_kept(const CombCandidate &ccand, int n_keep) const;
#endif
    bool long_step_to_layer(const TrackCand &tc) const;
    void offer_best_short(CombCandidate &ccand, const TrackCand &tc) const;
    // The direction-, layer- and candidate-dependent part of a layer step, filled
    // once per path root and carried down the tree.
    void fill_step_geometry(LayerStepFeatures &f, const PrimTCandRep &ptc, float log_rho) const;
    // Turn the Kalman results accumulated in m_sec_out into arena nodes, keeping
    // those that pass the chi2 cut. Their indices go to m_sec_harvested.
    void harvest_sec_nodes();
    // A slot for one node: from the free list if it has any, else appended.
    int sec_node_slot();
    // After a layer batch: select every CombCandidate of the batch that has no
    // candidate left in the queue, and free its nodes.
    void select_completed_ccreps(const LayerBatch &b);

    // The running variant of expand_in_layer(), InLayer::running_kalman. Depth 0
    // (the Hermite crossing is known, sPerp path) and deeper steps (the plane is
    // solved) need different propagation, so there are two batches. Each fires
    // only when full; rk_drain() empties both at end of layer.
    void expand_in_layer_running(LayerBatch &b);
    void rk_setup(KalmanOpArgs &koa, bool solve_plane);
    void rk_flush(KalmanOpArgs &koa);
    void rk_add_item(KalmanOpArgs &koa, PrimTCandRep &ptc, const PrimTCandRep::PQE &pqe, int parent_idx,
                     int hit_pos);
    void rk_expand_node(int ni);
    void rk_pump();
    void rk_drain();
    void load_hit_module(KalmanOpArgs &koa, const PrimTCandRep::PQE &pqe) const;

    //----------------------------------------------------------------------------
    // Job / batch-of-seeds control variables and globel references
    const MkJob *mp_job = nullptr;
    SteeringParams::iterator *mp_steeringparams_iter = nullptr;
    const Event *mp_event = nullptr;

    BatchManager m_batch_mgr;

    //----------------------------------------------------------------------------
    // Per-(di)layer state & control
    std::list<CCandRep> m_active_ccreps;
    // Current CombCand and next TrackCand to go through layer initialization, i.e.,
    // propagation to layer limits, Binnor creation and extraction of bin-indices, and
    // pre-selection of hits.
    std::list<CCandRep>::iterator m_active_ccreps_pos; // Current CombCand to be processed or is in processing.

    // Pre-selection queue -- list of pTcs to do initial prop + Binnor + hit extraction for.
    // Elements are slots in the pTC hot-tub.
    // A FIFO as a vector and a read position, reset when it drains: it keeps
    // its capacity, where a std::list allocated a node per candidate.
    std::vector<PrimTCandRep*> m_cand_queue;
    int m_cand_queue_head = 0;

    // Per-(di)layer geometrical state
    MkRZLimits m_rz_limits;

    // The in-layer combinatorial tree: one arena per finder, reached by index,
    // cleared with capacity kept at end of layer.
    std::vector<SecTCandRep> m_sec_arena;
    // Slots of nodes whose CombCandidate has been selected, reused before the
    // arena grows. A node is freed only with its whole CombCandidate, so no live
    // node can have a freed parent.
    std::vector<int> m_sec_free;
    // Arena indices appended by the last harvest_sec_nodes(): the next frontier.
    std::vector<int> m_sec_harvested;
    std::vector<int> m_sec_frontier;   // the depth being expanded
    // Running Kalman batches, and the nodes whose children are not yet queued.
    KalmanOpArgs m_koa_d0, m_koa_deep;
    std::vector<int> m_deep_queue;
    // CombCandidates whose last in-flight work finished outside their own layer
    // batch, checked by select_completed_ccreps().
    std::vector<CCandRep*> m_ccrep_ready;
    std::vector<CCandRep*> m_ccrep_touched;
    // Kalman outcomes of the depth currently being expanded, drained into the
    // arena by harvest_sec_nodes().
    std::vector<KalmanOpArgs::ItemOut> m_sec_out;

    // One competitor in the end-of-layer selection. node_idx >= 0 is an in-layer
    // path; otherwise the candidate declined the layer, and add_fake says whether
    // that costs it a HoT (a hole) or nothing at all (the layer was out of
    // reach, or it was already stopped).
    struct SelEntry {
      int   tcand_idx;
      int   node_idx;
      int   fake_hit;
      bool  add_fake;
      float score;
    };
    std::vector<SelEntry> m_sel;
    std::vector<TrackCand> m_new_cands;
  };

} // end namespace mkfit

#endif
