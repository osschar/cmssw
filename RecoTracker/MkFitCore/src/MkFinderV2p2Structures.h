#ifndef RecoTracker_MkFitCore_src_MkFinderV2p2Structures_h
#define RecoTracker_MkFitCore_src_MkFinderV2p2Structures_h

#include "RecoTracker/MkFitCore/interface/HitStructures.h"
#include "RecoTracker/MkFitCore/interface/TrackStructures.h"

#include "MiniPropagators.h"
#include "V2p2Score.h"

#include "V2p2Config.h"
#include <algorithm>
#include <cassert>
#include <new>
#include <type_traits>
#include <vector>

namespace mkfit {

#ifdef MKFIT_TRACE
    class Event;
#endif

    struct CCandRep;
    struct PrimTCandRep;

    //----------------------------------------------------------------------------

    // One node of the in-layer combinatorial tree: a candidate that has taken a
    // particular sequence of hits within this layer, with the Kalman updates
    // applied. Nodes live in MkFinderV2p2::m_sec_arena, are reached by index and
    // link only to their parent. See doc/MkFinderV2p2-DesignNotes.md,
    // "In-layer combinatorial search".
    struct SecTCandRep {
      PrimTCandRep *m_ptc;      // the candidate this path extends
      int   m_parent_idx;       // arena index of the predecessor; -1 = the PrimTCandRep itself
      int   m_hit_pos;          // position of this hit in m_ptc->m_layer_hits; the next
                                // hit a child may take is at m_hit_pos + 1
      HitOnTrack m_hot;
      TrackState m_state;       // AFTER the update at m_hot
      float m_chi2;             // of this hit alone
      // The whole path through this layer, accumulated parent -> child: one layer
      // step however many hits the path took.
      LayerStepFeatures m_feat;
#ifdef MKFIT_TRACE
      int m_tr_hitmatch_id = -1;
#endif
    };

    //----------------------------------------------------------------------------

    // A vector with fixed inline capacity and no construction of the unused
    // slots, for trivially copyable T. Per-candidate, per-layer storage that
    // std::vector would allocate and free on every layer.
    template <typename T, int Cap>
    class InlineVec {
      static_assert(std::is_trivially_copyable_v<T>);
      alignas(T) unsigned char m_buf[Cap * sizeof(T)];
      int m_size = 0;

    public:
      static constexpr int capacity() { return Cap; }
      int size() const { return m_size; }
      bool empty() const { return m_size == 0; }
      void clear() { m_size = 0; }
      T *begin() { return std::launder(reinterpret_cast<T *>(m_buf)); }
      T *end() { return begin() + m_size; }
      const T *begin() const { return std::launder(reinterpret_cast<const T *>(m_buf)); }
      const T *end() const { return begin() + m_size; }
      T &operator[](int i) { return begin()[i]; }
      const T &operator[](int i) const { return begin()[i]; }
      T &front() { return begin()[0]; }
      void push_back(const T &v) {
        assert(m_size < Cap);
        new (m_buf + m_size * sizeof(T)) T(v);
        ++m_size;
      }
      void pop_back() { --m_size; }
    };

    struct PrimTCandRep {
      CCandRep *mp_ccrep;
      int m_origin_tcand_index; // TrackCand index in CombCandidate

      struct PQE { // Priority-Queue Entry
        float score;
        unsigned int hit_orig_index;
        unsigned int hit_index;
        int layer;

#ifdef MKFIT_TRACE
        int tr_hitmatch_id;
#endif
        // float dalpha
        // Hit &
        // ModuleInfo &

        // state on new hit, inv_pt, inv_k, theta on hit on previous layer / last hit
        mini_propagators::InitialState mixed_state;

        // We want the worst / highest score (dphi) at the top -- so we can replace it.
        bool operator<(const PQE& o) const { return score < o.score; }
      };
      // One bounded pqueue per sub-layer, [0] primary and [1] secondary, so each
      // sensor has its own reduction budget. See doc/MkFinderV2p2-DesignNotes.md,
      // "Reduction and hit ordering".
      // Max-heaps on score, worst on top, run with std::push_heap / pop_heap
      // exactly as std::priority_queue<PQE> does, in inline storage.
      static constexpr int kPqCap = Config::V2p2::InLayer::max_presel_hits_limit;
      InlineVec<PQE, kPqCap> m_pqueue[2];
      int m_pqueue_size[2] = {0, 0};

      void pq_push(int sl, const PQE &e) {
        m_pqueue[sl].push_back(e);
        std::push_heap(m_pqueue[sl].begin(), m_pqueue[sl].end());
      }
      const PQE &pq_top(int sl) { return m_pqueue[sl].front(); }
      void pq_pop(int sl) {
        std::pop_heap(m_pqueue[sl].begin(), m_pqueue[sl].end());
        m_pqueue[sl].pop_back();
      }

      // Survivors of both queues, in path order (sorted by dir * dalpha).
      InlineVec<PQE, 2 * kPqCap> m_layer_hits;

      // Did this candidate produce any in-layer path? Set by expand_in_layer().
      bool m_has_sec_nodes = false;
      float m_log_rho = 0.0f;  // local hit density, ln(hits / cm^2), for the score

      PrimTCandRep(CCandRep *ccr, int orig_idx) {
        mp_ccrep = ccr;
        m_origin_tcand_index = orig_idx;
      }

      // Within-sensitive-region verdict for this candidate on this layer, set by
      // MkFinderV2p2::determine_wsr():
      //   WSR_Inside  -- the crossing is inside sensitive q; a missing hit is a hole;
      //   WSR_Edge    -- the crossing is near a boundary, or the track turns round
      //                  inside the layer; a missing hit is not counted;
      //   WSR_Outside -- the track does not reach the layer, which is skipped.
      WSR_Result m_wsr;

      CombCandidate& ccand();
      TrackCand& tcand();

      // these go to SecTCandRep:
      // score, params, chi2, new hits[4 or more for doulbe layers];
      // index in the CombCand -- thing we will update as needed

      // Hack for best hit
      TrackState bState;
      HitOnTrack bHot;
      float bChi2 = 999.999f;
      float bChi2Acc = 999.999f;  // acceptance chi2 of the bHot, Policy::hit_chi2_cut
#ifdef MKFIT_TRACE
      int b_tr_hitmatch_id = -1;
      // Diag::force_mc: the best-hit choice is made on bKey, normally bChi2. With
      // forcing on, an MC-matched hit gets a key below every other hit and bIsMc
      // lets it bypass the chi2 cut; bChi2 keeps the real chi2.
      float bKey = 999.999f;
      bool  bIsMc = false;
#endif
    }; // end struct PrimTCandRep

    //----------------------------------------------------------------------------

    struct CCandRep {
      CombCandidate &m_ccand;

      std::vector<PrimTCandRep> m_primTCs; // for now, could live in the shared arena

      // Arena indices of every in-layer path of this CombCandidate, across all its
      // PrimTCandReps, for the end-of-layer selection.
      std::vector<int> m_sec_nodes;

      // We could also keep track of the TrackCands that do not enter layer
      // processing at all -- either already stopped or missing this layer.
      // Store indices into m_ccand, the way m_primTCs does, rather than
      // pointers: indices stay valid regardless of how the CombCandidate grows.
      // Only worth doing together with SecTCandRep and the selection / merging
      // step -- the best-hit path has no use for it.
      // std::vector<int> m_otherTCs;

      // PrimTCandReps of this layer still in MkFinderV2p2::m_cand_queue. They are
      // queued together, so at zero the CombCandidate is through the layer and
      // can be selected before end of layer.
      int  m_n_queued = 0;
      // With the running Kalman batches: Kalman items in flight plus nodes
      // waiting to be expanded. Selection needs both counts at zero.
      int  m_n_inflight = 0;
      bool m_selected = false;  // select_and_materialise() already ran this layer
      // Beam width this layer, set in begin_next_Ccrep_in_layer(); 0 = the
      // CombCandidate's capacity. See InLayer::late_max_cands.
      int  m_max_cands = 0;

    #if defined(MKFIT_STANDALONE)
      // Tuning & Debugging. Managed in MkFinderV2p2 processing.
      int m_seed_mc_label = -1;
      int m_mc_layer_sequence = -1; // counts layers WITH mc hits
      int m_n_mc_hits_in_layer = -1;
      // int m_n_mc_hits_in_layer_sec = -1;
    #endif

      CCandRep(CombCandidate& ccand) :
        m_ccand(ccand)
      {
         // QQQQ reserve also in begin_next_Ccrep_in_layer()
         // QQQQ clear in end_layer() -- might want to reuse the objects more
        m_primTCs.reserve(ccand.capacity());
      }
    }; // end struct CCandRep

    inline CombCandidate& PrimTCandRep::ccand() { return mp_ccrep->m_ccand; }
    inline TrackCand& PrimTCandRep::tcand() { return mp_ccrep->m_ccand[m_origin_tcand_index]; }

    //----------------------------------------------------------------------------
    // BaseArgs

    struct BaseArgs {

      const PropagationConfig *prop_config = nullptr;

      mini_propagators::InitialStatePlex tsXyz;
      MPlexLS tsErr { 0.0f }; // input (on prev hit) and output (on current hit) [ts - track-state]
      MPlexLV tsPar { 0.0f }; // ""
      MPlexQI tsChg { 0 };    // "" Kalman update can flip it through curvature flip

      MPlexQF sPerp { 0.0f }; // path-length in transverse plane, calculated from alpha. p2plane really needs 3D path.

      MPlexLS propErr { 0.0f }; // intermediate: propagated from tsErr and used as input to Kalman
      MPlexLV propPar { 0.0f }; // input: pre-propagated as part of hit pre-selection

      MPlexQI outFailFlag { 0 }; // dummy, can be detected in pre-propagation, no other errors detected / reported

      int N_filled = 0;
    };

    //----------------------------------------------------------------------------
    // KalmanOpArgs

    struct KalmanOpArgs : public BaseArgs {

      PrimTCandRep *ptcp[NN];
      HitOnTrack    hot[NN];
      // Provenance of each lane, for the in-layer expansion: which arena node this
      // extension grew from (-1 = the PrimTCandRep's own incoming state) and which
      // entry of m_layer_hits it takes. Unused by the best-hit path.
      int           parent_idx[NN];
      int           hit_pos[NN];

      MPlexHS msErr { 0.0f }; // input measurement / hit [ms - measurement state]
      MPlexHV msPar { 0.0f }; // ""
      MPlexHV plNrm { 0.0f }; // input detector plane [pl - plane]
      MPlexHV plDir { 0.0f }; // ""
      MPlexHV plPnt { 0.0f }; // ""

      MPlexQF tsChi2 { 0.0f };   // output
      MPlexQF tsDetV { 0.0f };   // output: det of the 2x2 residual covariance
      MPlexQF tsChi2Acc { 0.0f };  // output: acceptance chi2, Policy::chi2_trk_fac

#ifdef MKFIT_TRACE
      int tr_hitmatch_ids[NN];
      void set_tr_hitmatch_id(int id) { tr_hitmatch_ids[N_filled] = id; }
      const Event *mp_event = nullptr;
#endif

      // One result per lane, appended when mp_out is set: the in-layer search keeps
      // every outcome, the best-hit path only the winner (PrimTCandRep::b*).
      struct ItemOut {
        PrimTCandRep *ptc;
        int   parent_idx;
        int   hit_pos;
        HitOnTrack hot;
        unsigned int hit_in_layer;   // index WITHIN the LayerOfHits, for hit_q_half_length()
        float chi2;
        float det_v;                 // det of the 2x2 residual covariance, for the score
        float chi2_acc;              // acceptance chi2, compared with Policy::hit_chi2_cut
        TrackState state;
#ifdef MKFIT_TRACE
        int   tr_hitmatch_id;
#endif
      };
      std::vector<ItemOut> *mp_out = nullptr;

      // Solve for the crossing with the module plane instead of taking sPerp and
      // propPar from the search. Set for depth >= 1 of the in-layer search, where
      // the Hermite crossings no longer apply.
      bool m_solve_plane = false;

      void reset() { N_filled = 0; }

      // There will be some more of this state, also secondary or who knows what.
      unsigned int hit_in_layer[NN];

      void item_begin(PrimTCandRep *ptc, HitOnTrack ht, int par_idx = -1, int hpos = -1,
                      unsigned int hil = 0) {
        ptcp[N_filled] = ptc; hot[N_filled] = ht;
        parent_idx[N_filled] = par_idx; hit_pos[N_filled] = hpos;
        hit_in_layer[N_filled] = hil;
      }
      bool item_finished() { return ++N_filled == NN; }

      // Start from a TrackState directly, for the solve-plane path: there is no
      // pre-solved crossing, so tsXyz is not filled and compute_pars() must be
      // skipped -- propPar is a pure output there, where on the sPerp path it is
      // an input carrying the already-known intersection.
      void load_state_err_chg(const TrackState &st) {
        tsPar.copyIn(N_filled, st.parArray());
        tsErr.copyIn(N_filled, st.errArray());
        tsChg[N_filled] = st.charge;
      }

      void load_state_err_chg(const mini_propagators::InitialState &params_on_hit, const TrackState &prev_state) {
        tsXyz.copyIn(N_filled, params_on_hit);
        tsPar.copyIn(N_filled, prev_state.parArray()); // propToPlane needs initial parameters, too
        tsErr.copyIn(N_filled, prev_state.errArray());
        tsChg[N_filled] = prev_state.charge;
      }

      void load_hit_module(const Hit &hit, const ModuleInfo & mi) {
        msErr.copyIn(N_filled, hit.errArray());
        msPar.copyIn(N_filled, hit.posArray());
        plNrm.copyIn(N_filled, mi.zdir.Array());
        plDir.copyIn(N_filled, mi.xdir.Array());
        plPnt.copyIn(N_filled, mi.pos.Array());
      }

      void compute_pars() {
        // Parameters are stored in StatePlex -- so we can vectorize translation to pt, phi, theta.
        // Some stuff could be passed over as it won't change before update: pt, theta, k_inv
        // They are passed in output-parameters as propagation also needs input pars.
        propPar.aij(0, 0) = tsXyz.x;
        propPar.aij(1, 0) = tsXyz.y;
        propPar.aij(2, 0) = tsXyz.z;
        propPar.aij(3, 0) = tsXyz.inv_pt;
        propPar.aij(4, 0) = Matriplex::fast_atan2(tsXyz.py, tsXyz.px);
        propPar.aij(5, 0) = tsXyz.theta;

        sPerp = tsXyz.dalpha / ( tsXyz.inv_pt * tsXyz.inv_k);
      }

      void do_kalman_stuff();
    };

}

#endif
