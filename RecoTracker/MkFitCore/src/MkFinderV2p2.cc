#include "MkFinderV2p2.h"

#include "RecoTracker/MkFitCore/interface/cms_common_macros.h"
// #include "FindingFoos.h"

#include "MkBins.h"

#include <algorithm>

#if defined(MKFIT_STANDALONE)
#include "RecoTracker/MkFitCore/standalone/Event.h"
#endif

#ifdef MKFIT_TRACE
#include "RecoTracker/MkFitCore/standalone/DataFormats/RntStructs.h"
#include "RecoTracker/MkFitCore/standalone/DataFormats/RntConversions.h"
#endif

//#define DEBUG
#include "Debug.h"

namespace mkfit {

  using namespace Config::V2p2;


  //------------------------------------------------------------------------------
  // Setup variables for full processing of a batch of CombCanditates

  void MkFinderV2p2::setup(const MkJob *job, EventOfCombCandidates &eoccs, int seed_begin, int seed_end,
                           SteeringParams::iterator &sp_it, const Event *ev) {
    mp_job = job;
    mp_steeringparams_iter = &sp_it;
    mp_event = ev;

    m_batch_mgr.setup(eoccs, seed_begin, seed_end);

    // Assure all ccands are dormant and have a single tcand (?)
    for (auto &ccand : m_batch_mgr) {
      if (ccand.state() != CombCandidate::Dormant)
        throw std::runtime_error("CombCandidate is expected to be dormant "
                                 "on entry into combinatorial search.");
      if (ccand.size() != 1)
        throw std::runtime_error("CombCandidate is expected to have a single "
                                 "TrackCand on entry into combinatorial search.");
    }
    m_batch_mgr.m_n_dormant = m_batch_mgr.n_total();

    // AAAA reserve based on iter-param max-cands and N-sub-layer / skipped hit combinatorials.
  }

  void MkFinderV2p2::release() {
    m_batch_mgr.release();
    mp_job = nullptr;
    mp_steeringparams_iter = nullptr;
    mp_event = nullptr;

    // XXXX this might need to be achieved through other means ... crashing when going to
    // XXXX multi-event setup as old ccands have stayed in.
    // See another XXXX about finalization below in end_layer()
    m_active_ccreps.clear();

    // Finders are reused across seed blocks and events, and the queue holds
    // pointers into CCandRep::m_primTCs. process_layer() asserts it is empty on exit.
    m_cand_queue.clear();
    m_cand_queue_head = 0;
  }

  //------------------------------------------------------------------------------
  // Helpers.

  int MkFinderV2p2::awaken_candidates() {
    // Awaken dormant ccands with hits on current layer for processing on the next layer.
    int count = 0;
    const LayerControl &lc = mp_steeringparams_iter->layer_control();
    for (auto &ccand : m_batch_mgr) {

      // Pick up at either sub-layer of a paired entry. A seed whose last hit is
      // in this layer needs no special case. See doc/MkFinderV2p2-DesignNotes.md,
      // "Candidate pickup and stopping cuts".

      if (ccand.state() == CombCandidate::Dormant &&
          (ccand.pickupLayer() == lc.m_layer || ccand.pickupLayer() == lc.m_layer_sec)) {
        ccand.setState(CombCandidate::Finding);
        // Start the in-flight score accumulator. Any constant seed term is
        // common to every candidate of this seed and so cancels in the per-seed
        // selection; the cross-seed comparison is done at the end, on the
        // official score.
        if (InLayer::comb)
          for (int ic = 0; ic < (int) ccand.size(); ++ic)
            ccand[ic].setScore(0.0f);
        // auto &: the MKFIT_STANDALONE fields below must land on the list element.
        [[maybe_unused]] auto &ccrep = m_active_ccreps.emplace_back(ccand);
        dprintf("MkFinderV2p2::awaken_candidates dummy printout N_TrackCands=%d\n",
               (int) ccand.size());
        ++count;

        #if defined(MKFIT_STANDALONE)
        // XXXX hmmh, i could really just get it where i need it.
        auto sifh = mp_event->simInfoForCurrentSeed(ccand.seed_origin_index());
        ccrep.m_seed_mc_label = sifh.label;
        ccrep.m_mc_layer_sequence = 0;
        // ccrep.m_n_mc_hits_in_current_layer set later
        #endif
      }
    }
    m_batch_mgr.m_n_dormant -= count;
    m_batch_mgr.m_n_finding += count;

    dprintf("MkFinderV2p2::awaken_candidates woke up %d cands\n", count);

    return count;
  }

  //------------------------------------------------------------------------------
  // stop_cuts_at_pickup() -- minPtCut and the looper stop, which need only the
  // candidate's own state. The looper test covers both wrap images of
  // |posPhi - momPhi|. See doc/MkFinderV2p2-DesignNotes.md,
  // "Candidate pickup and stopping cuts".
  //----------------------------------------------------------------------------

  TrackCand::StopReason_e MkFinderV2p2::stop_cuts_at_pickup(const TrackCand &tc) const {
    const SteeringParams::iterator &spi = *mp_steeringparams_iter;
    const auto &iter_params = (spi.type() == SteeringParams::IT_BkwSearch) ? mp_job->params_bks()
                                                                          : mp_job->params();

    if (tc.pT() < iter_params.minPtCut)
      return TrackCand::SR_MinPt;

    // Forward search only: going inward the track is leaving the turning region,
    // so a large angle there is not evidence that it is about to curl up.
    if (spi.type() == SteeringParams::IT_FwdSearch && tc.pT() < Policy::looper_max_pt &&
        tc.posRsq() > Policy::looper_min_r * Policy::looper_min_r) {
      const float max_ang = Policy::looper_max_angle;
      const float dphi = std::abs(tc.posPhi() - tc.momPhi());
      if (dphi > max_ang && dphi < Const::TwoPI - max_ang)
        return TrackCand::SR_Looper;
    }

    return TrackCand::SR_NotStopped;
  }

  //----------------------------------------------------------------------------
  // fake_hit_index() -- the fake HoT for a candidate that took no hit: the hole
  // limits choose miss or stop, then the WSR overrides with edge or gap, in
  // V1's order.
  //----------------------------------------------------------------------------

  int MkFinderV2p2::fake_hit_index(const TrackCand &tc, const WSR_Result &wsr) const {
    const SteeringParams::iterator &spi = *mp_steeringparams_iter;
    const auto &iter_params = (spi.type() == SteeringParams::IT_BkwSearch) ? mp_job->params_bks()
                                                                          : mp_job->params();

    int fake_hit_idx = Hit::kHitMissIdx;

    if (Policy::use_hole_limits &&
        (tc.nAllMinusOneHits() >= iter_params.maxHolesPerCand ||
         tc.nTailMinusOneHits() >= iter_params.maxConsecHoles))
      fake_hit_idx = Hit::kHitStopIdx;

    if (Policy::use_wsr) {
      if (wsr.m_wsr == WSR_Edge)
        fake_hit_idx = Hit::kHitEdgeIdx;
      else if (wsr.m_in_gap)
        fake_hit_idx = Hit::kHitInGapIdx;
    }

    return fake_hit_idx;
  }

  //------------------------------------------------------------------------------
  // Per layer initialization / cleanup tasks,

  void MkFinderV2p2::begin_layer() {
    // debug = true; // to be disabled at the end of end_layer()

    // Count number of awakend ccands and number of non-stopped tcands.
    // Well, I actually know n ccands.

    dprintf("MkFinderV2p2::begin_layer Expecting %d active ccands\n", (int) m_active_ccreps.size());

// #ifdef RNT_DUMP_MkF_SelHitIdcs
//     // clang-format off
//     const IterationConfig &IC = mp_job->m_iter_config;
//     SteeringParams::iterator &spi = *mp_steeringparams_iter;
//     const LayerOfHits &L = mp_job->m_event_of_hits[spi->m_layer];
//     const LayerInfo &LI = L.layer_info();
//     rnt_shi.ResetH();
//     rnt_shi.ResetF();
//     *rnt_shi.h = {mp_event->evtID(), IC.m_iteration_index, IC.m_track_algorithm,
//                   spi.region(), L.layer_id(),
//                   L.is_barrel() ? LI.rin() : LI.zmin(), LI.is_barrel() ? LI.rout() : LI.zmax(),
//                   L.is_barrel(), L.is_pixel(), L.is_stereo()};
//     *rnt_shi.f = *rnt_shi.h;

//     // Unlike in MkFinder, we process all active cands in one go.
//     // And actually do not use the inner indices as we assign CandInfo* into PrimTCand.
//     rnt_shi.InnerIdcsReset((int) m_active_ccreps.size());
//     // clang-format on
// #endif

#ifdef DEBUG
    { int i = 1;
      for (auto &ccrep : m_active_ccreps) {
        dprintf("  %2d. seed-idx=%d, n_tcands=%d\n", i,
               ccrep.m_ccand.seed_origin_index(), (int) ccrep.m_ccand.size());
        ++i;
      }
    }
#endif

    m_batch_mgr.reset_for_new_layer();
    m_active_ccreps_pos = m_active_ccreps.begin();

    { // Setup m_rz_limits.
      SteeringParams::iterator &spi = *mp_steeringparams_iter;

      const bool is_double_layer = spi->has_second_layer();
      const LayerInfo &LI_p = mp_job->m_trk_info[ spi->m_layer ];
      const LayerInfo &LI_s = mp_job->m_trk_info[ is_double_layer ? spi->m_layer_sec : 0 ];

      if (is_double_layer) {
        assert(LI_p.is_barrel() == LI_s.is_barrel() );
        m_rz_limits.setup(LI_p, LI_s, spi.is_outward());
      } else {
        m_rz_limits.setup(LI_p, spi.is_outward());
      }
    }

    // Hmmh, nothing to really do here, is it?

    // Initialize best short (?) - to what? How was it before? The seed itself?
    // Or let this happen after the first layer is processed?
    // Score is set to worst-possible elsewhere.
  }

  // Straight-line path length from the candidate's state (its last hit) to the
  // layer's near surface, against InLayer::late_max_cands_step_cm. Barrel: radial gap
  // times p / pT; endcap: z gap times p / |pz|. Curvature is ignored, so the
  // barrel path is underestimated for strongly bending tracks.
  // InLayer::sister_hole_drop: remove from m_sel the holes that a path of the same
  // TrackCand dominates by taking the sister sensor of its last hit. A sister pair
  // is two adjacent layers, the lower one stereo and the upper one not, and the
  // upper sensor's detid is the lower one's + 1; either may be crossed first.
  void MkFinderV2p2::drop_sister_dominated_holes(const CombCandidate &ccand) {
    const int lay = m_rz_limits.layer_info_1().layer_id();
    const TrackerInfo &ti = mp_job->m_trk_info;
    auto detid_of = [this](int layer, int idx) {
      const LayerOfHits &L = mp_job->m_event_of_hits[layer];
      return L.layer_info().module_info(L.refHit(idx).detIDinLayer()).detid;
    };
    int w = 0;
    for (int k = 0; k < (int) m_sel.size(); ++k) {
      const SelEntry e = m_sel[k];
      bool drop = false;
      if (e.node_idx < 0 && e.add_fake) {
        const TrackCand &tc = ccand[e.tcand_idx];
        const int pl = tc.getLastHitLyr(), pi = tc.getLastHitIdx();
        const int lo = std::min(pl, lay), hi = std::max(pl, lay);
        if (pi >= 0 && hi == lo + 1 && ti.layer(lo).is_stereo() && !ti.layer(hi).is_stereo()) {
          const unsigned int d_prev = detid_of(pl, pi);
          for (const SelEntry &f : m_sel) {
            if (f.node_idx < 0 || f.tcand_idx != e.tcand_idx)
              continue;
            int root = f.node_idx;
            while (m_sec_arena[root].m_parent_idx >= 0)
              root = m_sec_arena[root].m_parent_idx;
            const SecTCandRep &r = m_sec_arena[root];
            if (InLayer::sister_hole_chi2 > 0.0f && !(r.m_chi2 < InLayer::sister_hole_chi2))
              continue;
            const unsigned int d_this = detid_of(r.m_hot.layer, r.m_hot.index);
            if (pl < lay ? d_this == d_prev + 1 : d_prev == d_this + 1) {
              drop = true;
              break;
            }
          }
        }
      }
      if (drop)
        V2P2_COUNT(n_sister_hole_dropped);
      else
        m_sel[w++] = e;
    }
    m_sel.resize(w);
  }


  bool MkFinderV2p2::long_step_to_layer(const TrackCand &tc) const {
    if (InLayer::late_max_cands_step_cm <= 0.0f)
      return false;
    const LayerInfo &li = m_rz_limits.layer_info_1();
    float s;
    if (li.is_barrel()) {
      const float r = tc.posR();
      const float dr = std::max({li.rin() - r, r - li.rout(), 0.0f});
      s = dr * tc.p() / std::max(1e-3f, tc.pT());
    } else {
      const float z = tc.z();
      const float dz = std::max({li.zmin() - z, z - li.zmax(), 0.0f});
      s = dz * tc.p() / std::max(1e-3f, std::abs(tc.pz()));
    }
    if (s > InLayer::late_max_cands_step_cm) {
      V2P2_COUNT(n_late_max_cands_long);
      return true;
    }
    return false;
  }

  void MkFinderV2p2::begin_next_Ccrep_in_layer() {
    CCandRep &ccrep = * m_active_ccreps_pos;
    CombCandidate &ccand = ccrep.m_ccand;

    // QQQQ reserve to CombCand.capacity already done CCandRep ctor.
    // We will probably need a variable number per layer / layer pair.
    // But this will be first relevant / handled elsewhere.
    // Beam width for this layer, per seed: see InLayer::late_max_cands. The
    // TrackCands are in score order from the previous selection, so dropping
    // the tail here is the selection having kept fewer.
    // The capacity is reserved for the larger of the forward and backward widths
    // (MkJob::max_max_cands()); the search direction's own maxCandsPerSeed applies.
    const auto &iter_params = (mp_steeringparams_iter->type() == SteeringParams::IT_BkwSearch)
                                  ? mp_job->params_bks()
                                  : mp_job->params();
    ccrep.m_max_cands = std::min((int) ccand.capacity(), iter_params.maxCandsPerSeed);
    if (InLayer::late_max_cands > 0 && ccand.size() > 0 && ccand[0].nFoundHits() > InLayer::late_max_cands_hits &&
        !long_step_to_layer(ccand[0])) {
      ccrep.m_max_cands = std::min(ccrep.m_max_cands, (int) InLayer::late_max_cands);
      V2P2_COUNT(n_late_max_cands);
      if ((int) ccand.size() > ccrep.m_max_cands) {
        V2P2_COUNT_ADD(n_late_max_cands_dropped, (int) ccand.size() - ccrep.m_max_cands);
        ccand.resize(ccrep.m_max_cands);
      }
    }

    ccrep.m_primTCs.reserve(ccand.size());

    for (int ic = 0; ic < (int) ccand.size(); ++ic) {
      TrackCand &tcand = ccand[ic];

      if (tcand.getLastHitIdx() == Hit::kHitStopIdx)
        continue;

      // Stopping cuts, at pull-in: see stop_cuts_at_pickup().
      if (Policy::use_stop_cuts) {
        const TrackCand::StopReason_e sr = stop_cuts_at_pickup(tcand);
        if (sr != TrackCand::SR_NotStopped) {
          tcand.setStopReason(sr);
          tcand.addHitIdx(Hit::kHitStopIdx, m_rz_limits.layer_info_1().layer_id(), 0.0f);
          if (sr == TrackCand::SR_MinPt)
            V2P2_COUNT(n_stop_minpt);
          else
            V2P2_COUNT(n_stop_looper);
          continue;
        }
      }

      // XXXX Should one do rough "layer already passed" pre-check here?
      // Or in pre-select, where we engage MkBins ... but there I loose a vector slot.
      // Let's try.
      //
      // A failed check skips the layer without recording anything, as WSR_Outside
      // does, with a coarser test.

      if ( ! m_rz_limits.rz_quadrant_check(tcand.z(), tcand.pz())) {
        V2P2_COUNT(n_quadrant_skip);
        continue;
      }

      // Create and Register PrimTCandRep for processing.
      // The CCandRep vector has the capacity for N_max_cands.
      {
        PrimTCandRep &ptc = ccrep.m_primTCs.emplace_back( &ccrep, ic );
        m_cand_queue.push_back(&ptc);
      }
    }
    ccrep.m_n_queued = (int) ccrep.m_primTCs.size();
    ccrep.m_selected = false;
    ++m_active_ccreps_pos;

    #if defined(MKFIT_STANDALONE)
    {
      ccrep.m_n_mc_hits_in_layer = mp_event->countSimHitsInLayer(ccrep.m_seed_mc_label, mp_steeringparams_iter->layer());
      // ccrep.m_n_mc_hits_in_layer_sec = mp_event->countSimHitsInLayer(ccrep.m_seed_mc_label, mp_steeringparams_iter->layer_sec());
      if (ccrep.m_n_mc_hits_in_layer > 0)
        ++ccrep.m_mc_layer_sequence;
    }
    #endif
  }

  // void MkFinderV2p2::process_layer_batch() -- below in the "complex stuff" section

  void MkFinderV2p2::end_layer() {
    // Stop tracks -- pT / apogee / missing layers.
    // Choose best-short.
    // Figure out what to copy back to EventOfCombCandidates.

    // clear out ccands -- well, might keep them -- just flush the tcands out of hot-tub
    int count = 0;
    auto ai = m_active_ccreps.begin();
    while (ai != m_active_ccreps.end()) {

      if (InLayer::comb && ! ai->m_selected)
        select_and_materialise(*ai);

      // QQQQ should attempt to reuse the PrimTCandReps
      ai->m_primTCs.clear();
      ai->m_sec_nodes.clear();

      // A CombCandidate is finished once every TrackCand under it has stopped.
      bool is_finished = false;
      if (Policy::use_stop_cuts || Policy::use_hole_limits) {
        CombCandidate &cc = ai->m_ccand;
        is_finished = true;
        for (int ic = 0; ic < (int) cc.size(); ++ic) {
          if (cc[ic].getLastHitIdx() != Hit::kHitStopIdx) {
            is_finished = false;
            break;
          }
        }
        if (is_finished) {
          cc.setState(CombCandidate::Finished);
          V2P2_COUNT(n_ccand_retired);
        }
      }
      if (is_finished) {
        auto bi = ai++;
        m_active_ccreps.erase(bi);
        ++count;
        continue;
      }
      ++ai;
    }
    m_batch_mgr.m_n_finding -= count;
    m_batch_mgr.m_n_finished += count;

    // One rewind for the whole layer, keeping capacity: after a few layers the
    // arena is at high water and stops allocating. Its size here is this layer's
    // high water, since slots are reused but the vector never shrinks.
#if defined(MKFIT_STANDALONE)
    if (InLayer::comb) {
      const long hw = (long) m_sec_arena.size();
      V2P2_COUNT(n_arena_layers);
      V2P2_COUNT_ADD(n_arena_hw_sum, hw);
      long mx = g_v2p2_policy_counters.n_arena_hw_max.load();
      while (hw > mx && ! g_v2p2_policy_counters.n_arena_hw_max.compare_exchange_weak(mx, hw)) {}
    }
#endif
    m_sec_arena.clear();
    m_sec_free.clear();

    m_rz_limits.reset();

    dprintf("MkFinderV2p2::end_layer %d cands finished\n", count);

    if (m_batch_mgr.has_dormant_ccands())
      awaken_candidates();

    // debug = false;
  }

  //------------------------------------------------------------------------------
  // The main processing function -- process_layer()

  void MkFinderV2p2::process_layer() {

    // Pull CombCandidates into the layer until NN candidates are queued, then
    // drain the queue NN at a time. See doc/MkFinderV2p2-DesignNotes.md,
    // "Layer processing".
    while (any_Ccreps_to_begin()) {
      while ( ! enough_work_for_batch() && any_Ccreps_to_begin()) {
        begin_next_Ccrep_in_layer();
      }
      while (enough_work_for_batch() ||
             ( ! any_Ccreps_to_begin() && any_work_for_batch())) {
        process_layer_batch();
      }
    }

    // Drained by construction: the inner loop's second clause runs until empty
    // once there is nothing left to pull in.
    assert(cand_queue_size() == 0 && "pre-select queue not drained by process_layer()");

    // The running Kalman batches hold partial work across layer batches; empty
    // them before end_layer() selects whatever is left.
    if (InLayer::comb && InLayer::running_kalman)
      rk_drain();
  }

  //============================================================================
  // More complex functions -- to separate them from the main "logic" flow
  //============================================================================

  //----------------------------------------------------------------------------
  // process_layer_batch() -- one batch of up to NN candidates through the layer.
  // LayerBatch is NN candidates wide, HitBatch is NN (candidate, hit) pairs wide.
  // See doc/MkFinderV2p2-DesignNotes.md, "Layer processing".
  //----------------------------------------------------------------------------

  void MkFinderV2p2::process_layer_batch() {
    LayerBatch b;

    prop_to_layer_edges(b);       // stage 1: pop the queue, propagate to the layer edges
    determine_search_windows(b);  // covariance -> dphi/dq windows -> bin ranges
    select_hits(b);               // stage 2: walk the bins, pre-select, reduce in a pqueue
    prepare_kalman_workload(b);   // pqueue -> m_layer_hits, step-ordered
    if (InLayer::comb) {
      // No materialisation here: the paths are left in the arena and everything
      // competes at end of layer, in select_and_materialise().
      if (InLayer::running_kalman)
        expand_in_layer_running(b);
      else
        expand_in_layer(b);       // stage 3: grow the SecTCandRep tree over the ordered hits
      if (InLayer::early_select)
        select_completed_ccreps(b);
    } else {
      kalman_update(b);           // stage 3: propagate to each module plane + update
      process_kalman_results(b);  // best-hit acceptance into the TrackCand
    }
  }

  //----------------------------------------------------------------------------
  // Phase 1 -- which candidates, and where do they meet the layer.
  //----------------------------------------------------------------------------

  void MkFinderV2p2::prop_to_layer_edges(LayerBatch &b) {
    MkBins &B = b.B;
    PrimTCandRep **prim_tcand_ptrs = b.ptc;

    const int N_proc = b.N_proc = std::min(NN, cand_queue_size());
    B.m_n_proc = N_proc;   // MkBins used to be constructed with it

    dprintf("MkFinderV2p2::process_layer_batch work queue is %d, would process %d of them (NN=%d)\n",
            cand_queue_size(), N_proc, NN);

    MPlexQF phi(0.0f);
    MPlexQI chg(0);
    MPlexQF bf(Config::Bfield);
    const bool param_b = Window::mini_prop_param_b == 1 ||
                         (Window::mini_prop_param_b == 2 &&
                          mp_job->m_trk_info.prop_config().finding_inter_layer_pflags.use_param_b_field);

    for (int i = 0; i < N_proc; ++i) {
      PrimTCandRep &ptc = * m_cand_queue[m_cand_queue_head++];
      prim_tcand_ptrs[i] = & ptc;
      TrackCand &tc = ptc.tcand();
      --ptc.mp_ccrep->m_n_queued;

      // Copy in x, y,z, invpT, theta.
      B.m_isp.copyIn_partial_track_state(i, tc.state());
      phi[i] = tc.momPhi();
      chg[i] = tc.charge();
      if (param_b)
        bf[i] = Config::bFieldFromZR(tc.z(), std::hypot(tc.x(), tc.y()));
    }
    B.m_isp.init_momentum_vec_and_k(phi, chg, bf);
    if (m_cand_queue_head == (int) m_cand_queue.size()) {
      m_cand_queue.clear();
      m_cand_queue_head = 0;
    }

    // Propagation so point 1 is first edge hit, 2 the second
    B.prop_to_limits_in_order(m_rz_limits);

  }

  //----------------------------------------------------------------------------
  // Phase 2 -- the search windows. The position block of the track covariance is
  // transported to m_sp2 (MkBins::transport_position_cov()); from it come
  // dphi_track and dq_track, the WSR verdict, the binnor ranges and the Hermite
  // cubic. See
  // doc/MkFinderV2p2-DesignNotes.md, "Search window".
  //----------------------------------------------------------------------------

  void MkFinderV2p2::determine_search_windows(LayerBatch &b) {
    SteeringParams::iterator &spi = *mp_steeringparams_iter;
    const int N_proc = b.N_proc;
    PrimTCandRep **prim_tcand_ptrs = b.ptc;
    MkBins &B = b.B;
    MkBinTrackCovExtract &TCE = b.TCE;
    MkBinLimits &BL_p = b.BL_p;
    MkBinLimits &BL_s = b.BL_s;
    mini_propagators::Hermite3D &H = b.H;
#ifdef MKFIT_TRACE
    int *tr_layersearch_ids = b.tr_layersearch_ids;
#endif

    // Covariance position block at m_sp2, transported from the previous hit.
    MPlexLV par0;
    MPlexLS err0;
    for (int i = 0; i < N_proc; ++i) {
      const TrackCand &tc = prim_tcand_ptrs[i]->tcand();
      par0.copyIn(i, tc.posArray());
      err0.copyIn(i, tc.errArray());
    }
    for (int i = N_proc; i < NN; ++i) {
      par0.copyIn(i, par0, 0);
      err0.copyIn(i, err0, 0);
    }
    B.transport_position_cov(par0, err0, TCE);

    B.determine_bin_windows(TCE);

    // The WSR edge fuzz uses dq_track, so the verdict follows determine_bin_windows().
    determine_wsr(b);

    B.find_bin_ranges(mp_job->m_event_of_hits[spi->m_layer], BL_p);

    if (m_rz_limits.m_is_double) {
      B.find_bin_ranges(mp_job->m_event_of_hits[spi->m_layer_sec], BL_s);
    }

    for (int i = 0; i < N_proc; ++i) {
      dprintf("%d: BinCheck Prim %c %+8.6f %+8.6f | %3d %3d || %+8.6f %+8.6f | %2d %2d\n",
              i, m_rz_limits.m_is_barrel ? 'B' : 'E',
              B.m_phi_center[i], B.m_phi_delta[i], BL_p.p1[i], BL_p.p2[i],
              B.m_q_min[i], B.m_q_max[i], BL_p.q1[i], BL_p.q2[i]);
      if (m_rz_limits.m_is_double) {
        dprintf("%d: BinCheck Sec  %c %+8.6f %+8.6f | %3d %3d || %+8.6f %+8.6f | %2d %2d\n",
                i, m_rz_limits.m_is_barrel ? 'B' : 'E',
                B.m_phi_center[i], B.m_phi_delta[i], BL_s.p1[i], BL_s.p2[i],
                B.m_q_min[i], B.m_q_max[i], BL_s.q1[i], BL_s.q2[i]);
      }
    }

    // Hermite cubic through the two layer crossings, the trajectory model for the
    // per-hit plane solve. See the design notes, "Layer crossings and the Hermite
    // cubic".
    // This might belong better somewhere else, TPrimCanRep? Calculation into minipropagators.

    namespace mp = mini_propagators;

    H.calculate_coeffs(B.m_sp1, B.m_sp2, B.m_isp.inv_k);

#ifdef MKFIT_TRACE
    for (int i = 0; i < N_proc; ++i) {
      TrLayerSearch ls;
      ls.state_id   = prim_tcand_ptrs[i]->tcand().m_trace_state_id;
      ls.layer      = spi->m_layer;
      ls.layer_sec  = m_rz_limits.m_is_double ? spi->m_layer_sec : -1;
      ls.is_barrel  = m_rz_limits.m_is_barrel;
      ls.is_outward = m_rz_limits.m_is_outward;
      ls.wsr        = (signed char) prim_tcand_ptrs[i]->m_wsr.m_wsr;
      ls.wsr_in_gap = prim_tcand_ptrs[i]->m_wsr.m_in_gap;

      ls.prop_entry = statep2propinfo(B.m_sp1, i);
      ls.prop_exit  = statep2propinfo(B.m_sp2, i);

      ls.phi_center = B.m_phi_center[i];
      ls.phi_delta  = B.m_phi_delta[i];
      ls.q_center   = B.m_q_center[i];
      ls.q_min      = B.m_q_min[i];
      ls.q_max      = B.m_q_max[i];
      ls.dphi_track = B.m_dphi_track[i];
      ls.dq_track   = B.m_dq_track[i];

      ls.cov_xx = TCE.m_cov_0_0[i];
      ls.cov_xy = TCE.m_cov_0_1[i];
      ls.cov_yy = TCE.m_cov_1_1[i];
      ls.cov_zz = TCE.m_cov_2_2[i];

      ls.p1 = BL_p.p1[i];  ls.p2 = BL_p.p2[i];
      ls.q1 = BL_p.q1[i];  ls.q2 = BL_p.q2[i];
      if (m_rz_limits.m_is_double) {
        ls.p1_sec = BL_s.p1[i];  ls.p2_sec = BL_s.p2[i];
        ls.q1_sec = BL_s.q1[i];  ls.q2_sec = BL_s.q2[i];
      }

      tr_layersearch_ids[i] = mp_event->trace_layersearch(std::move(ls)).id;
    }
#endif

#ifdef MKFIT_TRACE_PROP_COMPARE
    // Hermite vs helix in mid-layer, where the cubic is furthest from the two
    // endpoints it interpolates (at t = 0 and 1 it reproduces them by
    // construction, so those are roundoff checks only and not worth taking).
    // Evaluate the cubic at t = 0.5, then propagate exactly to the bounding
    // surface -- r in the barrel, z in the endcap -- that the cubic's midpoint
    // landed on, so both points sit on the same surface and the difference is a
    // pure in-surface deviation.
    {
      MPlexQF hx, hy, hz;
      H.evaluate(0.5f, hx, hy, hz);
      mp::StatePlex sp_mid;
      if (m_rz_limits.m_is_barrel) {
        MPlexQF r_mid = Matriplex::hypot(hx, hy);
        B.m_isp.propagate_to_r(mp::PA_Exact, r_mid, sp_mid, false, N_proc);
      } else {
        B.m_isp.propagate_to_z(mp::PA_Exact, hz, sp_mid, false, N_proc);
      }
      for (int i = 0; i < N_proc; ++i) {
        EVec3 dev(hx[i] - sp_mid.x[i], hy[i] - sp_mid.y[i], hz[i] - sp_mid.z[i]);
        dprintf("%d: HermiteMid dev %+9.6f %+9.6f %+9.6f  (|d| = %.6f, dalpha %.4f -> %.4f, derfac %.5f)\n",
                i, dev[0], dev[1], dev[2], std::sqrt(dev.Mag2()),
                B.m_sp1.dalpha[i], B.m_sp2.dalpha[i], H.m_Hderfac[i]);
#ifdef MKFIT_TRACE
        mp_event->tr_layersearch(tr_layersearch_ids[i]).hermite_mid_dev = dev;
#endif
      }
    }
#endif

  }

  //----------------------------------------------------------------------------
  // Phase 2b -- within-sensitive-region verdict per candidate: does the segment
  // between the two layer crossings lie in sensitive material. The fail flags
  // of the two crossings answer the radial half; the segment's q extent,
  // widened by Policy::wsr_n_sigma, answers the rest. See
  // doc/MkFinderV2p2-DesignNotes.md, "Within-sensitive-region verdict".
  //----------------------------------------------------------------------------

  void MkFinderV2p2::determine_wsr(LayerBatch &b) {
    // m_dq_track is 3 sigma.
    const float dq_to_n_sigma = Policy::wsr_n_sigma / 3.0f;

    const int N_proc = b.N_proc;
    MkBins &B = b.B;
    PrimTCandRep **prim_tcand_ptrs = b.ptc;

    const bool is_barrel = m_rz_limits.m_is_barrel;
    const float q_lo = is_barrel ? m_rz_limits.m_zmin : m_rz_limits.m_rin;
    const float q_hi = is_barrel ? m_rz_limits.m_zmax : m_rz_limits.m_rout;
    const LayerInfo &LI = m_rz_limits.layer_info_1();

    for (int i = 0; i < N_proc; ++i) {
      const bool f_entry = B.m_sp1.fail_flag[i] != 0;
      const bool f_exit  = B.m_sp2.fail_flag[i] != 0;

      WSR_Result w;

      if (f_entry && f_exit) {
        // Clear miss: the track does not reach the layer at all.
        w = WSR_Result(WSR_Outside, false);
      } else {
        const float dq = dq_to_n_sigma * B.m_dq_track[i];
        const float qa = B.m_q_min[i], qb = B.m_q_max[i];

        if (qb < q_lo - dq || qa > q_hi + dq)
          w = WSR_Result(WSR_Outside, false);
        else if (qa > q_lo + dq && qb < q_hi - dq)
          w = WSR_Result(WSR_Inside, false);
        else
          w = WSR_Result(WSR_Edge, false);

        // Near miss -- in the layer, but turning round inside it.
        if ((f_entry || f_exit) && w.m_wsr == WSR_Inside)
          w.m_wsr = WSR_Edge;

        // Endcap discs have a hole at small r. Sitting wholly inside it is a
        // miss with a reason; touching it is an edge, and either way m_in_gap
        // records that the absence of a hit is explained.
        if (!is_barrel && LI.has_r_range_hole() && w.m_wsr != WSR_Outside) {
          const float h_lo = LI.hole_r_min(), h_hi = LI.hole_r_max();
          if (qa > h_lo + dq && qb < h_hi - dq) {
            w = WSR_Result(WSR_Outside, true);
          } else if (qb > h_lo - dq && qa < h_hi + dq) {
            w.m_wsr = WSR_Edge;
            w.m_in_gap = true;
          }
        }
      }

      prim_tcand_ptrs[i]->m_wsr = w;

      switch (w.m_wsr) {
        case WSR_Inside:  V2P2_COUNT(n_wsr_inside);  break;
        case WSR_Edge:    V2P2_COUNT(n_wsr_edge);    break;
        case WSR_Outside: V2P2_COUNT(n_wsr_outside); break;
        default: break;
      }
      if (w.m_in_gap)
        V2P2_COUNT(n_wsr_in_gap);

      dprintf("%d: WSR %d in_gap %d  (q %.3f - %.3f vs layer %.3f - %.3f, dq3 %.4f, fail %d/%d)\n",
              i, w.m_wsr, (int) w.m_in_gap, B.m_q_min[i], B.m_q_max[i], q_lo, q_hi,
              B.m_dq_track[i], B.m_sp1.fail_flag[i], B.m_sp2.fail_flag[i]);
    }
  }

  //----------------------------------------------------------------------------
  // Phase 3 -- walk the bin ranges and pre-select hits.
  //----------------------------------------------------------------------------

  void MkFinderV2p2::select_hits(LayerBatch &b) {
    SteeringParams::iterator &spi = *mp_steeringparams_iter;
    const int N_proc = b.N_proc;
    MkBins &B = b.B;
    MkBinLimits &BL_p = b.BL_p;
#ifdef MKFIT_TRACE
    int *tr_layersearch_ids = b.tr_layersearch_ids;
#endif

    // Prototype for extract hits

    HitBatch hb;
    int &fill_pos = hb.fill_pos;
    auto &is_plex = hb.is_plex;
    auto &prim_idcs = hb.prim_idcs;
    auto &hit_idcs = hb.hit_idcs;
    auto &hit_orig_idcs = hb.hit_orig_idcs;


    // Both sub-layers, each into its own reduction queue. See the design notes,
    // "Reduction and hit ordering".
    const int n_sub = m_rz_limits.m_is_double ? 2 : 1;
    for (int sub = 0; sub < n_sub; ++sub) {
      const bool is_sec_layer = (sub == 1);
      const int lay = is_sec_layer ? spi->m_layer_sec : spi->m_layer;
      const auto &L = mp_job->m_event_of_hits[lay];
      const auto &iteration_hit_mask = mp_job->get_mask_for_layer(lay);
      const auto &BL = is_sec_layer ? b.BL_s : BL_p;
      fill_pos = 0;

      for (int i = 0; i < N_proc; ++i) {

        using bidx_t = LayerOfHits::bin_index_t;
        using bcnt_t = LayerOfHits::bin_content_t;

        // The track does not reach this layer.
        if (Policy::use_wsr && b.ptc[i]->m_wsr.m_wsr == WSR_Outside)
          continue;

        // Line pre-cut, per candidate: the straight line m_sp1 -> m_sp2 in
        // (qbar, q) and (qbar, phi), and the hit-independent tolerance terms.
        // See doc/MkFinderV2p2-DesignNotes.md, "Line pre-cut".
        const bool barrel = m_rz_limits.m_is_barrel;
        bool pc_on = PreCut::q || PreCut::phi;
        float pc_qb1 = 0, pc_q1 = 0, pc_gq = 0, pc_phi1 = 0, pc_gphi = 0, pc_tol_q = 0, pc_tol_phi = 0;
        if (pc_on) {
          const float x1 = B.m_sp1.x[i], y1 = B.m_sp1.y[i], z1 = B.m_sp1.z[i];
          const float x2 = B.m_sp2.x[i], y2 = B.m_sp2.y[i], z2 = B.m_sp2.z[i];
          const float r1 = std::hypot(x1, y1), r2 = std::hypot(x2, y2);
          pc_qb1 = barrel ? r1 : z1;
          pc_q1 = barrel ? z1 : r1;
          const float dqb = (barrel ? r2 : z2) - pc_qb1;
          if (std::abs(dqb) < 1e-4f) {
            pc_on = false;  // no span to interpolate along
          } else {
            pc_gq = std::clamp(((barrel ? z2 : r2) - pc_q1) / dqb, -20.0f, 20.0f);
            pc_phi1 = vdt::fast_atan2f(y1, x1);
            pc_gphi = std::clamp(squashPhiGeneral(vdt::fast_atan2f(y2, x2) - pc_phi1) / dqb, -20.0f, 20.0f);
            pc_tol_q = PreCut::dq_slack * Window::dq_trk_fac * B.m_dq_track[i] * (1.0f + pc_gq * pc_gq);
            pc_tol_phi = PreCut::dphi_slack * Window::dphi_trk_fac * B.m_dphi_track[i];
          }
        }

        // The region rho counts over and divides by (Score::rho_region), see the end of
        // this function. Region 2 takes the cut's padding of this sub-layer, as the fetch.
        float rho_dphi = B.m_dphi_track[i], rho_dq = B.m_dq_track[i];
        if (Score::rho_region == 2) {
          const float phi_hit_term = Window::phi_per_hit ? Window::dphi_hit_fac * L.max_hit_phi_half_extent()
                                                         : Window::dphi_flat_rad;
          rho_dphi = std::min(Window::dphi_trk_fac * B.m_dphi_track[i] + phi_hit_term, 3.14f);
          rho_dq = Window::dq_trk_fac * B.m_dq_track[i] + Window::dq_hit_fac * L.max_hit_q_half_length();
          b.rho_dphi[i] = std::max(b.rho_dphi[i], rho_dphi);
          b.rho_dq[i] = std::max(b.rho_dq[i], rho_dq);
        }
        const float rho_phi_c = 0.5f * (B.m_phi_min[i] + B.m_phi_max[i]);
        const float rho_phi_half = 0.5f * (B.m_phi_max[i] - B.m_phi_min[i]) + rho_dphi;
        const float rho_q_lo = B.m_q_min[i] - rho_dq;
        const float rho_q_hi = B.m_q_max[i] + rho_dq;

        for (bidx_t qi = BL.q1[i]; qi != BL.q2[i]; ++qi) {
          for (bidx_t pi = BL.p1[i]; pi != BL.p2[i]; pi = L.phiMaskApply(pi + 1)) {

            // Dead regions -- Limit to central Q-bin ???
            // if (qi == qb && L.isBinDead(pi, qi) == true) {
            //   dprint("dead module for track in layer=" << L.layer_id() << " qb=" << qi << " pi=" << pi
            //                                            << " q=" << B.q_c[itrack] << " phi=" << B.phi_c[itrack]);
            //   m_XWsrResult[itrack].m_in_gap = true;
            // }

            auto pbi = L.phiQBinContent(pi, qi);
            for (bcnt_t hi = pbi.begin(); hi < pbi.end(); ++hi) {

              const unsigned int hi_orig = L.getOriginalHitIndex(hi);

              dprintf(" %d: P_HIT %3u %4u %5u [%5u]  %6.3f %6.3f %6.3f\n",
                i, pi, qi, hi, hi_orig, L.hit_phi(hi), L.hit_q(hi), L.hit_qbar(hi));

              if (Score::rho_region == 0 ||
                  (std::abs(squashPhiGeneral(L.hit_phi(hi) - rho_phi_c)) < rho_phi_half &&
                   L.hit_q(hi) > rho_q_lo && L.hit_q(hi) < rho_q_hi))
                ++b.n_scanned[i];
#ifdef MKFIT_TRACE
              ++mp_event->tr_layersearch(tr_layersearch_ids[i]).n_hits_scanned;
#endif

              if (iteration_hit_mask && (*iteration_hit_mask)[hi_orig]) {
                dprintf("Yay, denying masked hit on layer %u, hi %u, orig idx %u\n",
                        L.layer_info().layer_id(), hi, hi_orig);
#ifdef MKFIT_TRACE
                ++mp_event->tr_layersearch(tr_layersearch_ids[i]).n_hits_masked;
#endif
                continue;
              }

              if (pc_on) {
                const float dqb_h = L.hit_qbar(hi) - pc_qb1;
                const float qbar_term = barrel ? PreCut::qbar_fac * L.hit_qbar_half_extent(hi) : 0.0f;
                bool reject = false;
                if (PreCut::q) {
                  const float dev = std::abs(L.hit_q(hi) - (pc_q1 + dqb_h * pc_gq));
                  reject = dev > pc_tol_q + Window::dq_hit_fac * L.hit_q_half_length(hi) + std::abs(pc_gq) * qbar_term;
                }
                if (!reject && PreCut::phi) {
                  const float dev = std::abs(squashPhiGeneral(L.hit_phi(hi) - (pc_phi1 + dqb_h * pc_gphi)));
                  const float hit_term = Window::phi_per_hit ? Window::dphi_hit_fac * L.hit_phi_half_extent(hi)
                                                            : Window::dphi_flat_rad;
                  reject = dev > pc_tol_phi + hit_term + std::abs(pc_gphi) * qbar_term;
                }
                if (reject) {
#ifdef MKFIT_TRACE
                  ++mp_event->tr_layersearch(tr_layersearch_ids[i]).n_hits_precut;
#endif
                  continue;
                }
              }

              // Try preloading Hits for the next step ... probably not really relevant.
              _mm_prefetch(&L.refHit(hi_orig), _MM_HINT_T0);

              prim_idcs[fill_pos] = i;
              hit_idcs[fill_pos] = hi;
              hit_orig_idcs[fill_pos] = hi_orig;
              is_plex.copyIn(fill_pos, B.m_isp, i);

              if (++fill_pos == NN) {
                preselect_hit_batch(b, hb, L, NN, is_sec_layer);
                fill_pos = 0;
              }
            }
          }
        }

        // Done with one PrimTCandRep. We might have enough hits to go into full KalmanProp.
        // Or wait for a change in ccand.
        // But there will be more work to be done, the overlaps, the other layer ...
        // ... so let's see.
      }
      if (fill_pos > 0) {
        preselect_hit_batch(b, hb, L, fill_pos, is_sec_layer);
      }
    }

    // Local hit density, hits per cm^2, over the window that was actually walked.
    // The window rather than the pre-selected hits, so the estimate does not
    // depend on the cut it is about to feed. The phi extent becomes a length at
    // the crossing radius. Score::rho_region 0 counts every hit of the fetched bins
    // (a larger region set by the fetch) over the window's area; 1 and 2 count over the
    // same region they divide by, the window or the cut region.
    for (int i = 0; i < N_proc; ++i) {
      const float r = std::max(0.1f, std::hypot(B.m_sp2.x[i], B.m_sp2.y[i]));
      const bool cut_region = Score::rho_region == 2;
      const float pad_phi = cut_region ? b.rho_dphi[i] : B.m_dphi_track[i];
      const float pad_q   = cut_region ? b.rho_dq[i]   : B.m_dq_track[i];
      const float w_phi = (B.m_phi_max[i] - B.m_phi_min[i]) + 2.0f * pad_phi;
      const float w_q   = (B.m_q_max[i] - B.m_q_min[i]) + 2.0f * pad_q;
      const float area  = std::max(1e-4f, w_phi * r * w_q);
      // Over the window or the cut region the count can be 0 with a hit about to be
      // scored, so the scored hit is counted too.
      const float n_rho = Score::rho_region != 0 ? b.n_scanned[i] + 1.0f : (float) b.n_scanned[i];
      b.log_rho[i] = std::log(std::max(1e-6f, n_rho / area));
      b.ptc[i]->m_log_rho = b.log_rho[i];
#if defined(MKFIT_STANDALONE)
      // Diag::final_beam_purity: hits in the density region, the expected background
      // of this layer step, per (seed, layer), forward search only.
      if (Diag::final_beam_purity && spi.type() == SteeringParams::IT_FwdSearch) {
        const int seed = b.ptc[i]->ccand().seed_origin_index();
        v2p2_final_beam_bg_record(seed, spi->m_layer, b.n_scanned[i]);
        if (spi->m_layer_sec >= 0)
          v2p2_final_beam_bg_record(seed, spi->m_layer_sec, b.n_scanned[i]);
      }
#endif
    }
  }

  //----------------------------------------------------------------------------
  // surface_referenced_dq() -- the track's q error referenced to the hit's
  // module plane rather than to the fixed path length it was transported to:
  //     v = e_q - ((e_q.p^)/(n^.p^)) n^ ,   sigma_q^2 = v^T C v
  // with n^ the module normal. Returns 3 sigma_q, the convention of
  // MkBins::m_dq_track, which is also the fallback. See
  // doc/MkFinderV2p2-DesignNotes.md, "Search window".
  //----------------------------------------------------------------------------

  float MkFinderV2p2::surface_referenced_dq(float dq_track_fallback,
                                            const MkBinTrackCovExtract &TCE, int pi,
                                            const mini_propagators::StatePlex &h3_state, int h,
                                            const MPlex3V &module_norm, bool is_barrel) {
    const float px = h3_state.px[h], py = h3_state.py[h], pz = h3_state.pz[h];
    const float nx = module_norm(h, 0, 0), ny = module_norm(h, 1, 0), nz = module_norm(h, 2, 0);
    const float np = nx * px + ny * py + nz * pz;      // (n^.p) -- |p| cancels below
    // q direction: z in the barrel -- so e_q.p is just pz, with no dot product
    // and no hipo -- and r^ in the endcap.
    float ex = 0.0f, ey = 0.0f, ez = 1.0f, eqp = pz;
    bool ok = (np != 0.0f);
    if (ok && !is_barrel) {
      const float rr = hipo(h3_state.x[h], h3_state.y[h]);
      ok = (rr > 0.0f);
      if (ok) {
        ex = h3_state.x[h] / rr; ey = h3_state.y[h] / rr; ez = 0.0f;
        eqp = ex * px + ey * py;
      }
    }
    if (!ok)
      return dq_track_fallback;

    // f = (e_q.p^)/(n^.p^); |p| cancels, so no normalisation is needed. Clamped:
    // it diverges only at grazing incidence on the module, which is not an
    // operating point (the cluster is then too wide in phi to be a hit).
    const float f = std::clamp(eqp / np, -20.0f, 20.0f);
    const float v0 = ex - f * nx, v1 = ey - f * ny, v2 = ez - f * nz;
    const float var = v0 * v0 * TCE.m_cov_0_0[pi] + v1 * v1 * TCE.m_cov_1_1[pi] +
                      v2 * v2 * TCE.m_cov_2_2[pi] +
                      2.0f * (v0 * v1 * TCE.m_cov_0_1[pi] + v0 * v2 * TCE.m_cov_0_2[pi] +
                              v1 * v2 * TCE.m_cov_1_2[pi]);
    if (var <= 0.0f)
      return dq_track_fallback;
    return 3.0f * std::sqrt(var);
  }

  //----------------------------------------------------------------------------
  // Phase 3b -- pre-select one NN-wide batch of (candidate, hit) pairs.
  //
  // Solves the Hermite cubic onto each hit's own module plane, applies the
  // dq / dphi pre-selection cut, and pushes survivors into the candidate's
  // bounded priority queue. Was a lambda inside the monolith.
  //----------------------------------------------------------------------------

  void MkFinderV2p2::preselect_hit_batch(LayerBatch &b, HitBatch &hb,
                                         const LayerOfHits &L, int N_proc_hits,
                                         bool is_sec_layer) {
    SteeringParams::iterator &spi = *mp_steeringparams_iter;
    MkBins &B = b.B;
    MkBinTrackCovExtract &TCE = b.TCE;
    PrimTCandRep **prim_tcand_ptrs = b.ptc;
    mini_propagators::Hermite3D &H = b.H;
    namespace mp = mini_propagators;
    auto &is_plex = hb.is_plex;
    auto &prim_idcs = hb.prim_idcs;
    auto &hit_idcs = hb.hit_idcs;
    auto &hit_orig_idcs = hb.hit_orig_idcs;
#ifdef MKFIT_TRACE_PROP_COMPARE
    auto &h_plex = hb.h_plex;
#endif
#ifdef MKFIT_TRACE
    int *tr_layersearch_ids = b.tr_layersearch_ids;
#endif
    (void) spi; (void) TCE;

      MPlex3V module_pos;
      MPlex3V module_norm;
      mp::Hermite3D h3d;

#ifdef MKFIT_TRACE
      int tr_hitmatch_ids[NN];
      MPlex3V module_xdir, module_ydir;
#endif

      // Extract hit / target module data
      for (int h = 0; h < N_proc_hits; ++h) {
        const Hit &hit = L.refHit(hit_orig_idcs[h]);
        unsigned int mid = hit.detIDinLayer();
        const ModuleInfo &mi = L.layer_info().module_info(mid);
        module_pos.copyIn(h, mi.pos.Array());
        module_norm.copyIn(h, mi.zdir.Array());

        h3d.copyIn(h, H, prim_idcs[h]);

#ifdef MKFIT_TRACE
        module_xdir.copyIn(h, mi.xdir.Array());
        module_ydir.copyIn(h, mi.calc_ydir().Array());
        PrimTCandRep &ptc = * prim_tcand_ptrs[ prim_idcs[h] ];
        int sim_lbl = mp_event->simInfoForCurrentSeed(ptc.ccand().seed_origin_index()).label;
        const MCHitInfo &mchinfo = mp_event->simHitsInfo_[L.refHit(hit_orig_idcs[h]).mcHitID()];
        int hit_lbl = mchinfo.mcTrackID();

        // L.layer_id(), not spi->m_layer: with the secondary sub-layer scanned
        // these differ, and `layer` has to name the layer the hit is actually in.
        TrHitMatch hm { -1, ptc.tcand().m_trace_state_id, tr_layersearch_ids[ prim_idcs[h] ],
                        (int) L.layer_id(), (int) hit_orig_idcs[h], sim_lbl == hit_lbl };
        hm.is_sec_layer = is_sec_layer;
        tr_hitmatch_ids[h] = mp_event->trace_hitmatch(std::move(hm)).id;
#endif
      }

#ifdef MKFIT_TRACE_PROP_COMPARE
      is_plex.propagate_to_plane(mp::PA_Line, module_pos, module_norm, h_plex, true);
#endif

      mp::Hermite3DOnPlane h3dop;
      h3dop.init_coeffs(h3d, module_pos, module_norm);
#ifdef DEBUG
      // Distance-to-plane at the two layer bounding surfaces (d0, d1) and at the
      // Newton start point (t2, d2) -- debug-only, so keep them inside the ifdef
      // rather than trusting the optimiser to drop the calls.
      MPlexQF d0, d1, d2, t2;
      h3dop.evaluate(0.0f, d0);
      h3dop.evaluate(1.0f, d1);
      t2 = h3dop.m_T;
      h3dop.evaluate(h3dop.m_T, d2);
#endif
      // NOTE: solve() is a SINGLE Newton step, not an iteration to convergence.
      // d3 below is what says whether that one step was enough, so it goes into
      // the trace and not just into DEBUG.
      h3dop.solve();
#if defined(MKFIT_TRACE) || defined(DEBUG)
      MPlexQF d3;
      h3dop.evaluate(h3dop.m_T, d3); // residual distance to plane after the solve
#endif

      mp::StatePlex h3_state;
      h3d.evaluate(h3dop.m_T, h3_state);
      // h3_state.dalpha calculated below, as needed

      // Post-process hits into the per-candidate queues. The cut is taken on
      // h3_state, the cubic solved onto the hit's module plane.
      for (int h = 0; h < N_proc_hits; ++h) {
        PrimTCandRep &ptc = * prim_tcand_ptrs[ prim_idcs[h] ];
        float q, ddq, phi, ddphi;
        if (m_rz_limits.m_is_barrel) {
          q = h3_state.z[h];
        } else {
          q = hipo(h3_state.x[h], h3_state.y[h]);
        }
        ddq = std::abs(q - L.hit_q(hit_idcs[h]));
        phi = vdt::fast_atan2f(h3_state.y[h], h3_state.x[h]);
        ddphi = cdist(std::abs(phi - L.hit_phi(hit_idcs[h])));

        // q error referenced to this module's plane, see surface_referenced_dq().
        const float dq_trk = Window::surface_q
          ? surface_referenced_dq(B.m_dq_track[prim_idcs[h]], TCE, prim_idcs[h],
                                  h3_state, h, module_norm, m_rz_limits.m_is_barrel)
          : B.m_dq_track[prim_idcs[h]];

        // The pre-selection cut. See doc/MkFinderV2p2-DesignNotes.md, "Search window".
        const float dq_cut = Window::dq_trk_fac * dq_trk +
                             Window::dq_hit_fac * L.hit_q_half_length(hit_idcs[h]);
        const float dphi_cut = Window::dphi_trk_fac * B.m_dphi_track[prim_idcs[h]] +
                               (Window::phi_per_hit ? Window::dphi_hit_fac * L.hit_phi_half_extent(hit_idcs[h])
                                                    : Window::dphi_flat_rad);
        const bool dqdphi_presel = ddq < dq_cut && ddphi < dphi_cut;

        // To be moved down, only for hits that pass pre-selection, needed here for printout.
        // Could be vectorized if we repack binnor stuff.
        h3_state.dalpha[h] = B.m_sp1.dalpha[prim_idcs[h]] + h3dop.m_T[h]*(B.m_sp2.dalpha[prim_idcs[h]] - B.m_sp1.dalpha[prim_idcs[h]]);

#ifdef DEBUG
        // clang-format off
        dprintf("     SelHit %6.3f %6.3f %6.4f %7.5f   %6.4f   %s [dq = %d, dphi = %d]\n",
                L.hit_q(hit_idcs[h]), L.hit_phi(hit_idcs[h]),
                ddq, ddphi, h3_state.dalpha[h], dqdphi_presel ? "PASS" : "REJECT",
                (int) (ddq < dq_cut), (int) (ddphi < dphi_cut));
        dprintf("       ddq=%.3f, dq_track=%.4f, hit_q_half_len=%.4f, dq_cut=%.4f\n",
                ddq, B.m_dq_track[prim_idcs[h]], L.hit_q_half_length(hit_idcs[h]), dq_cut);

        dprintf("      H3 d0=%.4f d1=%.4f -> d2=%e t2=%e -> d3=%e t3=%e ... dalpha=%6.4f\n",
               d0[h], d1[h], d2[h], t2[h], d3[h], h3dop.m_T[h],
               h3_state.dalpha[h]);
        // The two cheap propagations onto this module plane, side by side.
        // t3 outside [0,1] means the Hermite is extrapolating past the layer.
        dprintf("      H3   pos %8.4f %8.4f %8.4f | mom %8.4f %8.4f %8.4f\n",
                h3_state.x[h], h3_state.y[h], h3_state.z[h],
                h3_state.px[h], h3_state.py[h], h3_state.pz[h]);
#ifdef MKFIT_TRACE_PROP_COMPARE
        dprintf("      LINE pos %8.4f %8.4f %8.4f | mom %8.4f %8.4f %8.4f\n",
                h_plex.x[h], h_plex.y[h], h_plex.z[h],
                h_plex.px[h], h_plex.py[h], h_plex.pz[h]);
        dprintf("      H3 - LINE %+9.6f %+9.6f %+9.6f -- |d| = %.6f\n",
                h3_state.x[h] - h_plex.x[h], h3_state.y[h] - h_plex.y[h], h3_state.z[h] - h_plex.z[h],
                hipo(hipo(h3_state.x[h] - h_plex.x[h], h3_state.y[h] - h_plex.y[h]),
                     h3_state.z[h] - h_plex.z[h]));
#endif
        // clang-format on
#endif

#ifdef MKFIT_TRACE
        auto &tr_hitmatch = mp_event->tr_hitmatch(tr_hitmatch_ids[h]);
        // h3_state -- the state the cuts just above were taken on.
        tr_hitmatch.kine_on_plane = statep2bivec3(h3_state, h);
        tr_hitmatch.t_hermite = h3dop.m_T[h];
        tr_hitmatch.d_plane_h3 = d3[h];
        tr_hitmatch.dphi = ddphi;
        tr_hitmatch.dq = ddq;
        tr_hitmatch.hit_q_half_len = L.hit_q_half_length(hit_idcs[h]);
        tr_hitmatch.passed_preselect = dqdphi_presel;

        auto &tr_ls = mp_event->tr_layersearch(tr_layersearch_ids[ prim_idcs[h] ]);
        if (dqdphi_presel)
          ++tr_ls.n_hits_presel;

#ifdef MKFIT_TRACE_PROP_COMPARE
        tr_hitmatch.kine_on_plane_cmp = statep2bivec3(h_plex, h);
#endif

        // residuals
        EVec3 res = EVec3(h3_state.x[h], h3_state.y[h], h3_state.z[h]) - hit2pos(L.refHit(hit_orig_idcs[h]));
        tr_hitmatch.residual_x = res.Dot( EVec3(module_xdir(h, 0, 0), module_xdir(h, 1, 0), module_xdir(h, 2, 0)) );
        tr_hitmatch.residual_y = res.Dot( EVec3(module_ydir(h, 0, 0), module_ydir(h, 1, 0), module_ydir(h, 2, 0)) );
        tr_hitmatch.residual_z = res.Dot( EVec3(module_norm(h, 0, 0), module_norm(h, 1, 0), module_norm(h, 2, 0)) );
#endif

        if (/*prop_fail || */ !dqdphi_presel)
          continue;

        // float dalpha = h_plex.dalpha[h];
        // is_plex might come from somewhere else, through another index.

        // The reduction is PER SUB-LAYER, so each sensor keeps its own budget.
        const int sl = is_sec_layer ? 1 : 0;
        auto do_pqueue_push = [&]() {
#ifdef MKFIT_TRACE
          ptc.pq_push(sl, { ddphi, hit_orig_idcs[h], hit_idcs[h], L.layer_id(), tr_hitmatch_ids[h], { h3_state, h, is_plex, h } } );
#else
          ptc.pq_push(sl, { ddphi, hit_orig_idcs[h], hit_idcs[h], L.layer_id(), { h3_state, h, is_plex, h } } );
#endif
        };

        if (ptc.m_pqueue_size[sl] < InLayer::max_presel_hits) {
          do_pqueue_push();
          ++ptc.m_pqueue_size[sl];
        } else if (ddphi < ptc.pq_top(sl).score) {
          ptc.pq_pop(sl);
          do_pqueue_push();
        }
      }
  }

  //----------------------------------------------------------------------------
  // Phase 4 -- move the pqueue survivors into the per-candidate hit list.
  //----------------------------------------------------------------------------

  void MkFinderV2p2::prepare_kalman_workload(LayerBatch &b) {
    const int N_proc = b.N_proc;
    PrimTCandRep **prim_tcand_ptrs = b.ptc;
#ifdef MKFIT_TRACE
    int *tr_layersearch_ids = b.tr_layersearch_ids;
#endif

    // Drain both sub-layer queues into one list, then order it by path length.
    // sub_rank (within a sub-layer) and full_rank (across both) are trace-only;
    // the pqueue pops worst first, hence the count-down.
    for (int i = 0; i < N_proc; ++i) {
      PrimTCandRep &ptc = * prim_tcand_ptrs[i];
#ifdef MKFIT_TRACE
      mp_event->tr_layersearch(tr_layersearch_ids[i]).n_hits_pqueue =
        ptc.m_pqueue_size[0] + ptc.m_pqueue_size[1];
#endif
      for (int sl = 0; sl < 2; ++sl) {
#ifdef MKFIT_TRACE
        int rank = ptc.m_pqueue_size[sl];
#endif
        while (ptc.m_pqueue_size[sl]) {
          --ptc.m_pqueue_size[sl];
          const auto &pqe = ptc.pq_top(sl);
          ptc.m_layer_hits.push_back( pqe );
#ifdef MKFIT_TRACE
          TrHitMatch &tr_hitmatch = mp_event->tr_hitmatch(pqe.tr_hitmatch_id);
          tr_hitmatch.sub_rank = rank--;
          tr_hitmatch.passed_pqueue = true;
#endif
          ptc.pq_pop(sl);
        }
      }

      // Path order for the in-layer search: dalpha is monotone in path length
      // for one candidate, but it is a signed turning angle, so its sign goes
      // with the charge. m_sp1 is the near bounding surface and m_sp2 the far
      // one in the search direction, so dalpha times the sign of its change
      // from m_sp1 to m_sp2 grows along the path, in either search direction.
      // Once per candidate is enough. See doc/MkFinderV2p2-DesignNotes.md,
      // "Reduction and hit ordering".
      const float dir = b.B.m_sp2.dalpha[i] >= b.B.m_sp1.dalpha[i] ? 1.0f : -1.0f;
      std::sort(ptc.m_layer_hits.begin(), ptc.m_layer_hits.end(),
                [dir](const PrimTCandRep::PQE &a, const PrimTCandRep::PQE &b) {
                  return dir * a.mixed_state.dalpha < dir * b.mixed_state.dalpha;
                });
    }

#ifdef MKFIT_TRACE
    // full_rank -- rank by score across both sub-layers. Counted rather than
    // sorted, so m_layer_hits keeps its path order.
    for (int i = 0; i < N_proc; ++i) {
      PrimTCandRep &ptc = * prim_tcand_ptrs[i];
      for (const auto &e : ptc.m_layer_hits) {
        int better = 0;
        for (const auto &o : ptc.m_layer_hits)
          better += (o.score < e.score);
        mp_event->tr_hitmatch(e.tr_hitmatch_id).full_rank = better + 1;
      }
    }
#endif
  }

  //----------------------------------------------------------------------------
  // Phase 5 -- propagate to each hit's module plane and run the Kalman update.
  //----------------------------------------------------------------------------

  void MkFinderV2p2::kalman_update(LayerBatch &b) {
    const int N_proc = b.N_proc;
    PrimTCandRep **prim_tcand_ptrs = b.ptc;

    auto do_kalman = [&](KalmanOpArgs& K) {
      K.compute_pars(); // Needed for conversion from mini_prop/bi-vec to std representation.
      K.do_kalman_stuff();
      K.reset();
    }; // end lambda do_kalman

    // KalmanProp, prim layer -- should also be done in-line as slots fill up.
    KalmanOpArgs koa;
    koa.prop_config = & mp_job->m_trk_info.prop_config();
#ifdef MKFIT_TRACE
    koa.mp_event = mp_event;
#endif

    for (int i = 0; i < N_proc; ++i) {
      PrimTCandRep &ptc = * prim_tcand_ptrs[i];
      int nlh = ptc.m_layer_hits.size();
      for (int lh = 0; lh < nlh; ++ lh) {
        const PrimTCandRep::PQE &pqe = ptc.m_layer_hits[lh];
        dprintf("scheduling %d %d %d %f\n", i, pqe.hit_index, pqe.hit_orig_index, pqe.mixed_state.dalpha);

        // Need to build all the Matriplexes for KalmanOperationPlane ... or some variant
        // There really needs to be an intermediate structure with packers so we can
        // set them up incrementally.

        TrackCand &tc = ptc.tcand();
        koa.item_begin(&ptc, { (int) pqe.hit_orig_index, pqe.layer });
        koa.load_state_err_chg(pqe.mixed_state, tc.state());
#ifdef MKFIT_TRACE
        koa.set_tr_hitmatch_id(pqe.tr_hitmatch_id);
#endif
        const auto &L = mp_job->m_event_of_hits[ pqe.layer ];
        const Hit &hit = L.refHit( pqe.hit_orig_index );
        unsigned int mid = hit.detIDinLayer();
        const ModuleInfo &mi = L.layer_info().module_info(mid);
        koa.load_hit_module(hit, mi);

        if (koa.item_finished()) {
          do_kalman(koa);
        }
      }
      // QQQQ hits stay in, somehow. Or it was leftover PrimTCands etc due to incomplete driver.
      ptc.m_layer_hits.clear();
    }
    if (koa.N_filled > 0) {
      do_kalman(koa);
    }

  }

  //----------------------------------------------------------------------------
  // Phase 6 -- best-hit acceptance, used with InLayer::comb off: take the
  // lowest-chi2 hit if it passes Policy::hit_chi2_cut, otherwise record the
  // fake HoT from fake_hit_index().
  //----------------------------------------------------------------------------

  void MkFinderV2p2::process_kalman_results(LayerBatch &b) {
    const int N_proc = b.N_proc;
    PrimTCandRep **prim_tcand_ptrs = b.ptc;

    // This, esp. the combinatorial part, should be done once prim-tcand is finished.
    // And, merging results, when ccand is finished.

    for (int i = 0; i < N_proc; ++i) {
      PrimTCandRep &ptc = * prim_tcand_ptrs[i];
      TrackCand &tc = ptc.tcand();

      // The track does not reach this layer: no HoT of any kind.
      if (Policy::use_wsr && ptc.m_wsr.m_wsr == WSR_Outside) {
        dprintf("Outside to tcand %d, layer skipped\n", i);
        V2P2_COUNT(n_layer_skipped);
        continue;
      }

#ifdef MKFIT_TRACE
      if (ptc.bChi2Acc < Policy::hit_chi2_cut || (Diag::force_mc && ptc.bIsMc)) {
#else
      if (ptc.bChi2Acc < Policy::hit_chi2_cut) {
#endif
        // XXXX Extra missed layer -- to check stuff / maxgrowth / scores etc
        // This is somewhat impure :)
        // Add a copy of the held-back candidate before adding the hit.
        if (ptc.bChi2 > 5.0f && ! ptc.mp_ccrep->m_ccand.is_full()) {
          dprintf("ExtraMissed to tcand %d\n", i);
          ptc.mp_ccrep->m_ccand.push_back(tc).addHitIdx(
              fake_hit_index(tc, ptc.m_wsr), m_rz_limits.layer_info_1().layer_id(), 0.0f);

#ifdef MKFIT_TRACE
          // QQQQQ the parent extraction will be different
          int pid = ptc.tcand().m_trace_state_id;
          int id = mp_event->trace_new_cand_state(pid, (*mp_steeringparams_iter)->m_layer, track2bivec3(ptc.tcand()), ptc.tcand().state());
          ptc.mp_ccrep->m_ccand.back().m_trace_state_id = id;
#endif
        }

        dprintf("Output to tcand %d, idx=%d, layer=%d, chi2=%f\n", i, ptc.bHot.index, ptc.bHot.layer, ptc.bChi2);
        tc.addHitIdx(ptc.bHot.index, ptc.bHot.layer, ptc.bChi2);
        tc.setState(ptc.bState);

#ifdef MKFIT_TRACE
        // QQQQQ the parent extraction will be different; also fix: step, proper state (what is it)

        // This is also best-hit hack
        int pid = ptc.tcand().m_trace_state_id;
        int id = mp_event->trace_new_cand_state(pid, (*mp_steeringparams_iter)->m_layer, track2bivec3(ptc.tcand()), ptc.tcand().state());

        auto &ku = mp_event->tr_kalmanupdate( mp_event->tr_hitmatch(ptc.b_tr_hitmatch_id).kalman_id );
        ku.accepted = true;
        ku.state_id_out = id;
        // At this point could also set chi2_trk ... but this might get tricky with multiple hits added per (multi-)layer.
        // local / (multi-)layer score might be more relevant.

        ptc.tcand().m_trace_state_id = id;
#endif
      } else {
        // XXXX Here, we need to handle double layers correctly.
        // For now we have singles.
        //
        // No hit taken: the fake HoT comes from fake_hit_index().
        const int fake_hit_idx = fake_hit_index(tc, ptc.m_wsr);
        dprintf("Missed to tcand %d, fake_hit_idx=%d\n", i, fake_hit_idx);
        switch (fake_hit_idx) {
          case Hit::kHitMissIdx:  V2P2_COUNT(n_hole);      break;
          case Hit::kHitEdgeIdx:  V2P2_COUNT(n_hot_edge);  break;
          case Hit::kHitInGapIdx: V2P2_COUNT(n_hot_gap);   break;
          case Hit::kHitStopIdx:  V2P2_COUNT(n_stop_holes); break;
          default: break;
        }
        if (fake_hit_idx == Hit::kHitStopIdx)
          tc.setStopReason(TrackCand::SR_TooManyHoles);
        tc.addHitIdx(fake_hit_idx, m_rz_limits.layer_info_1().layer_id(), 0.0f);

#ifdef MKFIT_TRACE
        // QQQQQ the parent extraction will be different; also fix: step, proper state (what is it)
        int pid = ptc.tcand().m_trace_state_id;
        int id = mp_event->trace_new_cand_state(pid, (*mp_steeringparams_iter)->m_layer, track2bivec3(ptc.tcand()), ptc.tcand().state());
        ptc.tcand().m_trace_state_id = id;
#endif
      }
    }

  }

  //----------------------------------------------------------------------------
  // In-layer combinatorial search: grow the tree of paths through the layer,
  // breadth-first by depth, over the path-ordered m_layer_hits. Nodes live in
  // m_sec_arena until end of layer. See doc/MkFinderV2p2-DesignNotes.md,
  // "In-layer combinatorial search".
  //----------------------------------------------------------------------------

  // InLayer::max_sec_depth caps the hits per path; InLayer::max_sec_depth_limit
  // sizes the chain array.

  int MkFinderV2p2::sec_node_slot() {
    if (m_sec_free.empty()) {
      m_sec_arena.emplace_back();
      return (int) m_sec_arena.size() - 1;
    }
    const int slot = m_sec_free.back();
    m_sec_free.pop_back();
    return slot;
  }

  void MkFinderV2p2::harvest_sec_nodes() {
    m_sec_harvested.clear();
    for (const auto &o : m_sec_out) {
      if ( ! (o.chi2_acc < Policy::hit_chi2_cut))   // also rejects NaN
        continue;
      SecTCandRep n;
      n.m_ptc        = o.ptc;
      n.m_parent_idx = o.parent_idx;
      n.m_hit_pos    = o.hit_pos;
      n.m_hot        = o.hot;
      n.m_state      = o.state;
      n.m_chi2       = o.chi2;

      // Features accumulate parent -> child. The parent is always already in the
      // arena: harvest runs once per depth, after the whole depth has flushed.
      const LayerStepFeatures *pf =
        (o.parent_idx >= 0) ? &m_sec_arena[o.parent_idx].m_feat : nullptr;
      LayerStepFeatures &f = n.m_feat;
      if (pf)
        f = *pf;
      else
        fill_step_geometry(f, *o.ptc, o.ptc->m_log_rho);
      if (pf)
        V2P2_COUNT(n_sec_deep);
      f.n_hits    = (pf ? pf->n_hits : 0) + 1;
      f.n_overlap = f.n_hits - 1;
      f.hole_kind = V2P2_NoHole;
      f.chi2_sum  = (pf ? pf->chi2_sum : 0.0f) + o.chi2;
      f.chi2_max  = std::max(pf ? pf->chi2_max : 0.0f, o.chi2);
      {
        const auto &L = mp_job->m_event_of_hits[o.hot.layer];
        const float qhl = L.hit_q_half_length(o.hit_in_layer);
        f.q_half_len_best = (pf && pf->n_hits > 0) ? std::min(pf->q_half_len_best, qhl) : qhl;
      }
      f.log_det_v_sum = (pf ? pf->log_det_v_sum : 0.0f) + std::log(std::max(1e-30f, o.det_v));
#ifdef MKFIT_TRACE
      n.m_tr_hitmatch_id = o.tr_hitmatch_id;
#endif
      o.ptc->m_has_sec_nodes = true;
      const int slot = sec_node_slot();   // may grow the arena: pf is dead by now
      m_sec_arena[slot] = n;              // a reused slot is overwritten whole
      o.ptc->mp_ccrep->m_sec_nodes.push_back(slot);
      m_sec_harvested.push_back(slot);
    }
    m_sec_out.clear();
    V2P2_COUNT_ADD(n_sec_nodes, (long) m_sec_harvested.size());
  }

  void MkFinderV2p2::select_completed_ccreps(const LayerBatch &b) {
    // The batch's PrimTCandReps stay alive: end_layer() clears m_primTCs.
    auto try_select = [this](CCandRep &ccrep) {
      if (ccrep.m_n_queued > 0 || ccrep.m_n_inflight > 0 || ccrep.m_selected)
        return;
      select_and_materialise(ccrep);
      ccrep.m_selected = true;
      V2P2_COUNT(n_early_selections);
      for (const int ni : ccrep.m_sec_nodes)
        m_sec_free.push_back(ni);
      ccrep.m_sec_nodes.clear();
    };
    for (int i = 0; i < b.N_proc; ++i)
      try_select(* b.ptc[i]->mp_ccrep);
    for (CCandRep *c : m_ccrep_ready)
      try_select(*c);
    m_ccrep_ready.clear();
  }

  //----------------------------------------------------------------------------
  // Running Kalman batches (InLayer::running_kalman). Same work as
  // expand_in_layer(), same nodes and the same per-node arithmetic, since every
  // Matriplex lane is independent; only which items share a batch changes.
  // Each CombCandidate counts its Kalman items in flight and its nodes still to
  // be expanded, so it can be selected as soon as both reach zero.
  //----------------------------------------------------------------------------

  void MkFinderV2p2::load_hit_module(KalmanOpArgs &koa, const PrimTCandRep::PQE &pqe) const {
    const auto &L = mp_job->m_event_of_hits[pqe.layer];
    const Hit &hit = L.refHit(pqe.hit_orig_index);
    koa.load_hit_module(hit, L.layer_info().module_info(hit.detIDinLayer()));
  }

  void MkFinderV2p2::rk_setup(KalmanOpArgs &koa, bool solve_plane) {
    koa.prop_config = &mp_job->m_trk_info.prop_config();
    koa.mp_out = &m_sec_out;
    koa.m_solve_plane = solve_plane;
#ifdef MKFIT_TRACE
    koa.mp_event = mp_event;
#endif
  }

  void MkFinderV2p2::rk_add_item(KalmanOpArgs &koa, PrimTCandRep &ptc, const PrimTCandRep::PQE &pqe,
                                 int parent_idx, int hit_pos) {
    koa.item_begin(&ptc, {(int) pqe.hit_orig_index, pqe.layer}, parent_idx, hit_pos, pqe.hit_index);
    if (parent_idx < 0)
      koa.load_state_err_chg(pqe.mixed_state, ptc.tcand().state());
    else
      koa.load_state_err_chg(m_sec_arena[parent_idx].m_state);
#ifdef MKFIT_TRACE
    koa.set_tr_hitmatch_id(pqe.tr_hitmatch_id);
#endif
    load_hit_module(koa, pqe);
    ++ptc.mp_ccrep->m_n_inflight;
    if (koa.item_finished())
      rk_flush(koa);
  }

  void MkFinderV2p2::rk_flush(KalmanOpArgs &koa) {
    if (koa.N_filled == 0)
      return;
    V2P2_COUNT(n_kalman_calls);
    V2P2_COUNT_ADD(n_kalman_lanes, koa.N_filled);
    if (!koa.m_solve_plane) {
      V2P2_COUNT(n_kalman_calls_d0);
      V2P2_COUNT_ADD(n_kalman_lanes_d0, koa.N_filled);
      koa.compute_pars();
    }
    koa.do_kalman_stuff();
    koa.reset();

    m_ccrep_touched.clear();
    for (const auto &o : m_sec_out)
      m_ccrep_touched.push_back(o.ptc->mp_ccrep);
    harvest_sec_nodes();
    for (const int ni : m_sec_harvested) {
      const SecTCandRep &n = m_sec_arena[ni];
      if (n.m_feat.n_hits < InLayer::max_sec_depth) {
        m_deep_queue.push_back(ni);
        ++n.m_ptc->mp_ccrep->m_n_inflight;
      }
    }
    // Items retire after their children are queued, so a count never passes
    // through zero while work remains.
    for (CCandRep *c : m_ccrep_touched)
      if (--c->m_n_inflight == 0 && c->m_n_queued == 0)
        m_ccrep_ready.push_back(c);
  }

  void MkFinderV2p2::rk_expand_node(int ni) {
    // By index throughout: a flush inside this loop may grow the arena. The node
    // itself cannot be freed meanwhile, its expansion still counts as in flight.
    PrimTCandRep &ptc = *m_sec_arena[ni].m_ptc;
    const int nlh = (int) ptc.m_layer_hits.size();
    for (int lh = m_sec_arena[ni].m_hit_pos + 1; lh < nlh; ++lh) {
      const PrimTCandRep::PQE &pqe = ptc.m_layer_hits[lh];
      bool same_module = false;
      for (int ci = ni; ci >= 0 && !same_module; ci = m_sec_arena[ci].m_parent_idx) {
        const SecTCandRep &an = m_sec_arena[ci];
        if (an.m_hot.layer != pqe.layer)
          continue;
        const auto &La = mp_job->m_event_of_hits[an.m_hot.layer];
        same_module = La.refHit(an.m_hot.index).detIDinLayer() == La.refHit(pqe.hit_orig_index).detIDinLayer();
      }
      if (same_module) {
        V2P2_COUNT(n_same_module_vetoed);
        continue;
      }
      rk_add_item(m_koa_deep, ptc, pqe, ni, lh);
    }
  }

  void MkFinderV2p2::rk_pump() {
    while (!m_deep_queue.empty()) {
      const int ni = m_deep_queue.back();
      m_deep_queue.pop_back();
      CCandRep *c = m_sec_arena[ni].m_ptc->mp_ccrep;
      rk_expand_node(ni);
      if (--c->m_n_inflight == 0 && c->m_n_queued == 0)
        m_ccrep_ready.push_back(c);
    }
  }

  void MkFinderV2p2::expand_in_layer_running(LayerBatch &b) {
    rk_setup(m_koa_d0, false);
    rk_setup(m_koa_deep, true);
    m_sec_out.clear();
    for (int i = 0; i < b.N_proc; ++i) {
      PrimTCandRep &ptc = *b.ptc[i];
      if (Policy::use_wsr && ptc.m_wsr.m_wsr == WSR_Outside)
        continue;
      const int nlh = (int) ptc.m_layer_hits.size();
      for (int lh = 0; lh < nlh; ++lh)
        rk_add_item(m_koa_d0, ptc, ptc.m_layer_hits[lh], -1, lh);
      rk_pump();
    }
  }

  void MkFinderV2p2::rk_drain() {
    rk_setup(m_koa_d0, false);
    rk_setup(m_koa_deep, true);
    while (m_koa_d0.N_filled > 0 || m_koa_deep.N_filled > 0 || !m_deep_queue.empty()) {
      rk_flush(m_koa_d0);
      rk_pump();
      rk_flush(m_koa_deep);
      rk_pump();
    }
    m_ccrep_ready.clear();   // end_layer() selects everything not yet selected
  }

  // The parts of a layer step that do not depend on which hits were taken: where
  // the step goes, which way, and what the candidate was when it arrived. Filled
  // once per path root and then carried down the tree.
  void MkFinderV2p2::fill_step_geometry(LayerStepFeatures &f, const PrimTCandRep &ptc,
                                        float log_rho) const {
    const TrackCand &tc = const_cast<PrimTCandRep &>(ptc).tcand();
    const LayerInfo &li = m_rz_limits.layer_info_1();
    // layer_from is the layer of the candidate's LAST FOUND HIT, so the pair
    // names the propagation that actually happens. A candidate that missed the
    // previous layer steps across two, and keying anything on plan adjacency
    // instead would attribute it to a step that never took place.
    f.layer_from = (short) tc.getLastFoundHitLyr();
    f.layer_to   = (short) li.layer_id();
    f.is_outward = m_rz_limits.is_outward();
    f.is_pixel   = li.is_pixel();
    f.is_barrel  = li.is_barrel();
    f.pt         = tc.pT();
    f.step       = (short) tc.nTotalHits();
    f.n_found_so_far = (short) tc.nFoundHits();
    f.n_holes_so_far = (short) tc.nAllMinusOneHits();
    f.log_rho = log_rho;
  }

  void MkFinderV2p2::expand_in_layer(LayerBatch &b) {
    const int N_proc = b.N_proc;
    PrimTCandRep **prim_tcand_ptrs = b.ptc;

    // NB: the arena is NOT rewound here. It has to live until end of layer, when
    // the whole CombCandidate's alternatives are selected between; end_layer()
    // rewinds it once for all of them.
    m_sec_out.clear();

    KalmanOpArgs koa;
    koa.prop_config = & mp_job->m_trk_info.prop_config();
    koa.mp_out = & m_sec_out;
#ifdef MKFIT_TRACE
    koa.mp_event = mp_event;
#endif

    auto load_module = [&](const PrimTCandRep::PQE &pqe) {
      const auto &L = mp_job->m_event_of_hits[ pqe.layer ];
      const Hit &hit = L.refHit( pqe.hit_orig_index );
      koa.load_hit_module(hit, L.layer_info().module_info(hit.detIDinLayer()));
    };
    auto flush = [&]() {
      if (koa.N_filled == 0)
        return;
      // Lane occupancy, split by depth: depth 0 is the batch shape the best-hit
      // path always had, deeper ones are new and are the ones that could run
      // ragged. Breadth-first BY DEPTH exists so they do not.
      V2P2_COUNT(n_kalman_calls);
      V2P2_COUNT_ADD(n_kalman_lanes, koa.N_filled);
      if ( ! koa.m_solve_plane) {
        V2P2_COUNT(n_kalman_calls_d0);
        V2P2_COUNT_ADD(n_kalman_lanes_d0, koa.N_filled);
        koa.compute_pars();   // propPar is an INPUT on the sPerp path, an output on the solve path
      }
      koa.do_kalman_stuff();
      koa.reset();
    };

    // Depth 0. The Hermite has already solved the crossing for every
    // pre-selected hit, so this is the cheap sPerp path and is exactly what the
    // best-hit code does -- it just keeps every outcome instead of the winner.
    koa.m_solve_plane = false;
    for (int i = 0; i < N_proc; ++i) {
      PrimTCandRep &ptc = * prim_tcand_ptrs[i];
      if (Policy::use_wsr && ptc.m_wsr.m_wsr == WSR_Outside)
        continue;
      const int nlh = (int) ptc.m_layer_hits.size();
      for (int lh = 0; lh < nlh; ++lh) {
        const PrimTCandRep::PQE &pqe = ptc.m_layer_hits[lh];
        koa.item_begin(&ptc, { (int) pqe.hit_orig_index, pqe.layer }, -1, lh, pqe.hit_index);
        koa.load_state_err_chg(pqe.mixed_state, ptc.tcand().state());
#ifdef MKFIT_TRACE
        koa.set_tr_hitmatch_id(pqe.tr_hitmatch_id);
#endif
        load_module(pqe);
        if (koa.item_finished())
          flush();
      }
    }
    flush();
    harvest_sec_nodes();
    std::vector<int> &frontier = m_sec_frontier;
    frontier.swap(m_sec_harvested);

    // Depths 1 and up. The starting state is now a node's UPDATED state, for
    // which no crossing has been solved, so propagate-to-plane solves it.
    for (int depth = 1; depth < InLayer::max_sec_depth && ! frontier.empty(); ++depth) {
      koa.m_solve_plane = true;
      for (const int ni : frontier) {
        // By index, not by reference: the arena grows under us only at harvest,
        // but the discipline is what keeps the indices the handles.
        PrimTCandRep &ptc = * m_sec_arena[ni].m_ptc;
        const int nlh = (int) ptc.m_layer_hits.size();
        for (int lh = m_sec_arena[ni].m_hit_pos + 1; lh < nlh; ++lh) {
          const PrimTCandRep::PQE &pqe = ptc.m_layer_hits[lh];

          // A second hit from a module already on this path is a split cluster or
          // another track's hit, not an overlap: skip it. The walk back is at most
          // InLayer::max_sec_depth long.
          bool same_module = false;
          for (int ci = ni; ci >= 0 && !same_module; ci = m_sec_arena[ci].m_parent_idx) {
            const SecTCandRep &an = m_sec_arena[ci];
            if (an.m_hot.layer != pqe.layer)
              continue;
            const auto &La = mp_job->m_event_of_hits[an.m_hot.layer];
            same_module = La.refHit(an.m_hot.index).detIDinLayer() ==
                          La.refHit(pqe.hit_orig_index).detIDinLayer();
          }
          if (same_module) {
            V2P2_COUNT(n_same_module_vetoed);
            continue;
          }

          koa.item_begin(&ptc, { (int) pqe.hit_orig_index, pqe.layer }, ni, lh, pqe.hit_index);
          koa.load_state_err_chg(m_sec_arena[ni].m_state);
#ifdef MKFIT_TRACE
          koa.set_tr_hitmatch_id(pqe.tr_hitmatch_id);
#endif
          load_module(pqe);
          if (koa.item_finished())
            flush();
        }
      }
      flush();
      harvest_sec_nodes();
      frontier.swap(m_sec_harvested);
    }

  }

  //----------------------------------------------------------------------------
  // End-of-layer selection: every in-layer path, every candidate as the hole it
  // would record, and every candidate that did not enter the layer compete in
  // one sort on the additive layer-step score; the top maxCandsPerSeed survive.
  // See doc/MkFinderV2p2-DesignNotes.md, "End-of-layer selection".
  //----------------------------------------------------------------------------

  // Keep the best-scoring stopped candidate of this CombCandidate. The score is
  // the in-flight accumulator, which is on one scale within a seed -- the only
  // comparison made here.
  void MkFinderV2p2::offer_best_short(CombCandidate &ccand, const TrackCand &tc) const {
    V2P2_COUNT(n_best_short_offered);
    if (ccand.refBestShortCand().combCandidate() == nullptr ||
        tc.score() > ccand.refBestShortCand().score()) {
      ccand.setBestShortCand(tc);
      V2P2_COUNT(n_best_short_taken);
    }
  }

  void MkFinderV2p2::select_and_materialise(CCandRep &ccrep) {
    CombCandidate &ccand = ccrep.m_ccand;
    const int max_cands = ccrep.m_max_cands > 0 ? ccrep.m_max_cands : ccand.capacity();

    m_sel.clear();

    // Which TrackCands got a PrimTCandRep this layer. The rest are entered
    // unchanged below.
    bool has_ptc[64] = {false};
    const int n_tc = (int) ccand.size();

    for (auto &ptc : ccrep.m_primTCs) {
      if (ptc.m_origin_tcand_index < 64)
        has_ptc[ptc.m_origin_tcand_index] = true;

      TrackCand &tc = ptc.tcand();

      // Every PrimTCandRep also competes as the hole it would record, including
      // those that found paths.
      // WSR_Outside is not a hole: the candidate competes unchanged.
      if (Policy::use_wsr && ptc.m_wsr.m_wsr == WSR_Outside) {
        V2P2_COUNT(n_layer_skipped);
        m_sel.push_back({ptc.m_origin_tcand_index, -1, 0, false, tc.score()});
      } else {
        const int fake = fake_hit_index(tc, ptc.m_wsr);
        LayerStepFeatures f;
        fill_step_geometry(f, ptc, 0.0f);
        f.hole_kind = v2p2_hole_kind(fake);
        m_sel.push_back({ptc.m_origin_tcand_index, -1, fake, true,
                         tc.score() + v2p2_layer_step_score(f)});
      }
    }

    for (const int ni : ccrep.m_sec_nodes) {
      const SecTCandRep &n = m_sec_arena[ni];
      m_sel.push_back({n.m_ptc->m_origin_tcand_index, ni, 0, false,
                       n.m_ptc->tcand().score() + v2p2_layer_step_score(n.m_feat)});
    }

    // Stopped or pull-in-skipped TrackCands, unchanged, at their own score. With
    // best-short on (outward only), stopped ones leave the beam instead and the
    // best is kept on the CombCandidate.
    const bool best_short = InLayer::best_short && m_rz_limits.is_outward();
    for (int ic = 0; ic < n_tc; ++ic) {
      if (ic < 64 && has_ptc[ic])
        continue;
      if (best_short && ccand[ic].getLastHitIdx() == Hit::kHitStopIdx) {
        offer_best_short(ccand, ccand[ic]);
        continue;
      }
      m_sel.push_back({ic, -1, 0, false, ccand[ic].score()});
    }

    if (m_sel.empty())
      return;

    if (InLayer::sister_hole_drop)
      drop_sister_dominated_holes(ccand);

    // Partial sort: only the top max_cands matter.
    const int n_keep = std::min((int) m_sel.size(), max_cands);
    std::partial_sort(m_sel.begin(), m_sel.begin() + n_keep, m_sel.end(),
                      [](const SelEntry &a, const SelEntry &b) { return a.score > b.score; });

    // Reserve a hole slot: keep the best candidate that declined the layer, even
    // if outranked, so one branch stays open to an earlier hit being wrong.
    if (InLayer::reserve_hole_slot && n_keep > 1) {
      bool kept_a_decliner = false;
      for (int k = 0; k < n_keep && !kept_a_decliner; ++k)
        kept_a_decliner = (m_sel[k].node_idx < 0);
      if (!kept_a_decliner) {
        int best_decliner = -1;
        for (int k = n_keep; k < (int) m_sel.size(); ++k)
          if (m_sel[k].node_idx < 0 &&
              (best_decliner < 0 || m_sel[k].score > m_sel[best_decliner].score))
            best_decliner = k;
        if (best_decliner >= 0) {
          std::swap(m_sel[n_keep - 1], m_sel[best_decliner]);
          V2P2_COUNT(n_hole_slot_reserved);
        }
      }
    }

#if defined(MKFIT_STANDALONE)
    count_dominated_kept(ccand, n_keep);
#endif

    // Build the survivors as copies BEFORE touching the CombCandidate: several
    // winners can descend from the same TrackCand -- that is what branching is --
    // so the sources have to stay readable while the copies are made.
    m_new_cands.clear();
    for (int k = 0; k < n_keep; ++k) {
      const SelEntry &e = m_sel[k];
      m_new_cands.push_back(ccand[e.tcand_idx]);
      TrackCand &nc = m_new_cands.back();

      if (e.node_idx >= 0) {
        int chain[InLayer::max_sec_depth_limit];
        int n_chain = 0;
        for (int ci = e.node_idx; ci >= 0; ci = m_sec_arena[ci].m_parent_idx)
          chain[n_chain++] = ci;
        for (int c = n_chain - 1; c >= 0; --c) {
          const SecTCandRep &n = m_sec_arena[chain[c]];
          nc.addHitIdx(n.m_hot.index, n.m_hot.layer, n.m_chi2);
          if (c != n_chain - 1) {
            nc.incOverlapCount();   // every hit past the first in one layer IS the overlap case
            // Counted by module: same module against a genuine overlap.
#if defined(MKFIT_STANDALONE)
            const SecTCandRep &pn = m_sec_arena[chain[c + 1]];
            const auto &Ln = mp_job->m_event_of_hits[n.m_hot.layer];
            const auto &Lp = mp_job->m_event_of_hits[pn.m_hot.layer];
            if (n.m_hot.layer == pn.m_hot.layer &&
                Ln.refHit(n.m_hot.index).detIDinLayer() == Lp.refHit(pn.m_hot.index).detIDinLayer())
              V2P2_COUNT(n_same_module);
            else
              V2P2_COUNT(n_diff_module);
#endif
          }
#ifdef MKFIT_TRACE
          int pid = nc.m_trace_state_id;
          int id = mp_event->trace_new_cand_state(pid, (*mp_steeringparams_iter)->m_layer,
                                                  track2bivec3(nc), n.m_state);
          if (n.m_tr_hitmatch_id >= 0) {
            auto &ku = mp_event->tr_kalmanupdate( mp_event->tr_hitmatch(n.m_tr_hitmatch_id).kalman_id );
            ku.accepted = true;
            ku.state_id_out = id;
          }
          nc.m_trace_state_id = id;
#endif
        }
        nc.setState(m_sec_arena[e.node_idx].m_state);
        V2P2_COUNT(n_path_taken);
        if (n_chain > 1)
          V2P2_COUNT_ADD(n_extra_hits, n_chain - 1);
      } else if (e.add_fake) {
        switch (e.fake_hit) {
          case Hit::kHitMissIdx:  V2P2_COUNT(n_hole);       break;
          case Hit::kHitEdgeIdx:  V2P2_COUNT(n_hot_edge);   break;
          case Hit::kHitInGapIdx: V2P2_COUNT(n_hot_gap);    break;
          case Hit::kHitStopIdx:  V2P2_COUNT(n_stop_holes); break;
          default: break;
        }
        if (e.fake_hit == Hit::kHitStopIdx)
          nc.setStopReason(TrackCand::SR_TooManyHoles);
        nc.addHitIdx(e.fake_hit, m_rz_limits.layer_info_1().layer_id(), 0.0f);
#ifdef MKFIT_TRACE
        int pid = nc.m_trace_state_id;
        nc.m_trace_state_id = mp_event->trace_new_cand_state(
            pid, (*mp_steeringparams_iter)->m_layer, track2bivec3(nc), nc.state());
#endif
      }

      // score_ holds the layer-step score during finding; track_score_func
      // overwrites it at the end of the search.
      nc.setScore(e.score);
#ifdef MKFIT_TRACE
      {
        const short int lay_sel = m_rz_limits.layer_info_1().layer_id();
        const short int rank = k + 1;
        if (rank > nc.m_rank_max) {
          nc.m_rank_max = rank;
          nc.m_rank_max_layer = lay_sel;
        }
        if (rank > Diag::rank_over && nc.m_rank_over_layer < 0)
          nc.m_rank_over_layer = lay_sel;
      }
#endif
    }

    // Candidates that stopped IN this layer leave the beam the same way.
    ccand.clear();
    for (auto &c : m_new_cands) {
      if (best_short && c.getLastHitIdx() == Hit::kHitStopIdx) {
        offer_best_short(ccand, c);
        continue;
      }
      ccand.push_back(c);
    }

    V2P2_COUNT(n_selections);
    V2P2_COUNT_ADD(n_sel_entries, (long) m_sel.size());
    V2P2_COUNT_ADD(n_sel_kept, n_keep);
  }

} // namespace mkfit
