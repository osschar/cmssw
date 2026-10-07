// The backward fit of MkFinderV2p2 tracks. A member of MkFinder, as it runs on
// MkFinder's Matriplex state and on bkFitInputTracks() / bkFitOutputTracks(),
// but kept here so that MkFinder::bkFitFitTracksProp2Plane(), used by V1 / V2,
// stays as it is upstream. MkBuilder::fit_cands() picks this one when the
// candidates come from MkBuilder::findTracksStandardv2p2().
//
// It walks the hits on track backwards with propagation to each hit's module
// plane and needs no layer plan.

#include "MkFinder.h"

#include "RecoTracker/MkFitCore/interface/cms_common_macros.h"
#include "RecoTracker/MkFitCore/interface/IterationConfig.h"
#include "KalmanUtilsMPlex.h"
#include "V2p2Config.h"

#include <algorithm>

//#define DEBUG
#include "Debug.h"

#if defined(MKFIT_STANDALONE)
#include "RecoTracker/MkFitCore/standalone/Event.h"
#endif

#ifdef MKFIT_TRACE
#include "RecoTracker/MkFitCore/standalone/DataFormats/RntStructs.h"
#include "RecoTracker/MkFitCore/standalone/DataFormats/RntConversions.h"
#endif

namespace mkfit {

  void MkFinder::bkFitFitTracksV2p2(const EventOfHits &eventofhits,
                                    const SteeringParams &st_par,
                                    const int N_proc,
                                    bool chiDebug) {
    // Prototyping final backward fit.
    // This works with track-finding indices, before remapping.
    //
    // Layers should be collected during track finding and list all layers that have actual hits.
    // Then we could avoid checking which layers actually do have hits.

    // bool debug = true;

    MPlexQF tmp_chi2{0.0f};
    MPlexQI done_flag(0);
    // The last hit fitted per lane, to skip a repeated hit on track (seed-hit
    // merging can put the same hit on twice).
    int fit_lay[NN], fit_idx[NN];
    for (int i = 0; i < NN; ++i) fit_lay[i] = fit_idx[i] = -1;

#ifdef MKFIT_TRACE
    // Per-lane bookkeeping for the trace: which hit each lane is on this pass
    // and how many steps it has taken. These existed for the old BkFitHook,
    // which is gone, so they are now trace-only and live inside the guard --
    // outside it they are set and never read, which is a -Werror build failure.
    int hk_layer[NN], hk_mcid[NN], hk_hit[NN], hk_step[NN], hk_live[NN];
    for (int i = 0; i < NN; ++i) { hk_layer[i] = hk_mcid[i] = hk_hit[i] = -1; hk_step[i] = 0; hk_live[i] = 0; }

    // Trace chain: one TrCandState per lane, advanced at every accepted hit.
    // Kept LOCAL rather than on TrackCand::m_trace_state_id, which belongs to
    // the search -- writing it here would make the search chain its own states
    // onto the fit's and inherit the wrong stage_id.
    int tr_state[NN];
    for (int i = 0; i < NN; ++i) tr_state[i] = -1;
    const bool bk_trace = (m_event != nullptr);
    if (bk_trace) {
      for (int i = 0; i < N_proc; ++i) {
        TrackCand *tc = m_TrkCand[i];
        if (tc == nullptr) continue;
        CombCandidate *cc = tc->combCandidate();
        if (cc == nullptr) continue;
        // Lazily created and idempotent, exactly as MkBuilder does it -- in HLT
        // the backward fit runs BEFORE any search, so it may well be first here.
        if (cc->m_trace_meta_id == -1)
          cc->m_trace_meta_id = m_event->trace_new_cand_meta(m_event->evtID(), cc->seed_origin_index());
        // The root state is the INPUT state, i.e. AFTER bkFitInputTracks has
        // inflated the covariance by 100 (10x in sigma).
        // Do not compare it against a search state without allowing for that.
        TrackState ts;
        m_Par[iC].copyOut(i, ts.parArray_nc());
        m_Err[iC].copyOut(i, ts.errArray_nc());
        ts.charge = m_Chg[i];
        const EBiVec3 kine { EVec3(ts.x(), ts.y(), ts.z()), EVec3(ts.px(), ts.py(), ts.pz()) };
        // layer -1: the root sits before any hit. stage 1 = BkwFit.
        auto [stage_id, state_id] =
            m_event->trace_new_cand_stage_and_state(cc->m_trace_meta_id, -1, 1, -1, kine, ts);
        cc->m_trace_stage_id = stage_id;
        TrCandMeta &cm = m_event->tr_candmeta(cc->m_trace_meta_id);
        cm.stage_ids[1] = stage_id;
        // MkBuilder fills this too, but in HLT the fit runs FIRST, so the meta
        // can exist without it. It is the truth-join key: seeds are relabelled
        // sequentially, so a track's label IS its index in seedTracks_.
        if (cm.global_seed == -1 && tc->label() >= 0)
          cm.global_seed = tc->label();
        tr_state[i] = state_id;
      }
    }
    const bool bk_rec = bk_trace;
#endif

    MPlexHV plNrm{0.0f};  // input detector plane [pl - plane]
    MPlexHV plDir{0.0f};  // ""
    MPlexHV plPnt{0.0f};  // ""

#if defined(DEBUG_PROP_UPDATE)
    const int DSLOT = 0;
    int DSLOT_layer;
    printf("bkfit-p2p entry, track in slot %d\n", DSLOT);
    print_par_err(iC, DSLOT);
#endif
#if defined(DEBUG_BACKWARD_FIT)
    const Hit *last_hit_ptr[NN];
    int last_layer[NN];
#endif

    // Skip the last hit (or two), ie, do not refit it (them).
    // If there are overlap hits in the same layer, they will still get processed.
    // A more proper thing to do might be to:
    // a) skip all hits on the last double layer; or
    // b) skip last hits that are closer than some ds, say, 5 cm.
    // for (int i = 0; i < N_proc; ++i) {
    //   m_CurNode[i] = m_HoTNodeArr[i][m_CurNode[i]].m_prev_idx;
    //   // m_CurNode[i] = m_HoTNodeArr[i][m_CurNode[i]].m_prev_idx;
    // }

    // Outlier rejection, as MkFinder::bkFitFitTracksProp2Plane(): the IterationConfig's
    // m_backward_fit_outlier_* (chi2 0 = off).
    const float outlier_chi2 = m_iteration_config->m_backward_fit_outlier_chi2;
    const int max_outliers = m_iteration_config->m_backward_fit_max_outliers;
    const float outlier_min_pt = m_iteration_config->m_backward_fit_outlier_min_pt;
    int n_outliers[NN] = {0};

    int done_count = 0;
    while (done_count != N_proc) {

      int fit_node[NN];  // HoT node fitted in this step, -1 if none
      std::fill_n(fit_node, NN, -1);
      int here_count = 0;
      for (int i = 0; i < N_proc; ++i) {
        if (done_flag[i])
          continue;

        // skip holes and a repeat of the hit just fitted
        while (m_CurNode[i] >= 0 && (m_HoTNodeArr[i][m_CurNode[i]].m_hot.index < 0 ||
                                     (m_HoTNodeArr[i][m_CurNode[i]].m_hot.layer == fit_lay[i] &&
                                      m_HoTNodeArr[i][m_CurNode[i]].m_hot.index == fit_idx[i]))) {
          m_CurNode[i] = m_HoTNodeArr[i][m_CurNode[i]].m_prev_idx;
        }

        if (m_CurNode[i] < 0) {
          // Mark as done and copy out.
          done_flag[i] = 1;
          ++done_count;

          TrackCand &trk = *m_TrkCand[i];
          m_Err[iC].copyOut(i, trk.errors_nc().Array());
          m_Par[iC].copyOut(i, trk.parameters_nc().Array());
          trk.setCharge(m_Chg[i]);
          trk.setChi2(m_Chi2[i]);
          if (isFinite(trk.chi2())) {
            trk.setScore(getScoreCand(m_steering_params->m_post_bkfit_track_scorer, trk));
          }
        } else {
          // Prepare the next hit and module info.

          // Every hit is fitted, overlap and sister hits included: MkFinderV2p2
          // puts a layer's hits on track in path order, so walking the chain
          // backwards is walking the path backwards. (The V1 / V2 fit keeps only
          // the first of consecutive hits in a layer, as their overlap hit is
          // placed after the primary one rather than in path order.)
          const int layer = m_HoTNodeArr[i][m_CurNode[i]].m_hot.layer;
          fit_lay[i] = layer;
          fit_idx[i] = m_HoTNodeArr[i][m_CurNode[i]].m_hot.index;

          const LayerOfHits &L = eventofhits[layer];
          const Hit &hit = L.refHit(m_HoTNodeArr[i][m_CurNode[i]].m_hot.index);
          const ModuleInfo &mi = L.layer_info().module_info(hit.detIDinLayer());

          m_msErr.copyIn(i, hit.errArray());
          m_msPar.copyIn(i, hit.posArray());
          plNrm.copyIn(i, mi.zdir.Array());
          plDir.copyIn(i, mi.xdir.Array());
          plPnt.copyIn(i, mi.pos.Array());

#ifdef MKFIT_TRACE
          if (bk_rec) {
            hk_layer[i] = layer;
            hk_mcid[i] = hit.mcHitID();
            hk_hit[i] = m_HoTNodeArr[i][m_CurNode[i]].m_hot.index;
            hk_live[i] = 1;
          }
#endif

          fit_node[i] = m_CurNode[i];
          ++here_count;

          m_CurNode[i] = m_HoTNodeArr[i][m_CurNode[i]].m_prev_idx;

#ifdef DEBUG_BACKWARD_FIT
          last_hit_ptr[i] = &hit;
          last_layer[i] = layer;
#endif
#if defined(DEBUG_PROP_UPDATE)
          DSLOT_layer = layer;
          printf("\nbkfit start layer %d, track in slot %d -- fail=%d, hit_xyz = (%g, %g, %g)\n\n",
             DSLOT_layer, DSLOT, m_FailFlag[DSLOT],
             m_msPar(DSLOT, 0, 0), m_msPar(DSLOT, 1, 0), m_msPar(DSLOT, 2, 0));
#endif
        }
      }

      if (done_count == N_proc)
        break;
      if (here_count == 0)
        continue;

      // ZZZ Could add missing hits here, only if there are any actual matches.

      clearFailFlag();

      // PROP-FAIL-ENABLE We do not check for pfailed propagation here
      // clang-format off

      m_FailFlag.setVal(0);
      const MPlexQI chg_prev = m_Chg;
      propagateHelixToPlaneMPlex(m_Err[iC], m_Par[iC], m_Chg, plPnt, plNrm,
                                 m_Err[iP], m_Par[iP], m_FailFlag,
                                 N_proc, m_prop_config->backward_fit_pflags, nullptr);
      kalmanOperationPlaneLocal(KFO_Calculate_Chi2 | KFO_Update_Params | KFO_Local_Cov,
                                m_Err[iP], m_Par[iP], m_Chg, m_msErr, m_msPar, plNrm, plDir, plPnt,
                                m_Err[iC], m_Par[iC], tmp_chi2, N_proc);
      kalmanCheckChargeFlip(m_Par[iC], m_Chg, N_proc);

      // Outlier rejection: the hit becomes a missing hit and
      // the propagated state is kept. The negated comparison leaves a NaN chi2 alone.
      if (outlier_chi2 > 0.0f) {
        for (int i = 0; i < N_proc; ++i) {
          if (fit_node[i] < 0 || !(tmp_chi2[i] > outlier_chi2) ||
              n_outliers[i] >= max_outliers || m_TrkCand[i]->pT() < outlier_min_pt)
            continue;
          m_Err[iC].copySlot(i, m_Err[iP]);
          m_Par[iC].copySlot(i, m_Par[iP]);
          m_Chg[i] = chg_prev[i];
          tmp_chi2[i] = 0.0f;
          TrackCand &trk = *m_TrkCand[i];
          trk.combCandidate()->hot_node_nc(fit_node[i]).m_hot.index = Hit::kHitMissIdx;
          trk.setNFoundHits(trk.nFoundHits() - 1);
          trk.setNMissingHits(trk.nMissingHits() + 1);
          ++n_outliers[i];
        }
      }

#ifdef MKFIT_TRACE
      // The whole block records TrBkFitUpdate and nothing else, so it is guarded
      // as a unit. It used to be compiled unconditionally with bk_rec forced
      // false, which left its locals unused and broke a non-tracing -Werror build.
      if (bk_rec) {
        for (int i = 0; i < N_proc; ++i) {
          if (!hk_live[i]) continue;
          hk_live[i] = 0;
          const int step = hk_step[i]++;
          // residual of the PROPAGATED state (iP) against the hit, in the
          // module frame. ydir = zdir x xdir.
          const float dx = m_Par[iP].constAt(i, 0, 0) - m_msPar.constAt(i, 0, 0);
          const float dy = m_Par[iP].constAt(i, 1, 0) - m_msPar.constAt(i, 1, 0);
          const float dz = m_Par[iP].constAt(i, 2, 0) - m_msPar.constAt(i, 2, 0);
          const float nx = plNrm.constAt(i, 0, 0), ny = plNrm.constAt(i, 1, 0), nz = plNrm.constAt(i, 2, 0);
          const float ux = plDir.constAt(i, 0, 0), uy = plDir.constAt(i, 1, 0), uz = plDir.constAt(i, 2, 0);
          const float vx = ny * uz - nz * uy, vy = nz * ux - nx * uz, vz = nx * uy - ny * ux;
          const float d_xdir = dx * ux + dy * uy + dz * uz;
          const float d_ydir = dx * vx + dy * vy + dz * vz;
          const float d_zdir = dx * nx + dy * ny + dz * nz;
          const float ipt = m_Par[iP].constAt(i, 3, 0);
          const float pt = ipt != 0.f ? 1.f / std::abs(ipt) : 0.f;
          const float theta = m_Par[iP].constAt(i, 5, 0);

          if (bk_trace && tr_state[i] >= 0) {
            const float phi_p = m_Par[iP].constAt(i, 4, 0);
            TrBkFitUpdate bu;
            bu.state_id_in = tr_state[i];
            bu.layer = hk_layer[i];
            bu.hit = hk_hit[i];
            bu.step = step;
            bu.fail = m_FailFlag[i];
            bu.chi2 = tmp_chi2.constAt(i, 0, 0);
            bu.chi2_cum = m_Chi2.constAt(i, 0, 0);
            bu.kine_on_plane = { EVec3(m_Par[iP].constAt(i, 0, 0),
                                       m_Par[iP].constAt(i, 1, 0),
                                       m_Par[iP].constAt(i, 2, 0)),
                                 EVec3(pt * std::cos(phi_p), pt * std::sin(phi_p),
                                       pt / std::tan(theta)) };
            bu.residual_x = d_xdir;
            bu.residual_y = d_ydir;
            bu.residual_z = d_zdir;
            // WHICH particle made this hit. Not compared to anything here: the
            // track's sim label is not available at this point (see the struct).
            if (hk_mcid[i] >= 0 && hk_mcid[i] < (int) m_event->simHitsInfo_.size())
              bu.mc_track_id = m_event->simHitsInfo_[hk_mcid[i]].mcTrackID();
            // Post-update state: iC, after kalmanOperationPlaneLocal above.
            TrackState ts;
            m_Par[iC].copyOut(i, ts.parArray_nc());
            m_Err[iC].copyOut(i, ts.errArray_nc());
            ts.charge = m_Chg[i];
#ifdef MKFIT_TRACE_KALMAN_DEBUG
            TrackState tp;
            m_Par[iP].copyOut(i, tp.parArray_nc());
            m_Err[iP].copyOut(i, tp.errArray_nc());
            tp.charge = m_Chg[i];
            bu.propagated_state = tp;
#endif
            const EBiVec3 kine { EVec3(ts.x(), ts.y(), ts.z()), EVec3(ts.px(), ts.py(), ts.pz()) };
            const int sid = m_event->trace_new_cand_state(tr_state[i], hk_layer[i], kine, ts);
            bu.state_id_out = sid;
            tr_state[i] = sid;
            m_event->trace_bkfitupdate(std::move(bu));
          }
        }
      }
#endif

#if defined(DEBUG_PROP_UPDATE)
      printf("\nbkfit at layer %d, track in slot %d -- fail=%d, hit_xyz = (%g, %g, %g)\n",
             DSLOT_layer, DSLOT, m_FailFlag[DSLOT],
             m_msPar(DSLOT, 0, 0), m_msPar(DSLOT, 1, 0), m_msPar(DSLOT, 2, 0));
      printf("Propagated:\n");
      print_par_err(iP, DSLOT);
      printf("Updated:\n");
      print_par_err(iC, DSLOT);
#endif
      // clang-format on

      // Fixup for failed propagation.
      // for (int i = 0; i < NN; ++i) {
      // PROP-FAIL-ENABLE The following to be enabled when propagation failure
      // detection is properly implemented in propagate-to-R/Z.
      // 1. The following code was only expecting barrel state to be restored.
      //      auto barrel_pf(m_prop_config->backward_fit_pflags);
      //      barrel_pf.copy_input_state_on_fail = true;
      // 2. There is also check on chi2, commented out to keep physics changes minimal.
      /*
        if (m_FailFlag[i] && LI.is_barrel()) {
          // Barrel pflags are set to include PF_copy_input_state_on_fail.
          // Endcap errors are immaterial here (relevant for fwd search), with prop error codes
          // one could do other things.
          // Are there also fail conditions in KalmanUpdate?
#ifdef DEBUG
          if (debug && g_debug) {
            dprintf("MkFinder::bkFitFitTracks prop fail: chi2=%f, layer=%d, label=%d. Recovering.\n",
                    tmp_chi2[i], LI.layer_id(), m_Label[i]);
            print_par_err(iC, i);
          }
#endif
          m_Err[iC].copySlot(i, m_Err[iP]);
          m_Par[iC].copySlot(i, m_Par[iP]);
        } else if (tmp_chi2[i] > 200 || tmp_chi2[i] < 0) {
#ifdef DEBUG
          if (debug && g_debug) {
            dprintf("MkFinder::bkFitFitTracks chi2 fail: chi2=%f, layer=%d, label=%d. Recovering.\n",
                    tmp_chi2[i], LI.layer_id(), m_Label[i]);
            print_par_err(iC, i);
          }
#endif
          // Go back to propagated state (at the current hit, the previous one is lost).
          m_Err[iC].copySlot(i, m_Err[iP]);
          m_Par[iC].copySlot(i, m_Par[iP]);
        }
        */
      // }

#if defined(DEBUG_BACKWARD_FIT)
      // clang-format off
      bool debug = true;
      const char beg_cur_sep = '/'; // set to ' ' root parsable printouts
      for (int i = 0; i < N_proc; ++i) {
        if (chiDebug && last_hit_ptr[i]) {
          TrackCand &bb = *m_TrkCand[i];
          int ti = iP;
          float chi = tmp_chi2.At(i, 0, 0);
          float chi_prnt = std::isfinite(chi) ? chi : -9;
          const int layer = last_layer[i];
          const LayerOfHits &L = eventofhits[layer];

#if defined(MKFIT_STANDALONE)
          const MCHitInfo &mchi = m_event->simHitsInfo_[last_hit_ptr[i]->mcHitID()];

          dprintf("BKF_OVERLAP %d %d %d %d %d %d %d "
                  "%f%c%f %f %f%c%f %f %f %f %d %d %d %d "
                  "%f %f %f %f %f\n",
              m_event->evtID(),
#else
          dprintf("BKF_OVERLAP %d %d %d %d %d %d "
                  "%f%c%f %f %f%c%f %f %f %f %d %d %d "
                  "%f %f %f %f %f\n",
#endif
              bb.label(), (int)bb.prodType(), bb.isFindable(),
              layer, L.is_stereo(), L.is_barrel(),
              bb.pT(), beg_cur_sep, 1.0f / m_Par[ti].At(i, 3, 0),
              bb.posEta(),
              bb.posPhi(), beg_cur_sep, std::atan2(m_Par[ti].At(i, 1, 0), m_Par[ti].At(i, 0, 0)),
              hipo(m_Par[ti].At(i, 0, 0), m_Par[ti].At(i, 1, 0)),
              m_Par[ti].At(i, 2, 0),
              chi_prnt,
              std::isnan(chi), std::isfinite(chi), chi > 0,
#if defined(MKFIT_STANDALONE)
              mchi.mcTrackID(),
#endif
              // The following three can get negative / prouce nans in e2s.
              // std::abs the args for FPE hunt.
              e2s(std::abs(m_Err[ti].At(i, 0, 0))),
              e2s(std::abs(m_Err[ti].At(i, 1, 1))),
              e2s(std::abs(m_Err[ti].At(i, 2, 2))),  // sx_t sy_t sz_t -- track errors
              1e4f * hipo(m_msPar.At(i, 0, 0) - m_Par[ti].At(i, 0, 0),
                                m_msPar.At(i, 1, 0) - m_Par[ti].At(i, 1, 0)),  // d_xy
              1e4f * (m_msPar.At(i, 2, 0) - m_Par[ti].At(i, 2, 0))             // d_z
          );
        }
      }
      // clang-format on
#endif

      // update chi2
      m_Chi2.add(tmp_chi2);
    }
  }

}  // namespace mkfit
