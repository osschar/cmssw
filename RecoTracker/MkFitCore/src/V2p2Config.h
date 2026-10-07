#ifndef RecoTracker_MkFitCore_src_V2p2Config_h
#define RecoTracker_MkFitCore_src_V2p2Config_h

#include "RecoTracker/MkFitCore/interface/Config.h"

//==============================================================================
// Tunable parameters of MkFinderV2p2 and of the v2p2 window in MkBins.
//
// Each knob is declared once, with its default, through MKFIT_V2P2_KNOB:
//   - CMSSW build:      constexpr;
//   - standalone build: a mutable variable, defined in V2p2Config.cc and set by
//                       mkFit.cc options and the shell drivers.
//
// The rationale for each default is in doc/MkFinderV2p2-DesignNotes.md, in the
// section named next to each group.
//==============================================================================

#if !defined(MKFIT_STANDALONE)
#define MKFIT_V2P2_KNOB(type, name, value) constexpr type name = value
#elif defined(MKFIT_V2P2_CONFIG_DEFINE)
#define MKFIT_V2P2_KNOB(type, name, value) type name = value
#else
#define MKFIT_V2P2_KNOB(type, name, value) extern type name
#endif

namespace mkfit {

  // Coefficients of the layer-step score (V2p2Score.h), one set per search
  // direction.
  struct V2p2ScoreParams {
    float hit_bonus = 3.0f;      // per hit taken; a hole wins above chi2 11, see the design notes
    float overlap_bonus = 0.0f;  // per hit beyond the first in a layer
    float chi2_weight = 1.0f;    // per unit of chi2
    float miss_penalty = 8.0f;   // a real hole
    float edge_penalty = 0.0f;   // crossing at the layer boundary
    float gap_penalty = 0.0f;    // crossing in an inactive region
    float stop_penalty = 8.0f;   // out of hole budget
    float hit_eff = 0.90f;       // per-layer hit efficiency, likelihood mode only
    // Per layer group, likelihood mode only; a negative value takes hit_eff. Groups, see
    // v2p2_hit_eff_group(): 0 pixel barrel, 1 pixel endcap, 2 strip barrel below layer 10
    // (phase-2 TBPS), 3 strip barrel from layer 10 (TB2S), 4 strip endcap (TEDD).
    float hit_eff_grp[5] = {-1.0f, -1.0f, -1.0f, -1.0f, -1.0f};
  };

  namespace Config::V2p2 {

    //--------------------------------------------------------------------------
    // Pre-selection window, design notes "Search window":
    //   ddq   < dq_trk_fac   * dq_trk     + dq_hit_fac   * hit_q_half_length
    //   ddphi < dphi_trk_fac * dphi_track + dphi_hit_fac * hit_phi_half_extent
    // dq_track and dphi_track are 3 sigma. The binnor fetch is derived from the
    // same factors, plus whole spare bins.
    namespace Window {
      // Field of the mini-propagators (window, edge states, Hermite path): 0 the uniform
      // Config::Bfield, 1 Config::bFieldFromZR() at each candidate's state, sampled once per
      // layer step at its start, as the Kalman propagation does with PF_use_param_b_field
      // (PropagationMPlexPlane.cc), 2 as 1 when the Kalman's finding_inter_layer_pflags has
      // use_param_b_field, else as 0. In TEDD the uniform field is 4-7 % above the parametrised.
      MKFIT_V2P2_KNOB(int, mini_prop_param_b, 2);
      MKFIT_V2P2_KNOB(float, dq_trk_fac, 1.5f);  // on dq_track
      MKFIT_V2P2_KNOB(float, dq_hit_fac, 1.2f);  // on hit_q_half_length; 1.0 just contains the hit
      MKFIT_V2P2_KNOB(bool, surface_q, true);    // reference dq_track to the hit's module plane
      MKFIT_V2P2_KNOB(int, q_extra_bins, 1);     // spare q bins in the fetch

      MKFIT_V2P2_KNOB(float, dphi_trk_fac, 2.0f);  // on dphi_track
      MKFIT_V2P2_KNOB(bool, phi_per_hit, true);    // hit term from the hit's own phi extent
      MKFIT_V2P2_KNOB(float, dphi_hit_fac, 3.0f);  // on hit_phi_half_extent
      MKFIT_V2P2_KNOB(float, dphi_flat_rad, 2.0f * 0.0123f);  // flat hit term [rad] when phi_per_hit is off
      MKFIT_V2P2_KNOB(int, phi_extra_bins, 0);     // spare phi bins in the fetch
    }  // namespace Window

    //--------------------------------------------------------------------------
    // Line pre-cut in MkFinderV2p2::select_hits(), design notes "Line pre-cut".
    // The slacks multiply the track terms of the real cut, so the pre-cut stays
    // looser than it.
    namespace PreCut {
      MKFIT_V2P2_KNOB(bool, q, true);
      MKFIT_V2P2_KNOB(bool, phi, true);
      MKFIT_V2P2_KNOB(float, dq_slack, 2.0f);
      MKFIT_V2P2_KNOB(float, dphi_slack, 1.5f);
      MKFIT_V2P2_KNOB(float, qbar_fac, 1.2f);  // on hit_qbar_half_extent, barrel only
    }  // namespace PreCut

    //--------------------------------------------------------------------------
    // In-layer combinatorial search and end-of-layer selection, design notes
    // "Reduction and hit ordering", "In-layer combinatorial search" and
    // "End-of-layer selection".
    namespace InLayer {
      MKFIT_V2P2_KNOB(bool, comb, true);         // off: one best hit per layer
      MKFIT_V2P2_KNOB(int, max_presel_hits, 6);  // reduction cap per sub-layer
      constexpr int max_presel_hits_limit = 12;   // inline storage; max_presel_hits must not exceed it
      MKFIT_V2P2_KNOB(int, max_sec_depth, 4);    // most hits one path may take in a layer
      constexpr int max_sec_depth_limit = 8;     // array size; max_sec_depth must not exceed it
      MKFIT_V2P2_KNOB(bool, reserve_hole_slot, false);  // keep one decliner in the beam
      MKFIT_V2P2_KNOB(bool, best_short, false);         // stopped candidates leave the beam (outward)
      // Select a CombCandidate as soon as its last candidate has been through a
      // layer batch, and return its tree nodes to the arena's free list, instead
      // of holding every node until end of layer.
      MKFIT_V2P2_KNOB(bool, early_select, true);
      // Keep the in-layer Kalman batches open across layer batches and fire them
      // only when full (or at end of layer), instead of flushing at every depth
      // boundary of every layer batch.
      MKFIT_V2P2_KNOB(bool, running_kalman, true);
      // Beam width per seed, decided when its CombCandidate is activated in a
      // layer: the CombCandidate keeps late_max_cands candidates instead of its
      // capacity (maxCandsPerSeed), and the surplus is dropped at activation,
      // unless its best candidate
      //  - has at most late_max_cands_hits found hits (the seed parameters are still
      //    poorly known), or
      //  - is more than late_max_cands_step_cm away from the layer, along a straight
      //    line from its last hit (a long step; 0 = this test off).
      // late_max_cands 0 = off.
      MKFIT_V2P2_KNOB(int, late_max_cands, 0);
      MKFIT_V2P2_KNOB(int, late_max_cands_hits, 10);
      MKFIT_V2P2_KNOB(float, late_max_cands_step_cm, 0.0f);
      // Sister-hole dominance: in the second-crossed sub-layer of an OT pair, the
      // hole of a TrackCand whose last hit is on the first sub-layer leaves the
      // selection when a path of the same TrackCand starts on the sister sensor
      // (detid + 1 of the lower sensor), with plain chi2 below sister_hole_chi2
      // (0 = any accepted hit).
      MKFIT_V2P2_KNOB(bool, sister_hole_drop, false);
      MKFIT_V2P2_KNOB(float, sister_hole_chi2, 0.0f);
    }  // namespace InLayer

    //--------------------------------------------------------------------------
    // Per-layer candidate policy and hit acceptance, design notes
    // "Within-sensitive-region verdict" and "Candidate pickup and stopping cuts".
    namespace Policy {
      MKFIT_V2P2_KNOB(bool, use_wsr, true);
      MKFIT_V2P2_KNOB(float, wsr_n_sigma, 5.0f);  // WSR edge fuzz, in sigma of the track's q error
      MKFIT_V2P2_KNOB(bool, use_hole_limits, true);  // maxHolesPerCand / maxConsecHoles
      MKFIT_V2P2_KNOB(bool, use_stop_cuts, true);    // minPtCut and the looper stop
      // The track scorer "v2p2:llh", for the final pick of each seed at the end of the
      // forward search (IterationConfig::m_final_pick_track_scorer_name): the summed
      // layer-step score of the search, minus final_pick_hole_penalty per inside hole
      // and final_pick_tail_penalty per tail hole.
      MKFIT_V2P2_KNOB(float, final_pick_hole_penalty, 8.0f);
      MKFIT_V2P2_KNOB(float, final_pick_tail_penalty, 3.0f);
      // Looper stop, forward search only: transverse angle between position and
      // momentum past looper_max_angle, for pT below looper_max_pt and r above
      // looper_min_r [cm].
      MKFIT_V2P2_KNOB(float, looper_max_angle, Const::PIOver2 - 0.2f);
      MKFIT_V2P2_KNOB(float, looper_max_pt, 1.2f);
      MKFIT_V2P2_KNOB(float, looper_min_r, 25.0f);
      // Kalman hit acceptance, both paths: chi2_acc < hit_chi2_cut, where chi2_acc
      // is the residual against chi2_trk_fac^2 * C_trk + C_hit in the module's
      // local frame. The hit covariance is taken as it is; chi2_trk_fac carries
      // the track-covariance deficit, so hit_chi2_cut is a probability: 9.21 is
      // the 99 % point of chi2 with 2 dof. The update and the score keep the
      // plain chi2. chi2_trk_fac = 1 is the plain chi2 and costs nothing extra;
      // the plain chi2 < 30 used before is hit_chi2_cut 30, chi2_trk_fac 1.
      MKFIT_V2P2_KNOB(float, hit_chi2_cut, 9.21f);
      MKFIT_V2P2_KNOB(float, chi2_trk_fac, 3.0f);
    }  // namespace Policy

    //--------------------------------------------------------------------------
    // Layer-step score, design notes "Layer-step score".
    namespace Score {
      MKFIT_V2P2_KNOB(int, mode, 1);  // 0 linear, 1 log-likelihood ratio
      // Forward: a phase-2 tuning (CMSSW MTV, initialStep with the offline seeds, cleaner and
      // refit, 2026-10): eps 0.10, TBPS 0.14, TB2S 0.20. It will move with the seeds and the
      // final fit. The initializer is parenthesised because its commas would split the macro.
      MKFIT_V2P2_KNOB(V2p2ScoreParams, fwd,
                      (V2p2ScoreParams{.hit_eff = 0.10f, .hit_eff_grp = {-1.0f, -1.0f, 0.14f, 0.20f, -1.0f}}));
      MKFIT_V2P2_KNOB(V2p2ScoreParams, bkw, {});
      // Term ablation of the likelihood: replace a term by a constant.
      MKFIT_V2P2_KNOB(bool, use_rho, true);
      MKFIT_V2P2_KNOB(float, rho_const, 0.0f);   // ln(hits/cm^2), used when use_rho is false
      MKFIT_V2P2_KNOB(bool, use_detv, true);
      MKFIT_V2P2_KNOB(float, detv_const, 0.0f);  // ln(det V) per hit, used when use_detv is false
      // The region rho counts hits over and divides by. 0: every hit of the fetched bins over
      // the track window's area (the count depends on the fetch: hit-extent padding, whole
      // bins, spare bins). 1: the track window for both, (n + 1) / area. 2: the cut region
      // for both, the window widened by the cut's track and largest-hit-extent terms as in
      // the fetch, without bin rounding or spare bins, (n + 1) / area.
      MKFIT_V2P2_KNOB(int, rho_region, 2);
    }  // namespace Score

#if defined(MKFIT_STANDALONE)
    //--------------------------------------------------------------------------
    // Diagnostic instruments, off by default.
    namespace Diag {
      // The MC-matched hit wins its layer and bypasses the chi2 cut. Best-hit
      // path, MKFIT_TRACE builds only. See PrimTCandRep::bKey.
      MKFIT_V2P2_KNOB(bool, force_mc, false);
      // Surface reference of dq_track with the layer's cylinder or disc normal,
      // in MkBins, so that the corrected value reaches the trace.
      MKFIT_V2P2_KNOB(bool, mkbins_surface_q, false);
      // Rank above which TrackCand::m_rank_over_layer records the layer, MKFIT_TRACE.
      MKFIT_V2P2_KNOB(int, rank_over, 6);
      // Standalone: classify each seed's beam at the end of the forward search by truth,
      // MkBuilder::v2p2_final_beam_diag_*() (single-threaded accumulation).
      MKFIT_V2P2_KNOB(bool, final_beam_purity, false);
      // Standalone: record each pair the pixel-seed duplicate cleaner decides, for
      // v2p2_dupclean_diag() to classify by truth (single-threaded accumulation).
      MKFIT_V2P2_KNOB(bool, dupclean_pairs, false);
    }  // namespace Diag

    // Running sums of the likelihood terms, for measuring the ablation
    // constants. Not thread safe: accumulate only in a serialised (trace) build.
    namespace ScoreStats {
      MKFIT_V2P2_KNOB(bool, accum, false);
      MKFIT_V2P2_KNOB(long, n_hits, 0);
      MKFIT_V2P2_KNOB(double, sum_log_rho, 0.0);  // weighted by n_hits
      MKFIT_V2P2_KNOB(double, sum_log_detv, 0.0);
    }  // namespace ScoreStats

    // Legacy compound dq setter (--v2p2-extra-dq, val_extra_dq): sets dq_trk_fac
    // to f and dq_hit_fac to 1.2 f, the ratio of the old single EXTRA_DQ factor.
    void set_extra_dq(float f);
#endif

  }  // namespace Config::V2p2

}  // namespace mkfit

#undef MKFIT_V2P2_KNOB

#endif
