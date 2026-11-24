#ifndef MkFitFWFF_Trace_TraceDigest_h
#define MkFitFWFF_Trace_TraceDigest_h

// A digest of the trace of one event: for each of the seven record vectors of a
// mkfit::TraceStore, its size and an order-dependent 64-bit FNV-1a hash of its
// contents, field by field, floats by their bits, padding skipped. Two traces
// have the same digest exactly when every record has the same values in the
// same order.
//
// The hashing is copied unchanged from mkFit-external's
// RecoTracker/MkFitCore/standalone/RdfTrace/TraceDigest.h (fwff-tracing
// afe5795), so a trace from cmsRun and one from the standalone are digested
// alike. Only the include (TraceStore.h in place of the standalone's
// RntStructs.h), the namespace and the printing differ. A field added to a
// record must be added to both.

#include "RecoTracker/MkFitCore/interface/TraceStore.h"

#include <array>
#include <cstdint>
#include <cstdio>
#include <type_traits>
#include <vector>

namespace mkfit_fwff::trace_digest {

  struct Hasher {
    static constexpr uint64_t kOffsetBasis = 1469598103934665603ull;
    static constexpr uint64_t kPrime = 1099511628211ull;

    uint64_t m_h = kOffsetBasis;

    void bytes(const void *p, size_t n) {
      const unsigned char *c = static_cast<const unsigned char *>(p);
      for (size_t i = 0; i < n; ++i)
        m_h = (m_h ^ c[i]) * kPrime;
    }

    // One field, by its object representation: a float by its bits.
    template <typename T>
    void operator()(T v) {
      static_assert(std::is_arithmetic_v<T>, "hash fields one by one, never a whole struct");
      bytes(&v, sizeof(T));
    }
  };

  inline void digest(Hasher &h, const EVec3 &v) {
    h(v.fX);
    h(v.fY);
    h(v.fZ);
  }

  inline void digest(Hasher &h, const EBiVec3 &v) {
    digest(h, v.pos);
    digest(h, v.mom);
  }

  inline void digest(Hasher &h, const PropInfo &p) {
    digest(h, static_cast<const EBiVec3 &>(p));
    h(p.dalpha);
    h(p.fail_flag);
  }

  inline void digest(Hasher &h, const mkfit::TrackState &s) {
    for (int i = 0; i < 6; ++i)
      h(s.parameters[i]);
    const float *e = s.errArray();
    for (unsigned int i = 0; i < mkfit::SMatrixSym66::rep_type::kSize; ++i)
      h(e[i]);
    h(s.charge);
    h(s.valid);
  }

  inline void digest(Hasher &h, const TrCandMeta &r) {
    h(r.id);
    h(r.event);
    h(r.seed);
    h(r.global_seed);
    h(r.sim);
    h(r.cand);
    for (int i = 0; i < 3; ++i)
      h(r.stage_ids[i]);
    h(r.fwd_rank_max);
    h(r.fwd_rank_max_layer);
    h(r.fwd_rank_over_layer);
    h(r.fwd_stop_reason);
    h(r.fwd_hits_hash);
  }

  inline void digest(Hasher &h, const TrCandStage &r) {
    h(r.id);
    h(r.meta_id);
    h(r.parent_stage_id);
    h(r.stage);
    h(r.iteration_idx);
    h(r.root_state_id);
    h(r.final_state_id);
  }

  inline void digest(Hasher &h, const TrCandState &r) {
    h(r.id);
    h(r.parent_id);
    h(r.meta_id);
    h(r.stage_id);
    h(r.layer);
    h(r.step);
    digest(h, r.kine);
    digest(h, r.state);
    h(r.search_id);
    h(r.has_children);
    h(r.on_final_path);
  }

  inline void digest(Hasher &h, const TrLayerSearch &r) {
    h(r.id);
    h(r.state_id);
    h(r.layer);
    h(r.layer_sec);
    h(r.is_barrel);
    h(r.is_outward);
    h(r.wsr);
    h(r.wsr_in_gap);
    digest(h, r.prop_entry);
    digest(h, r.prop_exit);
    h(r.phi_center);
    h(r.phi_delta);
    h(r.q_min);
    h(r.q_max);
    h(r.dphi_track);
    h(r.dq_track);
    h(r.cov_xx);
    h(r.cov_xy);
    h(r.cov_yy);
    h(r.cov_zz);
    h(r.p1);
    h(r.p2);
    h(r.q1);
    h(r.q2);
    h(r.p1_sec);
    h(r.p2_sec);
    h(r.q1_sec);
    h(r.q2_sec);
    h(r.n_hits_scanned);
    h(r.n_hits_masked);
    h(r.n_hits_precut);
    h(r.n_hits_presel);
    h(r.n_hits_pqueue);
    digest(h, r.hermite_mid_dev);
  }

  inline void digest(Hasher &h, const TrHitMatch &r) {
    h(r.id);
    h(r.state_id);
    h(r.search_id);
    h(r.layer);
    h(r.hit);
    h(r.mc_match);
    digest(h, r.kine_on_plane);
    h(r.t_hermite);
    h(r.d_plane_h3);
    h(r.dphi);
    h(r.dq);
    h(r.hit_q_half_len);
    h(r.score);
    h(r.passed_preselect);
    h(r.residual_x);
    h(r.residual_y);
    h(r.residual_z);
    h(r.is_sec_layer);
    h(r.sub_rank);
    h(r.full_rank);
    h(r.passed_pqueue);
    h(r.kalman_id);
#ifdef MKFIT_TRACE_PROP_COMPARE
    digest(h, r.kine_on_plane_cmp);
#endif
  }

  inline void digest(Hasher &h, const TrKalmanUpdate &r) {
    h(r.id);
    h(r.hit_match_id);
    h(r.state_id_in);
    h(r.state_id_out);
    h(r.chi2);
    h(r.chi2_trk);
    h(r.chi2_acc);
    h(r.accepted);
#ifdef MKFIT_TRACE_KALMAN_DEBUG
    digest(h, r.propagated_state);
    digest(h, r.updated_state);
#endif
  }

  inline void digest(Hasher &h, const TrBkFitUpdate &r) {
    h(r.id);
    h(r.state_id_in);
    h(r.state_id_out);
    h(r.layer);
    h(r.hit);
    h(r.mc_track_id);
    h(r.step);
    h(r.fail);
    h(r.chi2);
    h(r.chi2_cum);
    digest(h, r.kine_on_plane);
    h(r.residual_x);
    h(r.residual_y);
    h(r.residual_z);
#ifdef MKFIT_TRACE_KALMAN_DEBUG
    digest(h, r.propagated_state);
#endif
  }

  struct VecDigest {
    const char *name;
    size_t size;
    uint64_t hash;
  };

  template <typename Rec>
  VecDigest digest_vector(const char *name, const std::vector<Rec> &v) {
    Hasher h;
    for (const Rec &r : v)
      digest(h, r);
    return {name, v.size(), h.m_h};
  }

  template <typename Store>
  std::array<VecDigest, 7> digest_store(const Store &s) {
    return {{digest_vector("trCandMetas_", s.trCandMetas_),
             digest_vector("trCandStages_", s.trCandStages_),
             digest_vector("trCandStates_", s.trCandStates_),
             digest_vector("trLayerSearches_", s.trLayerSearches_),
             digest_vector("trHitMatches_", s.trHitMatches_),
             digest_vector("trKalmanUpdates_", s.trKalmanUpdates_),
             digest_vector("trBkFitUpdates_", s.trBkFitUpdates_)}};
  }

  // The record vectors in store order, by their names in TraceStore.
  constexpr std::array<const char *, 7> kVectorNames{{"trCandMetas_",
                                                      "trCandStages_",
                                                      "trCandStates_",
                                                      "trLayerSearches_",
                                                      "trHitMatches_",
                                                      "trKalmanUpdates_",
                                                      "trBkFitUpdates_"}};

}  // namespace mkfit_fwff::trace_digest

#endif
