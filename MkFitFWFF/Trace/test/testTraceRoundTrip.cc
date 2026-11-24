// testMkFitFWFFTraceRoundTrip [file]
//
// Writes synthetic trace stores with TraceRNTupleWriter (three events, one of them twice, one empty), reads the file
// back with readTraceFile() and compares each event's digest from the file with the digest of the store in memory.
// Every field of every record type gets a value that differs between records, so a field the I/O dropped or mangled
// changes a hash. Exit 0 when all digests and keys agree.

#include "MkFitFWFF/Trace/interface/TraceDigest.h"
#include "MkFitFWFF/Trace/interface/TraceFileDigest.h"
#include "MkFitFWFF/Trace/interface/TraceRNTupleWriter.h"
#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/interface/TraceStore.h"

#include <cstdio>
#include <string>
#include <vector>

namespace {
  float f(int i, int k) { return 0.001f * i + 0.37f * k - 1.5f; }

  EVec3 v3(int i, int k) { return EVec3(f(i, k), f(i, k + 1), f(i, k + 2)); }

  mkfit::TrackState ts(int i) {
    mkfit::TrackState s;
    for (int p = 0; p < 6; ++p)
      s.parameters[p] = f(i, p);
    float* e = s.errArray_nc();
    for (int p = 0; p < 21; ++p)
      e[p] = f(i, 10 + p);
    s.charge = (i % 2) ? 1 : -1;
    s.valid = (i % 3) != 0;
    return s;
  }

  void fill(mkfit::TraceStore& s, int n, int seed) {
    for (int i = 0; i < n; ++i) {
      const int j = i + 100 * seed;
      TrCandMeta& m = s.trace_candmeta({-1, seed, i});
      m.global_seed = 2 * j;
      m.sim = j - 1;
      m.cand = (i % 2) ? i / 2 : -1;
      m.stage_ids[0] = j;
      m.stage_ids[1] = j + 1;
      m.stage_ids[2] = -j;
      m.fwd_rank_max = j % 7;
      m.fwd_rank_max_layer = j % 11;
      m.fwd_rank_over_layer = -(j % 5);
      m.fwd_stop_reason = j % 3;
      m.fwd_hits_hash = 0x9e3779b97f4a7c15ull * (j + 1);

      auto [stage_id, state_id] =
          s.trace_new_cand_stage_and_state(m.id, -1, i % 3, j % 30, EBiVec3{v3(j, 0), v3(j, 3)}, ts(j));
      s.tr_candstage(stage_id).iteration_idx = j % 4;
      s.tr_candstage(stage_id).final_state_id = state_id;
      const int child = s.trace_new_cand_state(state_id, j % 29, EBiVec3{v3(j, 6), v3(j, 9)}, ts(j + 1));
      s.tr_candstate(child).search_id = j;
      s.tr_candstate(child).on_final_path = (i % 2) == 0;

      TrLayerSearch ls;
      ls.state_id = state_id;
      ls.layer = j % 50;
      ls.layer_sec = j % 3 - 1;
      ls.is_barrel = i % 2;
      ls.is_outward = i % 3;
      ls.wsr = (signed char)(j % 4 - 1);
      ls.wsr_in_gap = i % 5 == 0;
      static_cast<EBiVec3&>(ls.prop_entry) = EBiVec3{v3(j, 12), v3(j, 15)};
      ls.prop_entry.dalpha = f(j, 18);
      ls.prop_entry.fail_flag = j % 2;
      static_cast<EBiVec3&>(ls.prop_exit) = EBiVec3{v3(j, 19), v3(j, 22)};
      ls.prop_exit.dalpha = f(j, 25);
      ls.prop_exit.fail_flag = j % 3;
      ls.phi_center = f(j, 26);
      ls.phi_delta = f(j, 27);
      ls.q_min = f(j, 28);
      ls.q_max = f(j, 29);
      ls.dphi_track = f(j, 30);
      ls.dq_track = f(j, 31);
      ls.cov_xx = f(j, 32);
      ls.cov_xy = f(j, 33);
      ls.cov_yy = f(j, 34);
      ls.cov_zz = f(j, 35);
      ls.p1 = j % 100;
      ls.p2 = j % 100 + 1;
      ls.q1 = j % 13;
      ls.q2 = j % 13 + 2;
      ls.p1_sec = j % 17;
      ls.p2_sec = j % 17 + 3;
      ls.q1_sec = j % 19;
      ls.q2_sec = j % 19 + 4;
      ls.n_hits_scanned = j;
      ls.n_hits_masked = j / 2;
      ls.n_hits_precut = j / 3;
      ls.n_hits_presel = j / 4;
      ls.n_hits_pqueue = j / 5;
      ls.hermite_mid_dev = v3(j, 36);
      const int ls_id = s.trace_layersearch(std::move(ls)).id;

      TrHitMatch hm{-1, child, ls_id, j % 40, j * 3, (i % 2) == 1};
      hm.kine_on_plane = EBiVec3{v3(j, 40), v3(j, 43)};
      hm.t_hermite = f(j, 46);
      hm.d_plane_h3 = f(j, 47);
      hm.dphi = f(j, 48);
      hm.dq = f(j, 49);
      hm.hit_q_half_len = f(j, 50);
      hm.score = f(j, 51);
      hm.passed_preselect = i % 2;
      hm.residual_x = f(j, 52);
      hm.residual_y = f(j, 53);
      hm.residual_z = f(j, 54);
      hm.is_sec_layer = i % 3 == 1;
      hm.sub_rank = j % 6;
      hm.full_rank = j % 7;
      hm.passed_pqueue = i % 4 == 0;
      TrHitMatch& hmr = s.trace_hitmatch(std::move(hm));
      const int ku_id = s.trace_new_kalman_update(hmr.id, state_id, f(j, 55), f(j, 56));
      hmr.kalman_id = ku_id;
      s.tr_kalmanupdate(ku_id).chi2_acc = f(j, 57);
      s.tr_kalmanupdate(ku_id).accepted = i % 2;
      s.tr_kalmanupdate(ku_id).state_id_out = child;

      TrBkFitUpdate bu{-1, state_id, child, j % 30, j * 5, j - 7, i, j % 2};
      bu.chi2 = f(j, 58);
      bu.chi2_cum = f(j, 59);
      bu.kine_on_plane = EBiVec3{v3(j, 60), v3(j, 63)};
      bu.residual_x = f(j, 66);
      bu.residual_y = f(j, 67);
      bu.residual_z = f(j, 68);
      s.trace_bkfitupdate(std::move(bu));
    }
  }
}  // namespace

int main(int argc, char* argv[]) {
  const std::string file = argc > 1 ? argv[1] : "testMkFitFWFFTraceRoundTrip.root";
  using namespace mkfit_fwff;

  struct Ev {
    unsigned run, lumi;
    unsigned long long event;
    std::string module;
    int n;
  };
  const std::vector<Ev> evs = {{1, 3, 201, "initialStepTrackCandidatesMkFit", 50},
                               {1, 3, 202, "initialStepTrackCandidatesMkFit", 0},
                               {1, 3, 201, "initialStepTrackCandidatesMkFit", 50},
                               {1, 3, 201, "hltInitialStepTrackCandidatesMkFit", 7}};

  std::vector<std::string> memLines;
  {
    TraceRNTupleWriter w(file);
    for (const auto& e : evs) {
      mkfit::TraceStore s;
      fill(s, e.n, e.n);
      std::vector<mkfit::Track> tracks(3);
      for (int t = 0; t < 3; ++t)
        tracks[t].setLabel(2 * (100 * e.n + 2 * t + 1));
      mkfit::TraceEventInfo info;
      info.run = e.run;
      info.lumi = e.lumi;
      info.event = e.event;
      info.module = e.module;
      info.n_seeds_in = e.n + 1;
      info.out_tracks = &tracks;
      const int pass = w.write(s, info);
      for (auto const& d : trace_digest::digest_store(s))
        memLines.push_back(digestLine(e.run, e.lumi, e.event, e.module, pass, d));
    }
    w.close();
  }

  const auto fevs = readTraceFile(file);
  std::vector<std::string> fileLines;
  int bad = 0;
  for (const auto& fe : fevs) {
    for (int r = 0; r < 7; ++r) {
      fileLines.push_back(digestLine(fe.run, fe.lumi, fe.event, fe.module, fe.pass, fe.digest[r]));
      if ((std::size_t)fe.n_rec[r] != fe.digest[r].size)
        ++bad;
    }
    std::printf("event %u:%u:%llu %s pass=%d metas=%zu exported=%d out_matched=%d/%zu\n",
                fe.run,
                fe.lumi,
                (unsigned long long)fe.event,
                fe.module.c_str(),
                fe.pass,
                fe.digest[0].size,
                fe.n_exported,
                fe.n_out_matched,
                fe.out_labels.size());
  }
  if (fileLines.size() != memLines.size())
    ++bad;
  for (std::size_t i = 0; i < memLines.size() && i < fileLines.size(); ++i) {
    const bool eq = memLines[i] == fileLines[i];
    bad += !eq;
    std::printf("%s %s\n", eq ? "same" : "DIFF", memLines[i].c_str());
    if (!eq)
      std::printf("     file: %s\n", fileLines[i].c_str());
  }
  std::printf("%s: %d mismatches over %zu digest lines\n", bad ? "FAILED" : "PASSED", bad, memLines.size());
  return bad ? 1 : 0;
}
