// mkFitTraceDigest <trace file> [--no-check]
//
// Reads a trace file of MkFitTraceWriter and prints, per traced event in the order written, the digest line of each
// record vector (as MkFitTraceWriter prints them from memory with printDigest) and one line of counts:
//   MKFITTRACE check <run>:<lumi>:<event> <module> pass=<p> seeds_in= seeds= tracks= metas= exported= cands=
//              labeled= out_matched= stages=<fwd>/<bkfit>/<bkw> final=<fwd>/<bkfit>/<bkw> on_final_path= counts_ok=
// exported: TrCandMeta with cand >= 0; cands: their distinct cand values; out_matched: output tracks whose label is the
// global_seed of an exported TrCandMeta; counts_ok: the record counts in Events equal the rows read.
// Exit 0 if it read the file, 1 otherwise.

#include "MkFitFWFF/Trace/interface/TraceFileDigest.h"

#include <cstdio>
#include <cstring>
#include <exception>

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <trace file> [--no-check]\n", argv[0]);
    return 1;
  }
  const bool check = !(argc > 2 && std::strcmp(argv[2], "--no-check") == 0);
  try {
    const auto evs = mkfit_fwff::readTraceFile(argv[1]);
    std::printf("# %s: %zu traced events\n", argv[1], evs.size());
    for (const auto& ev : evs) {
      bool counts_ok = true;
      for (int r = 0; r < 7; ++r) {
        std::printf("%s\n",
                    mkfit_fwff::digestLine(ev.run, ev.lumi, ev.event, ev.module, ev.pass, ev.digest[r]).c_str());
        counts_ok = counts_ok && (std::size_t)ev.n_rec[r] == ev.digest[r].size;
      }
      if (check)
        std::printf(
            "MKFITTRACE check %u:%u:%llu %s pass=%d seeds_in=%d seeds=%d tracks=%d metas=%zu exported=%d cands=%d "
            "labeled=%d out_matched=%d stages=%d/%d/%d final=%d/%d/%d on_final_path=%d counts_ok=%d\n",
            ev.run,
            ev.lumi,
            (unsigned long long)ev.event,
            ev.module.c_str(),
            ev.pass,
            ev.n_seeds_in,
            ev.n_seeds,
            ev.n_tracks,
            ev.digest[0].size,
            ev.n_exported,
            ev.n_cand_distinct,
            ev.n_labeled,
            ev.n_out_matched,
            ev.n_stage[0],
            ev.n_stage[1],
            ev.n_stage[2],
            ev.n_stage_final[0],
            ev.n_stage_final[1],
            ev.n_stage_final[2],
            ev.n_on_final_path,
            counts_ok ? 1 : 0);
    }
  } catch (std::exception const& e) {
    std::fprintf(stderr, "mkFitTraceDigest: %s\n", e.what());
    return 1;
  }
  return 0;
}
