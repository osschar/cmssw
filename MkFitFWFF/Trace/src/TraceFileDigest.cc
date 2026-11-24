#include "MkFitFWFF/Trace/interface/TraceFileDigest.h"
#include "MkFitFWFF/Trace/interface/TraceRNTupleWriter.h"

#include "ROOT/RNTupleReader.hxx"

#include <cstdio>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>

namespace mkfit_fwff {

  namespace {

    using Key = std::tuple<std::uint32_t, std::uint32_t, std::uint64_t, std::string, int>;

    // Reads the RNTuple of record type R and digests it per event. extra(ev, rec) is called for every record.
    template <class R, class F>
    void readRecords(const std::string& fileName,
                     const char* name,
                     int vecIdx,
                     std::vector<TraceFileEvent>& evs,
                     const std::map<Key, std::size_t>& index,
                     F&& extra) {
      auto reader = ROOT::RNTupleReader::Open(name, fileName);
      auto vRun = reader->GetView<std::uint32_t>("run");
      auto vLumi = reader->GetView<std::uint32_t>("lumi");
      auto vEvent = reader->GetView<std::uint64_t>("event");
      auto vModule = reader->GetView<std::string>("module");
      auto vPass = reader->GetView<std::int32_t>("pass");
      auto vRec = reader->GetView<R>("rec");

      std::vector<trace_digest::Hasher> hashers(evs.size());
      std::vector<std::size_t> counts(evs.size(), 0);
      // The records of one event are contiguous: look the key up only when it changes.
      Key last;
      std::size_t idx = evs.size();
      for (auto i : reader->GetEntryRange()) {
        Key k{vRun(i), vLumi(i), vEvent(i), vModule(i), vPass(i)};
        if (idx == evs.size() || k != last) {
          auto it = index.find(k);
          if (it == index.end())
            throw std::runtime_error(std::string("readTraceFile: a record in ") + name + " has no entry in Events");
          idx = it->second;
          last = k;
        }
        const R& rec = vRec(i);
        trace_digest::digest(hashers[idx], rec);
        ++counts[idx];
        extra(evs[idx], rec);
      }
      for (std::size_t e = 0; e < evs.size(); ++e)
        evs[e].digest[vecIdx] = {trace_digest::kVectorNames[vecIdx], counts[e], hashers[e].m_h};
    }

    struct NoExtra {
      template <class R>
      void operator()(TraceFileEvent&, const R&) const {}
    };

  }  // namespace

  std::vector<TraceFileEvent> readTraceFile(const std::string& fileName) {
    std::vector<TraceFileEvent> evs;
    std::map<Key, std::size_t> index;
    {
      auto reader = ROOT::RNTupleReader::Open(TraceRNTupleWriter::kEventsName, fileName);
      auto vRun = reader->GetView<std::uint32_t>("run");
      auto vLumi = reader->GetView<std::uint32_t>("lumi");
      auto vEvent = reader->GetView<std::uint64_t>("event");
      auto vModule = reader->GetView<std::string>("module");
      auto vPass = reader->GetView<std::int32_t>("pass");
      auto vSeedsIn = reader->GetView<std::int32_t>("n_seeds_in");
      auto vSeeds = reader->GetView<std::int32_t>("n_seeds");
      auto vTracks = reader->GetView<std::int32_t>("n_tracks");
      const char* recNames[7] = {"n_TrCandMeta",
                                 "n_TrCandStage",
                                 "n_TrCandState",
                                 "n_TrLayerSearch",
                                 "n_TrHitMatch",
                                 "n_TrKalmanUpdate",
                                 "n_TrBkFitUpdate"};
      std::vector<ROOT::RNTupleView<std::int32_t>> vRec;
      for (int r = 0; r < 7; ++r)
        vRec.push_back(reader->GetView<std::int32_t>(recNames[r]));
      auto vOut = reader->GetView<std::vector<std::int32_t>>("out_labels");
      for (auto i : reader->GetEntryRange()) {
        TraceFileEvent ev;
        ev.run = vRun(i);
        ev.lumi = vLumi(i);
        ev.event = vEvent(i);
        ev.module = vModule(i);
        ev.pass = vPass(i);
        ev.n_seeds_in = vSeedsIn(i);
        ev.n_seeds = vSeeds(i);
        ev.n_tracks = vTracks(i);
        for (int r = 0; r < 7; ++r)
          ev.n_rec[r] = vRec[r](i);
        const auto& out = vOut(i);
        ev.out_labels.assign(out.begin(), out.end());
        if (!index.emplace(Key{ev.run, ev.lumi, ev.event, ev.module, ev.pass}, evs.size()).second)
          throw std::runtime_error("readTraceFile: a key appears twice in Events");
        evs.push_back(std::move(ev));
      }
    }

    // Exported metas' global_seed, per event, for matching the output tracks' labels.
    std::vector<std::multiset<int>> exportedSeeds(evs.size());
    std::vector<std::set<int>> cands(evs.size());
    readRecords<TrCandMeta>(fileName, "TrCandMeta", 0, evs, index, [&](TraceFileEvent& ev, const TrCandMeta& r) {
      const std::size_t e = &ev - evs.data();
      if (r.cand >= 0) {
        ++ev.n_exported;
        exportedSeeds[e].insert(r.global_seed);
        cands[e].insert(r.cand);
      }
      if (r.global_seed >= 0)
        ++ev.n_labeled;
    });
    readRecords<TrCandStage>(fileName, "TrCandStage", 1, evs, index, [](TraceFileEvent& ev, const TrCandStage& r) {
      if (r.stage >= 0 && r.stage < 3) {
        ++ev.n_stage[r.stage];
        if (r.final_state_id >= 0)
          ++ev.n_stage_final[r.stage];
      }
    });
    readRecords<TrCandState>(fileName, "TrCandState", 2, evs, index, [](TraceFileEvent& ev, const TrCandState& r) {
      if (r.on_final_path)
        ++ev.n_on_final_path;
    });
    readRecords<TrLayerSearch>(fileName, "TrLayerSearch", 3, evs, index, NoExtra());
    readRecords<TrHitMatch>(fileName, "TrHitMatch", 4, evs, index, NoExtra());
    readRecords<TrKalmanUpdate>(fileName, "TrKalmanUpdate", 5, evs, index, NoExtra());
    readRecords<TrBkFitUpdate>(fileName, "TrBkFitUpdate", 6, evs, index, NoExtra());

    for (std::size_t e = 0; e < evs.size(); ++e) {
      auto& ev = evs[e];
      ev.n_cand_distinct = cands[e].size();
      for (int l : ev.out_labels) {
        auto it = exportedSeeds[e].find(l);
        if (it != exportedSeeds[e].end()) {
          ++ev.n_out_matched;
          exportedSeeds[e].erase(it);
        }
      }
    }
    return evs;
  }

  std::string digestLine(std::uint32_t run,
                         std::uint32_t lumi,
                         std::uint64_t event,
                         std::string_view module,
                         int pass,
                         const trace_digest::VecDigest& d) {
    char buf[512];
    std::snprintf(buf,
                  sizeof(buf),
                  "MKFITTRACE digest %u:%u:%llu %.*s pass=%d %-17s %9zu %016llx",
                  run,
                  lumi,
                  (unsigned long long)event,
                  (int)module.size(),
                  module.data(),
                  pass,
                  d.name,
                  d.size,
                  (unsigned long long)d.hash);
    return buf;
  }

}  // namespace mkfit_fwff
