#include "MkFitFWFF/Trace/interface/TraceRNTupleWriter.h"

#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/interface/TraceStore.h"

#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleWriteOptions.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "TDirectory.h"
#include "TFile.h"

#include <cstdint>
#include <map>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace mkfit_fwff {

  namespace {

    // The key fields every RNTuple of the file starts with.
    struct KeyFields {
      std::shared_ptr<std::uint32_t> run, lumi;
      std::shared_ptr<std::uint64_t> event;
      std::shared_ptr<std::string> module;
      std::shared_ptr<std::int32_t> pass;

      void make(ROOT::RNTupleModel& m) {
        run = m.MakeField<std::uint32_t>("run");
        lumi = m.MakeField<std::uint32_t>("lumi");
        event = m.MakeField<std::uint64_t>("event");
        module = m.MakeField<std::string>("module", "label of the MkFitProducer");
        pass = m.MakeField<std::int32_t>("pass", "how often this (run, lumi, event, module) was written before");
      }
      void set(const mkfit::TraceEventInfo& info, int p) {
        *run = info.run;
        *lumi = info.lumi;
        *event = info.event;
        module->assign(info.module);
        *pass = p;
      }
    };

    ROOT::RNTupleWriteOptions writeOptions() {
      ROOT::RNTupleWriteOptions opt;
      // A hang was seen in cmsRun with implicit MT in ROOT I/O called from a framework thread
      // (mkfit-notes hgcal-showers, "Several events in flight"); one small file gains nothing from it.
      opt.SetUseImplicitMT(ROOT::RNTupleWriteOptions::EImplicitMT::kOff);
      return opt;
    }

    template <class R>
    struct RecordTuple {
      KeyFields key;
      std::shared_ptr<R> rec;
      std::unique_ptr<ROOT::RNTupleWriter> writer;

      void open(const char* name, TDirectory& dir) {
        auto model = ROOT::RNTupleModel::Create();
        key.make(*model);
        rec = model->template MakeField<R>("rec");
        writer = ROOT::RNTupleWriter::Append(std::move(model), name, dir, writeOptions());
      }
      void fill(const mkfit::TraceEventInfo& info, int p, const std::vector<R>& v) {
        key.set(info, p);
        for (const R& r : v) {
          *rec = r;
          writer->Fill();
        }
      }
    };

  }  // namespace

  struct TraceRNTupleWriter::Impl {
    std::unique_ptr<TFile> file;

    RecordTuple<TrCandMeta> metas;
    RecordTuple<TrCandStage> stages;
    RecordTuple<TrCandState> states;
    RecordTuple<TrLayerSearch> searches;
    RecordTuple<TrHitMatch> hitMatches;
    RecordTuple<TrKalmanUpdate> kalmanUpdates;
    RecordTuple<TrBkFitUpdate> bkFitUpdates;

    KeyFields evKey;
    std::shared_ptr<std::int32_t> nSeedsIn, nSeeds, nTracks;
    std::shared_ptr<std::int32_t> nRec[7];
    std::shared_ptr<std::vector<std::int32_t>> outLabels;
    std::unique_ptr<ROOT::RNTupleWriter> events;

    std::map<std::tuple<std::uint32_t, std::uint32_t, std::uint64_t, std::string>, int> passes;
  };

  TraceRNTupleWriter::TraceRNTupleWriter(const std::string& fileName) : m_impl(std::make_unique<Impl>()) {
    TDirectory::TContext context;  // TFile::Open makes the new file gDirectory; restore it on return
    m_impl->file.reset(TFile::Open(fileName.c_str(), "RECREATE"));
    if (!m_impl->file || m_impl->file->IsZombie())
      throw std::runtime_error("TraceRNTupleWriter: cannot create " + fileName);
    TDirectory& dir = *m_impl->file;

    m_impl->metas.open("TrCandMeta", dir);
    m_impl->stages.open("TrCandStage", dir);
    m_impl->states.open("TrCandState", dir);
    m_impl->searches.open("TrLayerSearch", dir);
    m_impl->hitMatches.open("TrHitMatch", dir);
    m_impl->kalmanUpdates.open("TrKalmanUpdate", dir);
    m_impl->bkFitUpdates.open("TrBkFitUpdate", dir);

    auto model = ROOT::RNTupleModel::Create();
    m_impl->evKey.make(*model);
    m_impl->nSeedsIn = model->MakeField<std::int32_t>("n_seeds_in", "seeds given to mkFit");
    m_impl->nSeeds = model->MakeField<std::int32_t>("n_seeds", "seeds after cleaning; TrCandMeta::seed indexes them");
    m_impl->nTracks = model->MakeField<std::int32_t>("n_tracks", "mkFit's output tracks, after the duplicate cleaner");
    const char* recNames[7] = {"n_TrCandMeta",
                               "n_TrCandStage",
                               "n_TrCandState",
                               "n_TrLayerSearch",
                               "n_TrHitMatch",
                               "n_TrKalmanUpdate",
                               "n_TrBkFitUpdate"};
    for (int i = 0; i < 7; ++i)
      m_impl->nRec[i] = model->MakeField<std::int32_t>(recNames[i]);
    m_impl->outLabels =
        model->MakeField<std::vector<std::int32_t>>("out_labels", "labels (input seed indices) of the output tracks");
    m_impl->events = ROOT::RNTupleWriter::Append(std::move(model), kEventsName, dir, writeOptions());
  }

  TraceRNTupleWriter::~TraceRNTupleWriter() { close(); }

  int TraceRNTupleWriter::write(const mkfit::TraceStore& s, const mkfit::TraceEventInfo& info) {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (!m_impl->file)
      throw std::runtime_error("TraceRNTupleWriter::write after close");
    TDirectory::TContext context;
    Impl& I = *m_impl;

    const int p = I.passes[{info.run, info.lumi, info.event, std::string(info.module)}]++;

    I.metas.fill(info, p, s.trCandMetas_);
    I.stages.fill(info, p, s.trCandStages_);
    I.states.fill(info, p, s.trCandStates_);
    I.searches.fill(info, p, s.trLayerSearches_);
    I.hitMatches.fill(info, p, s.trHitMatches_);
    I.kalmanUpdates.fill(info, p, s.trKalmanUpdates_);
    I.bkFitUpdates.fill(info, p, s.trBkFitUpdates_);

    I.evKey.set(info, p);
    *I.nSeedsIn = info.n_seeds_in;
    *I.nSeeds = info.seeds ? info.seeds->size() : 0;
    *I.nTracks = info.out_tracks ? info.out_tracks->size() : 0;
    const std::size_t sizes[7] = {s.trCandMetas_.size(),
                                  s.trCandStages_.size(),
                                  s.trCandStates_.size(),
                                  s.trLayerSearches_.size(),
                                  s.trHitMatches_.size(),
                                  s.trKalmanUpdates_.size(),
                                  s.trBkFitUpdates_.size()};
    for (int i = 0; i < 7; ++i)
      *I.nRec[i] = sizes[i];
    I.outLabels->clear();
    if (info.out_tracks) {
      for (const auto& t : *info.out_tracks)
        I.outLabels->push_back(t.label());
    }
    I.events->Fill();

    ++m_n_events;
    return p;
  }

  void TraceRNTupleWriter::close() {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (!m_impl || !m_impl->file)
      return;
    TDirectory::TContext context;
    Impl& I = *m_impl;
    // Destroying a writer commits its RNTuple into the file.
    I.metas.writer.reset();
    I.stages.writer.reset();
    I.states.writer.reset();
    I.searches.writer.reset();
    I.hitMatches.writer.reset();
    I.kalmanUpdates.writer.reset();
    I.bkFitUpdates.writer.reset();
    I.events.reset();
    I.file->Close();
    I.file.reset();
  }

}  // namespace mkfit_fwff
