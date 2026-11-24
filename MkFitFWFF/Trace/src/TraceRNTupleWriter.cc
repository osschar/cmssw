#include "MkFitFWFF/Trace/interface/TraceRNTupleWriter.h"

#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/interface/TraceStore.h"

#include "ROOT/RNTupleFillStatus.hxx"
#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleWriteOptions.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RRawPtrWriteEntry.hxx"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <exception>
#include <map>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

namespace mkfit_fwff {

  namespace {

    using Clock = std::chrono::steady_clock;

    double msSince(Clock::time_point t0) {
      return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    }

    // One event as the RNTuples see it: the key, the Events fields, and where its records are.
    struct Batch {
      // The records, swapped out of the caller's store (threads > 0). Unused when writing in the caller.
      mkfit::TraceStore recs;
      // The records to write: &recs, or the caller's store.
      const mkfit::TraceStore* src = nullptr;

      std::uint32_t run = 0, lumi = 0;
      std::uint64_t event = 0;
      std::string module;
      std::int32_t pass = 0;
      std::int32_t nSeedsIn = 0, nSeeds = 0, nTracks = 0;
      std::int32_t nRec[7] = {};
      std::vector<std::int32_t> outLabels;

      Clock::time_point tCall;
      double takeMs = 0;
      mutable std::atomic<long long> writeNs{0};  // summed over the RNTuples, which fill a const batch
      std::atomic<int> pending{0};                // writer threads still to fill it

      void release() {
        // Keeps the capacity: the vectors go back into a store at a later write().
        recs.clear_trace();
        outLabels.clear();
        src = nullptr;
        writeNs = 0;
      }
    };

    ROOT::RNTupleWriteOptions writeOptions(const TraceRNTupleWriter::Options& o) {
      ROOT::RNTupleWriteOptions opt;
      // A hang was seen in cmsRun with implicit MT in ROOT I/O called from a framework thread
      // (mkfit-notes hgcal-showers, "Several events in flight"). The writer threads give the parallelism instead.
      opt.SetUseImplicitMT(ROOT::RNTupleWriteOptions::EImplicitMT::kOff);
      opt.SetCompression(o.compression);
      opt.SetMaxUnzippedPageSize(o.maxUnzippedPageSize);
      return opt;
    }

    // The key fields every RNTuple of the file starts with. The entry binds them to the batch's values.
    void makeKeyFields(ROOT::RNTupleModel& m) {
      m.MakeField<std::uint32_t>("run");
      m.MakeField<std::uint32_t>("lumi");
      m.MakeField<std::uint64_t>("event");
      m.MakeField<std::string>("module", "label of the MkFitProducer");
      m.MakeField<std::int32_t>("pass", "how often this (run, lumi, event, module) was written before");
    }

    void bindKey(ROOT::Detail::RRawPtrWriteEntry& e, const Batch& b) {
      e.BindRawPtr("run", &b.run);
      e.BindRawPtr("lumi", &b.lumi);
      e.BindRawPtr("event", &b.event);
      e.BindRawPtr("module", &b.module);
      e.BindRawPtr("pass", &b.pass);
    }

    // One RNTuple of the file, with one entry reused for every row.
    class Tuple {
    public:
      virtual ~Tuple() = default;

      // Fills the batch's rows. Writes a full cluster into the file under fileMutex.
      void fill(const Batch& b, std::mutex& fileMutex) {
        const auto t0 = Clock::now();
        fillRows(b, fileMutex);
        b.writeNs += std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - t0).count();
      }
      // Compresses the pages not yet compressed; no file access.
      void flushColumns() {
        if (m_writer)
          m_writer->FlushColumns();
      }
      // Writes the rest and commits the RNTuple into the file. Call under the file mutex.
      void commit() { m_writer.reset(); }

    protected:
      void open(std::unique_ptr<ROOT::RNTupleModel> model,
                const char* name,
                TDirectory& dir,
                const TraceRNTupleWriter::Options& o) {
        m_writer = ROOT::RNTupleWriter::Append(std::move(model), name, dir, writeOptions(o));
        m_entry = m_writer->CreateRawPtrWriteEntry();
      }
      void fillOne(std::mutex& fileMutex) {
        m_writer->FillNoFlush(*m_entry, m_status);
        if (m_status.ShouldFlushCluster()) {
          m_writer->FlushColumns();  // compression, outside the lock
          std::lock_guard<std::mutex> lock(fileMutex);
          TDirectory::TContext context;
          m_writer->FlushCluster();
        }
      }
      virtual void fillRows(const Batch& b, std::mutex& fileMutex) = 0;

      std::unique_ptr<ROOT::RNTupleWriter> m_writer;
      std::unique_ptr<ROOT::Detail::RRawPtrWriteEntry> m_entry;
      ROOT::RNTupleFillStatus m_status;
    };

    template <class R>
    class RecordTuple : public Tuple {
    public:
      using Member = std::vector<R> mkfit::TraceStore::*;

      RecordTuple(const char* name, Member member, TDirectory& dir, const TraceRNTupleWriter::Options& o)
          : m_member(member) {
        auto model = ROOT::RNTupleModel::Create();
        makeKeyFields(*model);
        model->template MakeField<R>("rec");
        open(std::move(model), name, dir, o);
        m_rec = m_entry->GetToken("rec");
      }

    private:
      void fillRows(const Batch& b, std::mutex& fileMutex) override {
        bindKey(*m_entry, b);
        for (const R& r : b.src->*m_member) {
          m_entry->BindRawPtr(m_rec, &r);
          fillOne(fileMutex);
        }
      }

      Member m_member;
      ROOT::RFieldToken m_rec;
    };

    class EventsTuple : public Tuple {
    public:
      static constexpr const char* kRecNames[7] = {"n_TrCandMeta",
                                                   "n_TrCandStage",
                                                   "n_TrCandState",
                                                   "n_TrLayerSearch",
                                                   "n_TrHitMatch",
                                                   "n_TrKalmanUpdate",
                                                   "n_TrBkFitUpdate"};

      EventsTuple(TDirectory& dir, const TraceRNTupleWriter::Options& o) {
        auto model = ROOT::RNTupleModel::Create();
        makeKeyFields(*model);
        model->MakeField<std::int32_t>("n_seeds_in", "seeds given to mkFit");
        model->MakeField<std::int32_t>("n_seeds", "seeds after cleaning; TrCandMeta::seed indexes them");
        model->MakeField<std::int32_t>("n_tracks", "mkFit's output tracks, after the duplicate cleaner");
        for (int i = 0; i < 7; ++i)
          model->MakeField<std::int32_t>(kRecNames[i]);
        model->MakeField<std::vector<std::int32_t>>("out_labels", "labels (input seed indices) of the output tracks");
        open(std::move(model), TraceRNTupleWriter::kEventsName, dir, o);
      }

    private:
      void fillRows(const Batch& b, std::mutex& fileMutex) override {
        bindKey(*m_entry, b);
        m_entry->BindRawPtr("n_seeds_in", &b.nSeedsIn);
        m_entry->BindRawPtr("n_seeds", &b.nSeeds);
        m_entry->BindRawPtr("n_tracks", &b.nTracks);
        for (int i = 0; i < 7; ++i)
          m_entry->BindRawPtr(kRecNames[i], &b.nRec[i]);
        m_entry->BindRawPtr("out_labels", &b.outLabels);
        fillOne(fileMutex);
      }
    };

  }  // namespace

  struct TraceRNTupleWriter::Impl {
    Options opt;
    EventDoneFn onDone;
    std::unique_ptr<TFile> file;

    // The RNTuples in the order of their size per event, which the writer threads take them in, round robin.
    std::vector<std::unique_ptr<Tuple>> tuples;
    // The same in the file's order: created and committed in it.
    std::vector<Tuple*> fileOrder;

    // mutex guards everything below; fileMutex the writing into the file.
    std::mutex mutex;
    std::mutex fileMutex;
    std::condition_variable cvWork;  // a batch queued, or stopping
    std::condition_variable cvFree;  // a batch released, or an error
    std::map<std::tuple<std::uint32_t, std::uint32_t, std::uint64_t, std::string>, int> passes;
    long long nEvents = 0;
    bool closed = false;
    bool stopping = false;
    std::exception_ptr error;

    // Writer threads (threads > 0).
    std::vector<std::thread> workers;
    std::vector<std::vector<Tuple*>> workerTuples;
    std::vector<std::deque<std::shared_ptr<Batch>>> queues;
    std::vector<std::shared_ptr<Batch>> freeBatches;
    int nBatches = 0;

    void setError(std::exception_ptr e) {
      std::lock_guard<std::mutex> lock(mutex);
      if (!error)
        error = e;
      cvFree.notify_all();
    }

    void eventDone(Batch& b) {
      if (!onDone)
        return;
      EventDone d;
      d.run = b.run;
      d.lumi = b.lumi;
      d.event = b.event;
      d.module = b.module;
      d.pass = b.pass;
      for (int i = 0; i < 7; ++i)
        d.n_records += b.nRec[i];
      d.take_ms = b.takeMs;
      d.write_ms = b.writeNs * 1e-6;
      d.done_ms = msSince(b.tCall);
      onDone(*b.src, d);
    }

    void workerLoop(int w) {
      for (;;) {
        std::shared_ptr<Batch> b;
        bool failed;
        {
          std::unique_lock<std::mutex> lock(mutex);
          cvWork.wait(lock, [&] { return !queues[w].empty() || stopping; });
          if (queues[w].empty())
            break;
          b = std::move(queues[w].front());
          queues[w].pop_front();
          failed = bool(error);
        }
        // After an error the batches are still taken and released, so that write() does not wait for ever.
        if (!failed) {
          try {
            for (Tuple* t : workerTuples[w])
              t->fill(*b, fileMutex);
          } catch (...) {
            setError(std::current_exception());
          }
        }
        // The last thread to finish a batch reports it. That thread fills the next batch only after the report, so
        // the reports come one at a time and in the order of the write() calls.
        if (--b->pending == 0) {
          try {
            eventDone(*b);
          } catch (...) {
            setError(std::current_exception());
          }
          b->release();
          std::lock_guard<std::mutex> lock(mutex);
          freeBatches.push_back(std::move(b));
          cvFree.notify_all();
        }
      }
      // The queue is empty and no more batches come: compress the last pages in parallel with the other threads.
      try {
        for (Tuple* t : workerTuples[w])
          t->flushColumns();
      } catch (...) {
        setError(std::current_exception());
      }
    }
  };

  TraceRNTupleWriter::TraceRNTupleWriter(const std::string& fileName, const Options& options, EventDoneFn onDone)
      : m_impl(std::make_unique<Impl>()) {
    Impl& I = *m_impl;
    I.opt = options;
    I.onDone = std::move(onDone);
    if (I.opt.threads > 0)
      ROOT::EnableThreadSafety();

    TDirectory::TContext context;  // TFile::Open makes the new file gDirectory; restore it on return
    I.file.reset(TFile::Open(fileName.c_str(), "RECREATE"));
    if (!I.file || I.file->IsZombie())
      throw std::runtime_error("TraceRNTupleWriter: cannot create " + fileName);
    TDirectory& dir = *I.file;

    // The RNTuples are created, and committed at close, in the order the file has always had them. The list for the
    // threads orders them by size.
    using TS = mkfit::TraceStore;
    auto metas = std::make_unique<RecordTuple<TrCandMeta>>("TrCandMeta", &TS::trCandMetas_, dir, I.opt);
    auto stages = std::make_unique<RecordTuple<TrCandStage>>("TrCandStage", &TS::trCandStages_, dir, I.opt);
    auto states = std::make_unique<RecordTuple<TrCandState>>("TrCandState", &TS::trCandStates_, dir, I.opt);
    auto searches = std::make_unique<RecordTuple<TrLayerSearch>>("TrLayerSearch", &TS::trLayerSearches_, dir, I.opt);
    auto hits = std::make_unique<RecordTuple<TrHitMatch>>("TrHitMatch", &TS::trHitMatches_, dir, I.opt);
    auto kalman = std::make_unique<RecordTuple<TrKalmanUpdate>>("TrKalmanUpdate", &TS::trKalmanUpdates_, dir, I.opt);
    auto bkfit = std::make_unique<RecordTuple<TrBkFitUpdate>>("TrBkFitUpdate", &TS::trBkFitUpdates_, dir, I.opt);
    auto events = std::make_unique<EventsTuple>(dir, I.opt);
    I.fileOrder = {
        metas.get(), stages.get(), states.get(), searches.get(), hits.get(), kalman.get(), bkfit.get(), events.get()};
    I.tuples.push_back(std::move(searches));
    I.tuples.push_back(std::move(states));
    I.tuples.push_back(std::move(hits));
    I.tuples.push_back(std::move(kalman));
    I.tuples.push_back(std::move(bkfit));
    I.tuples.push_back(std::move(stages));
    I.tuples.push_back(std::move(metas));
    I.tuples.push_back(std::move(events));

    if (I.opt.threads > 0) {
      const int n = std::min<int>(I.opt.threads, I.tuples.size());
      I.workerTuples.resize(n);
      I.queues.resize(n);
      for (std::size_t i = 0; i < I.tuples.size(); ++i)
        I.workerTuples[i % n].push_back(I.tuples[i].get());
      for (int w = 0; w < n; ++w)
        I.workers.emplace_back([this, w] { m_impl->workerLoop(w); });
    }
  }

  TraceRNTupleWriter::~TraceRNTupleWriter() {
    try {
      close();
    } catch (...) {
    }
  }

  long long TraceRNTupleWriter::nEventsWritten() const {
    std::lock_guard<std::mutex> lock(m_impl->mutex);
    return m_impl->nEvents;
  }

  int TraceRNTupleWriter::write(const mkfit::TraceStore& s, const mkfit::TraceEventInfo& info) {
    Impl& I = *m_impl;
    const auto t0 = Clock::now();
    std::unique_lock<std::mutex> lock(I.mutex);
    if (I.closed)
      throw std::runtime_error("TraceRNTupleWriter::write after close");
    if (I.error)
      std::rethrow_exception(I.error);

    std::shared_ptr<Batch> sb;
    Batch* b;
    Batch local;  // when writing here
    const bool async = !I.workers.empty();
    if (async) {
      I.cvFree.wait(
          lock, [&] { return !I.freeBatches.empty() || I.nBatches < I.opt.maxEventsInFlight || I.error || I.closed; });
      if (I.closed)
        throw std::runtime_error("TraceRNTupleWriter::write after close");
      if (I.error)
        std::rethrow_exception(I.error);
      if (!I.freeBatches.empty()) {
        sb = std::move(I.freeBatches.back());
        I.freeBatches.pop_back();
      } else {
        sb = std::make_shared<Batch>();
        ++I.nBatches;
      }
      b = sb.get();
    } else {
      b = &local;
    }

    b->tCall = t0;
    b->run = info.run;
    b->lumi = info.lumi;
    b->event = info.event;
    b->module.assign(info.module);
    b->pass = I.passes[{info.run, info.lumi, info.event, b->module}]++;
    b->nSeedsIn = info.n_seeds_in;
    b->nSeeds = info.seeds ? info.seeds->size() : 0;
    b->nTracks = info.out_tracks ? info.out_tracks->size() : 0;
    const std::size_t sizes[7] = {s.trCandMetas_.size(),
                                  s.trCandStages_.size(),
                                  s.trCandStates_.size(),
                                  s.trLayerSearches_.size(),
                                  s.trHitMatches_.size(),
                                  s.trKalmanUpdates_.size(),
                                  s.trBkFitUpdates_.size()};
    for (int i = 0; i < 7; ++i)
      b->nRec[i] = sizes[i];
    b->outLabels.clear();
    if (info.out_tracks) {
      for (const auto& t : *info.out_tracks)
        b->outLabels.push_back(t.label());
    }
    ++I.nEvents;
    const int pass = b->pass;

    if (async) {
      // Take the records; the store gets the batch's empty vectors, with their capacity.
      b->recs.trCandMetas_.swap(s.trCandMetas_);
      b->recs.trCandStages_.swap(s.trCandStages_);
      b->recs.trCandStates_.swap(s.trCandStates_);
      b->recs.trLayerSearches_.swap(s.trLayerSearches_);
      b->recs.trHitMatches_.swap(s.trHitMatches_);
      b->recs.trKalmanUpdates_.swap(s.trKalmanUpdates_);
      b->recs.trBkFitUpdates_.swap(s.trBkFitUpdates_);
      b->src = &b->recs;
      b->pending = I.queues.size();
      b->takeMs = msSince(t0);
      for (auto& q : I.queues)
        q.push_back(sb);
      lock.unlock();
      I.cvWork.notify_all();
      return pass;
    }

    // No writer threads: fill every RNTuple here, under the lock (write() is thread safe).
    b->src = &s;
    for (auto& t : I.tuples)
      t->fill(*b, I.fileMutex);
    b->takeMs = msSince(t0);
    I.eventDone(*b);
    return pass;
  }

  void TraceRNTupleWriter::close() {
    if (!m_impl)
      return;
    Impl& I = *m_impl;
    {
      std::lock_guard<std::mutex> lock(I.mutex);
      if (I.closed)
        return;
      I.closed = true;
      I.stopping = true;
    }
    I.cvWork.notify_all();
    I.cvFree.notify_all();
    for (auto& w : I.workers)
      w.join();
    I.workers.clear();

    {
      std::lock_guard<std::mutex> lock(I.fileMutex);
      TDirectory::TContext context;
      // Destroying a writer commits its RNTuple into the file.
      for (Tuple* t : I.fileOrder)
        t->commit();
      I.fileOrder.clear();
      I.tuples.clear();
      I.file->Close();
      I.file.reset();
    }
    if (I.error)
      std::rethrow_exception(I.error);
  }

}  // namespace mkfit_fwff
