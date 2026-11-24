// -*- C++ -*-
//
// Package:    MkFitFWFF/Trace
// Class:      MkFitTraceWriter
//
// A service that writes the mkFit V2p2 trace of every MkFitProducer in the job into one ROOT file of plain RNTuples
// (TraceRNTupleWriter.h has the layout). It registers mkFit's trace sink (MkFitCore/interface/TraceStore.h), which
// each MkFitProducer calls at the end of its event with its stream's store; the sink appends the store under a mutex.
// The file is closed at the end of the job.
//
// Needs a build with -DMKFIT_TRACE in MkFitCore, MkFitCMS and RecoTracker/MkFit; without it the service refuses to
// start.
//
// With printDigest, it also prints for each event and record vector a digest line computed from the store in
// memory, in the format mkFitTraceDigest prints from the file:
//   MKFITTRACE digest <run>:<lumi>:<event> <module> pass=<p> <vector> <size> <hash>
// and one line per event with the writing time:
//   MKFITTRACE event <run>:<lumi>:<event> <module> pass=<p> records=<n> write_ms=<t>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "MkFitFWFF/Trace/interface/TraceFileDigest.h"
#include "MkFitFWFF/Trace/interface/TraceRNTupleWriter.h"
#include "RecoTracker/MkFitCore/interface/TraceStore.h"

#include <chrono>
#include <cstdio>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>

namespace mkfit_fwff {

  class MkFitTraceWriter {
  public:
    MkFitTraceWriter(edm::ParameterSet const& iPS, edm::ActivityRegistry& iAR);
    ~MkFitTraceWriter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void write(const mkfit::TraceStore& store, const mkfit::TraceEventInfo& info);
    void finish();

    const std::string fileName_;
    const bool printDigest_;
    std::unique_ptr<TraceRNTupleWriter> writer_;
    std::mutex printMutex_;
  };

  MkFitTraceWriter::MkFitTraceWriter(edm::ParameterSet const& iPS, [[maybe_unused]] edm::ActivityRegistry& iAR)
      : fileName_(iPS.getUntrackedParameter<std::string>("fileName")),
        printDigest_(iPS.getUntrackedParameter<bool>("printDigest")) {
#ifndef MKFIT_TRACE
    throw cms::Exception("Configuration")
        << "MkFitTraceWriter: this build has no mkFit trace. Build MkFitCore, MkFitCMS, RecoTracker/MkFit and "
           "MkFitFWFF/Trace with -DMKFIT_TRACE (e.g. scram b USER_CXXFLAGS=-DMKFIT_TRACE in an area of its own).";
#else
    try {
      writer_ = std::make_unique<TraceRNTupleWriter>(fileName_);
    } catch (std::exception const& e) {
      throw cms::Exception("Configuration") << "MkFitTraceWriter: " << e.what();
    }
    if (mkfit::trace_sink())
      throw cms::Exception("Configuration") << "MkFitTraceWriter: another mkFit trace sink is registered";
    mkfit::set_trace_sink([this](const mkfit::TraceStore& s, const mkfit::TraceEventInfo& i) { write(s, i); });
    iAR.watchPostEndJob(this, &MkFitTraceWriter::finish);
    edm::LogInfo("MkFitTraceWriter") << "writing the mkFit trace to " << fileName_;
#endif
  }

  MkFitTraceWriter::~MkFitTraceWriter() { finish(); }

  void MkFitTraceWriter::write(const mkfit::TraceStore& store, const mkfit::TraceEventInfo& info) {
    const auto t0 = std::chrono::steady_clock::now();
    const int pass = writer_->write(store, info);
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();

    const std::size_t nRec = store.trCandMetas_.size() + store.trCandStages_.size() + store.trCandStates_.size() +
                             store.trLayerSearches_.size() + store.trHitMatches_.size() +
                             store.trKalmanUpdates_.size() + store.trBkFitUpdates_.size();
    std::ostringstream os;
    if (printDigest_) {
      for (auto const& d : trace_digest::digest_store(store))
        os << digestLine(info.run, info.lumi, info.event, info.module, pass, d) << "\n";
    }
    char buf[256];
    std::snprintf(buf,
                  sizeof(buf),
                  "MKFITTRACE event %u:%u:%llu %.*s pass=%d records=%zu write_ms=%.1f",
                  info.run,
                  info.lumi,
                  (unsigned long long)info.event,
                  (int)info.module.size(),
                  info.module.data(),
                  pass,
                  nRec,
                  ms);
    os << buf << "\n";
    std::lock_guard<std::mutex> lock(printMutex_);
    std::cout << os.str() << std::flush;
  }

  void MkFitTraceWriter::finish() {
    if (!writer_)
      return;
#ifdef MKFIT_TRACE
    mkfit::set_trace_sink(nullptr);
#endif
    writer_->close();
    edm::LogInfo("MkFitTraceWriter") << "wrote the mkFit trace of " << writer_->nEventsWritten() << " events to "
                                     << fileName_;
    writer_.reset();
  }

  void MkFitTraceWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.addUntracked<std::string>("fileName", "mkfit-trace.root")->setComment("the ROOT file of the trace RNTuples");
    desc.addUntracked<bool>("printDigest", false)
        ->setComment("print a digest line per event and record vector, computed from the store in memory");
    descriptions.add("MkFitTraceWriter", desc);
  }

}  // namespace mkfit_fwff

using mkfit_fwff::MkFitTraceWriter;
DEFINE_FWK_SERVICE(MkFitTraceWriter);
