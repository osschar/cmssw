// -*- C++ -*-
//
// Package:    MkFitFWFF/Trace
// Class:      MkFitTraceWriter
//
// A service that writes the mkFit V2p2 trace of every MkFitProducer in the job into one ROOT file of plain RNTuples
// (TraceRNTupleWriter.h has the layout). It registers mkFit's trace sink (MkFitCore/interface/TraceStore.h), which
// each MkFitProducer calls at the end of its event with its stream's store. With writerThreads > 0 the sink takes the
// store's record vectors by swap and returns; the writer's threads write them (TraceRNTupleWriter.h). The producer
// clears its store at the start of each event and does not read it after the sink. The file is closed at the end of
// the job, after the queued events are written.
//
// Needs a build with -DMKFIT_TRACE in MkFitCore, MkFitCMS and RecoTracker/MkFit; without it the service refuses to
// start.
//
// With printDigest, it also prints for each event and record vector a digest line computed from the store in
// memory, in the format mkFitTraceDigest prints from the file:
//   MKFITTRACE digest <run>:<lumi>:<event> <module> pass=<p> <vector> <size> <hash>
// and one line per event, once it is written:
//   MKFITTRACE event <run>:<lumi>:<event> <module> pass=<p> records=<n> sink_ms=<s> write_ms=<w> done_ms=<d>
// sink_ms: the time the producer spent in the sink; write_ms: the time spent writing the event, summed over the
// RNTuples; done_ms: from the sink's call until the event was written.

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
    void eventDone(const mkfit::TraceStore& records, const TraceRNTupleWriter::EventDone& done);
    void finish();

    const std::string fileName_;
    const bool printDigest_;
    TraceRNTupleWriter::Options options_;
    std::unique_ptr<TraceRNTupleWriter> writer_;
    std::mutex printMutex_;
  };

  MkFitTraceWriter::MkFitTraceWriter(edm::ParameterSet const& iPS, [[maybe_unused]] edm::ActivityRegistry& iAR)
      : fileName_(iPS.getUntrackedParameter<std::string>("fileName")),
        printDigest_(iPS.getUntrackedParameter<bool>("printDigest")) {
    options_.compression = iPS.getUntrackedParameter<int>("compression");
    options_.maxUnzippedPageSize = iPS.getUntrackedParameter<unsigned int>("maxPageSize");
    options_.threads = iPS.getUntrackedParameter<int>("writerThreads");
#ifndef MKFIT_TRACE
    throw cms::Exception("Configuration")
        << "MkFitTraceWriter: this build has no mkFit trace. Build MkFitCore, MkFitCMS, RecoTracker/MkFit and "
           "MkFitFWFF/Trace with -DMKFIT_TRACE (e.g. scram b USER_CXXFLAGS=-DMKFIT_TRACE in an area of its own).";
#else
    try {
      writer_ = std::make_unique<TraceRNTupleWriter>(
          fileName_, options_, [this](const mkfit::TraceStore& r, const TraceRNTupleWriter::EventDone& d) {
            eventDone(r, d);
          });
    } catch (std::exception const& e) {
      throw cms::Exception("Configuration") << "MkFitTraceWriter: " << e.what();
    }
    if (mkfit::trace_sink())
      throw cms::Exception("Configuration") << "MkFitTraceWriter: another mkFit trace sink is registered";
    mkfit::set_trace_sink([this](const mkfit::TraceStore& s, const mkfit::TraceEventInfo& i) { writer_->write(s, i); });
    iAR.watchPostEndJob(this, &MkFitTraceWriter::finish);
    edm::LogInfo("MkFitTraceWriter") << "writing the mkFit trace to " << fileName_ << " (compression "
                                     << options_.compression << ", pages up to " << options_.maxUnzippedPageSize
                                     << " B, " << options_.threads << " writer threads)";
#endif
  }

  MkFitTraceWriter::~MkFitTraceWriter() {
    try {
      finish();
    } catch (std::exception const& e) {
      edm::LogError("MkFitTraceWriter") << "closing " << fileName_ << ": " << e.what();
    }
  }

  void MkFitTraceWriter::eventDone(const mkfit::TraceStore& records, const TraceRNTupleWriter::EventDone& d) {
    std::ostringstream os;
    if (printDigest_) {
      for (auto const& dg : trace_digest::digest_store(records))
        os << digestLine(d.run, d.lumi, d.event, d.module, d.pass, dg) << "\n";
    }
    char buf[320];
    std::snprintf(buf,
                  sizeof(buf),
                  "MKFITTRACE event %u:%u:%llu %.*s pass=%d records=%zu sink_ms=%.1f write_ms=%.1f done_ms=%.1f",
                  d.run,
                  d.lumi,
                  (unsigned long long)d.event,
                  (int)d.module.size(),
                  d.module.data(),
                  d.pass,
                  d.n_records,
                  d.take_ms,
                  d.write_ms,
                  d.done_ms);
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
    desc.addUntracked<int>("compression", 501)
        ->setComment("ROOT compression setting, 100 * algorithm + level; 501 is zstd at level 1, 505 ROOT's default");
    desc.addUntracked<unsigned int>("maxPageSize", 4 * 1024 * 1024)
        ->setComment("the size in bytes up to which an uncompressed RNTuple page grows; ROOT's default is 1 MiB");
    desc.addUntracked<int>("writerThreads", 4)
        ->setComment("threads that write the trace in the background; 0 writes in the producer's call");
    descriptions.add("MkFitTraceWriter", desc);
  }

}  // namespace mkfit_fwff

using mkfit_fwff::MkFitTraceWriter;
DEFINE_FWK_SERVICE(MkFitTraceWriter);
