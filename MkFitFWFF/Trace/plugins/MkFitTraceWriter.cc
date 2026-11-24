// Package:    MkFitFWFF/Trace
// Class:      MkFitTraceWriter
//
// The service; its description is in MkFitTraceWriter.h.

#include "MkFitFWFF/Trace/plugins/MkFitTraceWriter.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "MkFitFWFF/Trace/interface/TraceFileDigest.h"
#include "RecoTracker/MkFitCore/interface/TraceStore.h"

#include <cstdio>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>

namespace mkfit_fwff {

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

  int MkFitTraceWriter::writeTruth(const TraceTruthTables& truth) {
    if (!writer_)
      throw cms::Exception("LogicError") << "MkFitTraceWriter: the trace file is not open";
    return writer_->writeTruth(truth);
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
