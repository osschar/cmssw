#ifndef MkFitFWFF_Trace_plugins_MkFitTraceWriter_h
#define MkFitFWFF_Trace_plugins_MkFitTraceWriter_h

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
//
// Modules in this plugin library reach the service as edm::Service<mkfit_fwff::MkFitTraceWriter>: MkFitTraceTruth
// gives it the truth of each traced event (writeTruth()), MkFitTraceInputs mkFit's inputs of it (writeInputs():
// the hits, the beam spot, the dead bins, the geometry), and MkFitTraceHitOrigin the origin of each hit's charge
// (writeHitOrigin()), which it writes into the same file.

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"

#include "MkFitFWFF/Trace/interface/TraceRNTupleWriter.h"

#include <memory>
#include <mutex>
#include <string>

namespace mkfit_fwff {

  class MkFitTraceWriter {
  public:
    MkFitTraceWriter(edm::ParameterSet const& iPS, edm::ActivityRegistry& iAR);
    ~MkFitTraceWriter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    // Writes the truth of one traced event into the trace file (TraceRNTupleWriter::writeTruth()). Thread safe.
    // Throws when the writer is not open (a build without the trace).
    int writeTruth(const TraceTruthTables& truth);

    // Writes mkFit's inputs of one traced event into the trace file (TraceRNTupleWriter::writeInputs()). Thread
    // safe. Throws when the writer is not open (a build without the trace).
    int writeInputs(const TraceInputTables& inputs, bool* geomWritten = nullptr);

    // Writes the origin of each hit's charge of one traced event into the trace file
    // (TraceRNTupleWriter::writeHitOrigin()). Thread safe. Throws when the writer is not open.
    int writeHitOrigin(const TraceHitOriginTables& origin);

  private:
    void eventDone(const mkfit::TraceStore& records, const TraceRNTupleWriter::EventDone& done);
    void finish();

    const std::string fileName_;
    const bool printDigest_;
    TraceRNTupleWriter::Options options_;
    std::unique_ptr<TraceRNTupleWriter> writer_;
    std::mutex printMutex_;
  };

}  // namespace mkfit_fwff

#endif
