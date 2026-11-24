// Two ED modules for testing ModuleChanger (Schedule::changeModule): each instance counts the beginJob, endJob,
// beginStream and endStream calls it gets and prints one FWFFCHECK line per call. When the instance is destroyed it
// prints its counts and "OK" when every beginJob and beginStream had its matching end, "UNMATCHED" otherwise.
// Nothing throws: the lines are the result.
//
//   FWFFGlobalStreamCacheChecker: an edm::global::EDProducer with a StreamCache, like MkFitProducer. It wants stream
//     run transitions but no global ones.
//   FWFFStreamChecker: an edm::stream::EDProducer with a GlobalCache, which carries the counts of the instance (one
//     global cache per module instance; beginJob and endJob reach it as globalBeginJob and globalEndJob). It wants
//     neither stream nor global run transitions.
//
// The parameter 'value' does nothing; an edit can change it. With 'throwInConstructor' True the constructor throws,
// to test that a failed rebuild leaves the old module in place.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include <atomic>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

namespace fwfftest {
  inline std::atomic<int>& nextInstance() {
    static std::atomic<int> s_next{0};
    return s_next;
  }

  struct Counts {
    Counts(std::string kind, edm::ParameterSet const& pset)
        : kind_(std::move(kind)),
          label_(pset.getParameter<std::string>("@module_label")),
          instance_(nextInstance()++),
          value_(pset.getParameter<int>("value")) {
      print("construct");
    }
    ~Counts() {
      int const bj = beginJob_, ej = endJob_, bs = beginStream_, es = endStream_;
      std::ostringstream os;
      os << "beginJob=" << bj << " endJob=" << ej << " beginStream=" << bs << " endStream=" << es << " events=" << events_
         << " " << (bj == ej && bs == es ? "OK" : "UNMATCHED");
      print("destroy", os.str());
    }
    void print(char const* what, std::string const& extra = "") const {
      std::ostringstream os;
      os << "FWFFCHECK " << kind_ << " label=" << label_ << " instance=" << instance_ << " value=" << value_ << " "
         << what;
      if (!extra.empty())
        os << " " << extra;
      os << "\n";
      std::cout << os.str() << std::flush;
    }

    std::string const kind_;
    std::string const label_;
    int const instance_;
    int const value_;
    mutable std::atomic<int> beginJob_{0}, endJob_{0}, beginStream_{0}, endStream_{0}, events_{0};
  };

  void fillCheckerDescription(edm::ParameterSetDescription& desc) {
    desc.add<int>("value", 0)->setComment("does nothing; an edit can change it");
    desc.add<bool>("throwInConstructor", false)->setComment("True: the constructor throws");
  }

  void maybeThrow(edm::ParameterSet const& pset, Counts const& counts) {
    if (pset.getParameter<bool>("throwInConstructor")) {
      counts.print("throwInConstructor");
      throw cms::Exception("FWFFCheck") << "throwInConstructor is True for '" << counts.label_ << "'";
    }
  }
}  // namespace fwfftest

class FWFFGlobalStreamCacheChecker : public edm::global::EDProducer<edm::StreamCache<int>> {
public:
  explicit FWFFGlobalStreamCacheChecker(edm::ParameterSet const& pset) : counts_("global", pset) {
    fwfftest::maybeThrow(pset, counts_);
  }
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    fwfftest::fillCheckerDescription(desc);
    descriptions.add("fwffGlobalStreamCacheChecker", desc);
  }

  void beginJob() override {
    ++counts_.beginJob_;
    counts_.print("beginJob");
  }
  void endJob() override {
    ++counts_.endJob_;
    counts_.print("endJob");
  }
  std::unique_ptr<int> beginStream(edm::StreamID id) const override {
    ++counts_.beginStream_;
    counts_.print("beginStream", "stream=" + std::to_string(id.value()));
    return std::make_unique<int>(counts_.instance_);
  }
  void endStream(edm::StreamID id) const override {
    ++counts_.endStream_;
    counts_.print("endStream", "stream=" + std::to_string(id.value()));
  }
  void produce(edm::StreamID id, edm::Event&, edm::EventSetup const&) const override {
    ++counts_.events_;
    int const* cache = streamCache(id);
    if (cache == nullptr || *cache != counts_.instance_)
      counts_.print("produce", "BAD-STREAM-CACHE");
  }

private:
  fwfftest::Counts counts_;
};

class FWFFStreamChecker : public edm::stream::EDProducer<edm::GlobalCache<fwfftest::Counts>> {
public:
  FWFFStreamChecker(edm::ParameterSet const& pset, fwfftest::Counts const* counts) : stream_(-1) {
    fwfftest::maybeThrow(pset, *counts);
  }
  static std::unique_ptr<fwfftest::Counts> initializeGlobalCache(edm::ParameterSet const& pset) {
    return std::make_unique<fwfftest::Counts>("stream", pset);
  }
  static void globalBeginJob(fwfftest::Counts* counts) {
    ++counts->beginJob_;
    counts->print("beginJob");
  }
  static void globalEndJob(fwfftest::Counts* counts) {
    ++counts->endJob_;
    counts->print("endJob");
  }
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    fwfftest::fillCheckerDescription(desc);
    descriptions.add("fwffStreamChecker", desc);
  }

  void beginStream(edm::StreamID id) override {
    stream_ = static_cast<int>(id.value());
    ++globalCache()->beginStream_;
    globalCache()->print("beginStream", "stream=" + std::to_string(stream_));
  }
  void endStream() override {
    ++globalCache()->endStream_;
    globalCache()->print("endStream", "stream=" + std::to_string(stream_));
  }
  void produce(edm::Event&, edm::EventSetup const&) override { ++globalCache()->events_; }

private:
  int stream_;
};

DEFINE_FWK_MODULE(FWFFGlobalStreamCacheChecker);
DEFINE_FWK_MODULE(FWFFStreamChecker);
