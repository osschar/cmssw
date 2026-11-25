// Package:    MkFitFWFF/Trace
// Class:      MkFitTraceCrossingFrameToVector<T>
//
// Copies the objects of the MixingModule's crossing frames (CrossingFrame<T>, made with the crossing frames on) into
// plain vectors, which can be written to a file and read by a later job. MkFitTraceHitOrigin reads the pileup sim
// hits, SimTracks and SimVertices from them.
//
// For each InputTag in src, a std::vector<T> with the InputTag's instance as its instance: first the signal objects,
// then the pileup objects, in the crossing frame's order. That order is the order in which the MixingModule
// accumulated the pileup events (bunch crossings from minBunch to maxBunch, in each the events in the order of the
// playback record), and in each event the order of its own collection. The MixingModule has given every pileup object
// its event id (EncodedEventId: bunch crossing, and the event's number in it from 1) and shifted its time by the bunch
// crossing; a SimTrack keeps the index of its vertex in its own event's SimVertices. The signal objects are unchanged
// (event id 0).
//
// So for a tracker sim hit collection, an object's index in the vector is the index the digitizer counted for the
// digi-sim links (Phase2TrackerDigitizer::accumulate_local(): the signal's hits, then each pileup event's hits after
// all those before it). Instances of the module: MkFitTraceCrossingFramePSimHits, MkFitTraceCrossingFrameSimTracks,
// MkFitTraceCrossingFrameSimVertices.
//
// Prints one line per event and src:
//   MKFITCROSSINGFRAME <run>:<lumi>:<event> <module label> <src> signal=<n> pileup=<n> bcr_offsets=<the crossing
//     frame's offset of each bunch crossing's first pileup object>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace mkfit_fwff {

  template <typename T>
  class MkFitTraceCrossingFrameToVector : public edm::global::EDProducer<> {
  public:
    explicit MkFitTraceCrossingFrameToVector(const edm::ParameterSet& iConfig) {
      for (const edm::InputTag& t : iConfig.getParameter<std::vector<edm::InputTag>>("src")) {
        tags_.push_back(t);
        tokens_.push_back(consumes<CrossingFrame<T>>(t));
        produces<std::vector<T>>(t.instance());
      }
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<std::vector<edm::InputTag>>("src", {})->setComment(
          "crossing frames (CrossingFrame<T>); each gives a vector with the InputTag's instance");
      descriptions.addDefault(desc);
    }

    void produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const override {
      std::ostringstream os;
      for (std::size_t i = 0; i < tokens_.size(); ++i) {
        const CrossingFrame<T>& cf = iEvent.get(tokens_[i]);
        const std::vector<const T*>& sig = cf.getSignal();
        const std::vector<const T*>& pu = cf.getPileups();
        auto out = std::make_unique<std::vector<T>>();
        out->reserve(sig.size() + pu.size());
        for (const T* p : sig)
          out->push_back(*p);
        for (const T* p : pu)
          out->push_back(*p);
        os << "MKFITCROSSINGFRAME " << iEvent.id().run() << ":" << iEvent.id().luminosityBlock() << ":"
           << iEvent.id().event() << " " << moduleDescription().moduleLabel() << " " << tags_[i].encode()
           << " signal=" << sig.size() << " pileup=" << pu.size() << " bcr_offsets=";
        const std::vector<unsigned int>& off = cf.getPileupOffsetsBcr();
        for (std::size_t k = 0; k < off.size(); ++k)
          os << (k ? "," : "") << off[k];
        os << "\n";
        iEvent.put(std::move(out), tags_[i].instance());
      }
      std::cout << os.str() << std::flush;
    }

  private:
    std::vector<edm::InputTag> tags_;
    std::vector<edm::EDGetTokenT<CrossingFrame<T>>> tokens_;
  };

}  // namespace mkfit_fwff

using MkFitTraceCrossingFramePSimHits = mkfit_fwff::MkFitTraceCrossingFrameToVector<PSimHit>;
using MkFitTraceCrossingFrameSimTracks = mkfit_fwff::MkFitTraceCrossingFrameToVector<SimTrack>;
using MkFitTraceCrossingFrameSimVertices = mkfit_fwff::MkFitTraceCrossingFrameToVector<SimVertex>;
DEFINE_FWK_MODULE(MkFitTraceCrossingFramePSimHits);
DEFINE_FWK_MODULE(MkFitTraceCrossingFrameSimTracks);
DEFINE_FWK_MODULE(MkFitTraceCrossingFrameSimVertices);
