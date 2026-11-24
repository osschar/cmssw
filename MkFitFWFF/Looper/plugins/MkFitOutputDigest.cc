// -*- C++ -*-
//
// Package:    MkFitFWFF/Looper
// Class:      MkFitOutputDigest
//
// Prints one line per event that identifies mkFit's output: the number of track candidates, their total number of
// hits and an order-independent 64-bit hash of them; optionally the same for a track collection. Two jobs that
// produce the same output print the same counts and hashes.
//
//   MKFITDIGEST tag=<tag> loop=<n> event=<run>:<lumi>:<event> cands=<n> candHits=<n> candHash=<hex>
//               tracks=<n> trackHits=<n> trackHash=<hex>
//
// A candidate hashes its seed reference index and its hits in their order, each hit as (DetId, hit type, cluster
// raw index); a track hashes its five parameters rounded to 1e-6, its charge and ndof, and its hits the same way
// when its TrackExtra is available. The collection hash is the sum of the per-element hashes, so it does not depend
// on the order of the collection. 'loop' is FWFFMkFitLooper's loop number, -1 without the looper. Fields of a
// collection that is not configured print '-'.

#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/trackerHitRTTI.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWFFLoopCounter.h"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {
  // splitmix64 finalizer
  uint64_t mix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
  }

  // Order-dependent hash of a sequence of 64-bit values.
  struct Hasher {
    uint64_t h = 0x6a09e667f3bcc909ULL;
    void add(uint64_t v) { h = mix64(h ^ v); }
  };

  uint64_t rounded(double x) {
    if (!std::isfinite(x))
      return 0xfff8000000000000ULL;
    return static_cast<uint64_t>(std::llround(x * 1e6));
  }

  void addHit(Hasher& h, TrackingRecHit const& hit) {
    h.add(hit.geographicalId().rawId());
    h.add(static_cast<uint64_t>(hit.getType()));
    if (!hit.isValid())
      return;
    // Hits from clusters of a real detector (not FastSim, not multi-hits): single, projected and matched.
    auto const rtti = trackerHitRTTI::rtti(hit);
    if (rtti >= trackerHitRTTI::single && rtti <= trackerHitRTTI::match) {
      h.add(static_cast<BaseTrackerRecHit const&>(hit).firstClusterRef().rawIndex());
      if (rtti == trackerHitRTTI::match)
        h.add(static_cast<SiStripMatchedRecHit2D const&>(hit).stereoClusterRef().rawIndex());
    }
  }

  std::string hex16(uint64_t v) {
    std::ostringstream os;
    os << std::hex << std::setw(16) << std::setfill('0') << v;
    return os.str();
  }
}  // namespace

class MkFitOutputDigest : public edm::global::EDAnalyzer<> {
public:
  explicit MkFitOutputDigest(edm::ParameterSet const& iConfig);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void analyze(edm::StreamID, edm::Event const& iEvent, edm::EventSetup const&) const override;

private:
  std::string const tag_;
  edm::EDGetTokenT<TrackCandidateCollection> candToken_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
};

MkFitOutputDigest::MkFitOutputDigest(edm::ParameterSet const& iConfig)
    : tag_{iConfig.getParameter<std::string>("tag")} {
  auto const cands = iConfig.getParameter<edm::InputTag>("candidates");
  if (!cands.label().empty())
    candToken_ = consumes(cands);
  auto const tracks = iConfig.getParameter<edm::InputTag>("tracks");
  if (!tracks.label().empty())
    trackToken_ = consumes(tracks);
}

void MkFitOutputDigest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("candidates", edm::InputTag("initialStepTrackCandidates"))
      ->setComment("TrackCandidateCollection of the mkFit output converter; empty label: not digested");
  desc.add<edm::InputTag>("tracks", edm::InputTag(""))
      ->setComment("reco::TrackCollection, e.g. initialStepTracks; empty label: not digested");
  desc.add<std::string>("tag", "")->setComment("printed on every line, to tell jobs apart");
  descriptions.addWithDefaultLabel(desc);
}

void MkFitOutputDigest::analyze(edm::StreamID, edm::Event const& iEvent, edm::EventSetup const&) const {
  std::ostringstream os;
  os << "MKFITDIGEST tag=" << tag_ << " loop=" << fwff::currentLoop().load() << " event=" << iEvent.id().run() << ":"
     << iEvent.id().luminosityBlock() << ":" << iEvent.id().event();

  if (!candToken_.isUninitialized()) {
    auto const& cands = iEvent.get(candToken_);
    uint64_t sum = 0;
    size_t nHits = 0;
    for (auto const& cand : cands) {
      Hasher h;
      h.add(cand.seedRef().isNonnull() ? static_cast<uint64_t>(cand.seedRef().key()) : ~0ULL);
      h.add(cand.nRecHits());
      for (auto const& hit : cand.recHits())
        addHit(h, hit);
      nHits += cand.nRecHits();
      sum += h.h;
    }
    os << " cands=" << cands.size() << " candHits=" << nHits << " candHash=" << hex16(sum);
  } else {
    os << " cands=- candHits=- candHash=-";
  }

  if (!trackToken_.isUninitialized()) {
    auto const& tracks = iEvent.get(trackToken_);
    uint64_t sum = 0;
    size_t nHits = 0;
    for (auto const& track : tracks) {
      Hasher h;
      for (unsigned int i = 0; i < reco::TrackBase::dimension; ++i)
        h.add(rounded(track.parameter(i)));
      h.add(static_cast<uint64_t>(track.charge()));
      h.add(rounded(track.ndof()));
      if (track.extra().isAvailable()) {
        h.add(track.recHitsSize());
        for (auto const& hit : track.recHits())
          addHit(h, *hit);
        nHits += track.recHitsSize();
      } else {
        h.add(track.hitPattern().numberOfAllHits(reco::HitPattern::TRACK_HITS));
        nHits += track.hitPattern().numberOfAllHits(reco::HitPattern::TRACK_HITS);
      }
      sum += h.h;
    }
    os << " tracks=" << tracks.size() << " trackHits=" << nHits << " trackHash=" << hex16(sum);
  } else {
    os << " tracks=- trackHits=- trackHash=-";
  }

  os << "\n";
  std::cout << os.str() << std::flush;
}

DEFINE_FWK_MODULE(MkFitOutputDigest);
