// Package:    MkFitFWFF/Trace
// Class:      MkFitTraceTruth
//
// The truth of one MkFitProducer's traced event, in the way of the tracking ntuple (Validation/RecoTrack's
// TrackingNtuple) and of mkFit-external's writeMemoryFile, which makes the standalone's samples from it. It gives the
// tables to the MkFitTraceWriter service, which writes them into the trace file next to the records
// (TraceRNTupleWriter.h: TruthHits, SimTracks, TruthSeeds). The loader of the standalone joins them to the records.
//
// - The sim tracks are the TrackingParticles of trackingParticles, the ntuple's selection trackingParticlesIntime
//   (in time, charged and neutral, no other cut), in its order; a sim track's label is its index. As writeMemoryFile
//   makes them: the state at the production vertex, with the placeholder covariance (each variance the square of its
//   value) in mkFit's CCS parameters; the production type from the event id. chi2 holds 0: writeMemoryFile puts the
//   number of layers of the matched CMSSW reco track there, which this job does not have.
// - A mkFit hit is (layer, index): the index is the cluster's index in its cluster collection (MkFitHitWrapper), as
//   the records keep it. Its sim tracks are the TrackingParticles that ClusterTPAssociation (tpClusterProducer)
//   gives for its cluster, restricted to the sim tracks above. The hit is on the hit list of each of them. Its own
//   sim track (TruthHits sim, the standalone's MCHitInfo::mcTrackID) is the one with the largest share of the
//   cluster's charge in the digi-sim links (pixel: the links' fractions weighted by the pixels' ADC, as the ntuple's
//   chargeFraction(); outer tracker: the links' fractions), the first in the association's order on a tie.
//   writeMemoryFile chooses by the momentum of the sim hits instead, and drops a hit whose sim hit has less than
//   half the TrackingParticle's momentum or less than 50 MeV; that needs the sim hits of the pileup as well, which
//   the RAW file does not carry (the ntuple re-runs the mixing for them).
// - A sim track's hits are in the order of their distance from the production vertex (writeMemoryFile: the order of
//   the sim hits).
// - Each input seed: sim, sim_frac as the ntuple's see_simTrkIdx[0] and see_simTrkShareFrac[0] (the track
//   associator on the seeds as tracks); best_sim, best_frac as its see_bestSimTrkIdx and see_bestSimTrkShareFrac
//   (the sim track with most clusters of the seed, at least 3, the innermost on a tie).
//
// Prints one line per event:
//   MKFITTRUTH <run>:<lumi>:<event> <module> hits=<registered mkFit hits> hits_sim=<with a sim track>
//     hits_shared=<with more than one> sims=<sim tracks> sims_hits=<with a hit> seeds=<input seeds>
//     seeds_sim=<with sim> seeds_best=<with best_sim> pass=<the truth's pass> ms=<time>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/SiPixelDetId/interface/PixelChannelIdentifier.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackAssociation/interface/trackHitsToClusterRefs.h"
#include "SimTracker/TrackerHitAssociation/interface/ClusterTPAssociation.h"

#include "RecoTracker/MkFit/interface/MkFitClusterIndexToHit.h"
#include "RecoTracker/MkFit/interface/MkFitGeometry.h"
#include "RecoTracker/MkFit/interface/MkFitHitWrapper.h"
#include "RecoTracker/MkFit/interface/MkFitOutputWrapper.h"
#include "RecoTracker/MkFitCore/interface/Hit.h"
#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"

#include "MkFitFWFF/Trace/plugins/MkFitTraceWriter.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <limits>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace mkfit_fwff {

  class MkFitTraceTruth : public edm::global::EDAnalyzer<> {
  public:
    explicit MkFitTraceTruth(const edm::ParameterSet& iConfig);

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void analyze(edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

    using SimLinks = edm::DetSetVector<PixelDigiSimLink>;

    const std::string module_;
    const edm::EDGetTokenT<MkFitOutputWrapper> mkFitOutputToken_;
    const edm::EDGetTokenT<MkFitHitWrapper> pixelHitsToken_;
    const edm::EDGetTokenT<MkFitHitWrapper> stripHitsToken_;
    const edm::EDGetTokenT<MkFitClusterIndexToHit> pixelClusterIndexToken_;
    const edm::EDGetTokenT<MkFitClusterIndexToHit> stripClusterIndexToken_;
    const edm::EDGetTokenT<std::vector<int>> pixelLayerToken_;
    const edm::EDGetTokenT<std::vector<int>> stripLayerToken_;
    const edm::EDGetTokenT<TrajectorySeedCollection> seedsToken_;
    const edm::EDGetTokenT<edm::View<reco::Track>> seedTracksToken_;
    const edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> associatorToken_;
    const edm::EDGetTokenT<TrackingParticleRefVector> trackingParticlesToken_;
    const edm::EDGetTokenT<ClusterTPAssociation> clusterTPToken_;
    const edm::EDGetTokenT<SimLinks> pixelSimLinkToken_;
    const edm::EDGetTokenT<SimLinks> phase2OTSimLinkToken_;
    const edm::ESGetToken<MkFitGeometry, TrackerRecoGeometryRecord> mkFitGeomToken_;
  };

  namespace {
    // A SimTrack of one event: its id and the encoded event id.
    std::uint64_t simTrackKey(unsigned int trackId, EncodedEventId eid) {
      return (std::uint64_t(eid.rawId()) << 32) | trackId;
    }

    // The sum of the links' fractions on a channel of the det set, per sim track (index), weighted by w.
    void addLinks(const edm::DetSet<PixelDigiSimLink>& ds,
                  unsigned int channel,
                  double w,
                  const std::unordered_map<std::uint64_t, int>& simOfTrack,
                  std::unordered_map<int, double>& weight) {
      for (const PixelDigiSimLink& l : ds.data) {
        if (l.channel() != channel)
          continue;
        auto f = simOfTrack.find(simTrackKey(l.SimTrackId(), l.eventId()));
        if (f != simOfTrack.end())
          weight[f->second] += w * l.fraction();
      }
    }
  }  // namespace

  MkFitTraceTruth::MkFitTraceTruth(const edm::ParameterSet& iConfig)
      : module_(iConfig.getParameter<std::string>("mkFitProducer")),
        mkFitOutputToken_(consumes(edm::InputTag(module_))),
        pixelHitsToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelHits"))),
        stripHitsToken_(consumes(iConfig.getParameter<edm::InputTag>("stripHits"))),
        pixelClusterIndexToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelHits"))),
        stripClusterIndexToken_(consumes(iConfig.getParameter<edm::InputTag>("stripHits"))),
        pixelLayerToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelHits"))),
        stripLayerToken_(consumes(iConfig.getParameter<edm::InputTag>("stripHits"))),
        seedsToken_(consumes(iConfig.getParameter<edm::InputTag>("seeds"))),
        seedTracksToken_(consumes(iConfig.getParameter<edm::InputTag>("seedTracks"))),
        associatorToken_(consumes(iConfig.getParameter<edm::InputTag>("trackAssociator"))),
        trackingParticlesToken_(consumes(iConfig.getParameter<edm::InputTag>("trackingParticles"))),
        clusterTPToken_(consumes(iConfig.getParameter<edm::InputTag>("clusterTPMap"))),
        pixelSimLinkToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelSimLink"))),
        phase2OTSimLinkToken_(consumes(iConfig.getParameter<edm::InputTag>("phase2OTSimLink"))),
        mkFitGeomToken_(esConsumes()) {}

  void MkFitTraceTruth::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("mkFitProducer", "initialStepTrackCandidatesMkFit")
        ->setComment(
            "label of the traced MkFitProducer: the key of the truth, and its output orders this module after it");
    desc.add<edm::InputTag>("pixelHits", edm::InputTag("mkFitSiPixelHits"));
    desc.add<edm::InputTag>("stripHits", edm::InputTag("mkFitSiPhase2Hits"));
    desc.add<edm::InputTag>("seeds", edm::InputTag("initialStepSeeds"))->setComment("the producer's input seeds");
    desc.add<edm::InputTag>("seedTracks", edm::InputTag("seedTracksinitialStepSeeds"))
        ->setComment("the seeds as tracks (TrackFromSeedProducer), for the track associator");
    desc.add<edm::InputTag>("trackAssociator", edm::InputTag("quickTrackAssociatorByHits"));
    desc.add<edm::InputTag>("trackingParticles", edm::InputTag("trackingParticlesIntime"))
        ->setComment("TrackingParticleRefVector: the sim tracks, in order");
    desc.add<edm::InputTag>("clusterTPMap", edm::InputTag("tpClusterProducer"));
    desc.add<edm::InputTag>("pixelSimLink", edm::InputTag("simSiPixelDigis", "Pixel"));
    desc.add<edm::InputTag>("phase2OTSimLink", edm::InputTag("simSiPixelDigis", "Tracker"));
    descriptions.add("mkFitTraceTruth", desc);
  }

  void MkFitTraceTruth::analyze(edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
    const auto t0 = std::chrono::steady_clock::now();
    edm::Service<MkFitTraceWriter> writer;
    if (!writer.isAvailable())
      throw cms::Exception("Configuration") << "MkFitTraceTruth needs the MkFitTraceWriter service";

    // Only to run after the producer, whose trace the writer has then taken.
    iEvent.get(mkFitOutputToken_);

    const auto& trackerInfo = iSetup.getData(mkFitGeomToken_).trackerInfo();
    const TrackingParticleRefVector& tps = iEvent.get(trackingParticlesToken_);
    const ClusterTPAssociation& clusterTP = iEvent.get(clusterTPToken_);
    const SimLinks& pixelLinks = iEvent.get(pixelSimLinkToken_);
    const SimLinks& phase2Links = iEvent.get(phase2OTSimLinkToken_);

    std::unordered_map<unsigned int, int> simOfKey;     // TrackingParticle key -> sim track index
    std::unordered_map<std::uint64_t, int> simOfTrack;  // SimTrack -> sim track index
    for (int i = 0; i < (int)tps.size(); ++i) {
      simOfKey.emplace(tps[i].key(), i);
      const EncodedEventId eid = tps[i]->eventId();
      for (const SimTrack& g4 : tps[i]->g4Tracks())
        simOfTrack.emplace(simTrackKey(g4.trackId(), eid), i);
    }

    TraceTruthTables T;
    T.run = iEvent.id().run();
    T.lumi = iEvent.id().luminosityBlock();
    T.event = iEvent.id().event();
    T.module = module_;
    T.layers.resize(trackerInfo.n_layers());
    for (int l = 0; l < (int)T.layers.size(); ++l)
      T.layers[l].layer = l;

    // The hits of each sim track: (distance from its production vertex, layer, index).
    std::vector<std::vector<std::tuple<float, int, int>>> simHits(tps.size());
    long long nHits = 0, nHitsSim = 0, nHitsShared = 0;

    auto doHits = [&](const MkFitHitWrapper& wrapper,
                      const MkFitClusterIndexToHit& clusterIndex,
                      const std::vector<int>& layerOf,
                      const SimLinks& links,
                      bool pixel) {
      std::vector<int> cands;
      std::unordered_map<int, double> weight;
      for (int i = 0; i < (int)layerOf.size(); ++i) {
        const int ilay = layerOf[i];
        if (ilay < 0)
          continue;
        ++nHits;
        auto& L = T.layers.at(ilay);
        L.n_hits = std::max(L.n_hits, i + 1);
        const TrackingRecHit* rh = clusterIndex.hits()[i];
        const OmniClusterRef& cl = static_cast<const BaseTrackerRecHit*>(rh)->firstClusterRef();
        cands.clear();
        auto range = clusterTP.equal_range(cl);
        for (auto it = range.first; it != range.second; ++it) {
          auto f = simOfKey.find(it->second.key());
          if (f != simOfKey.end())
            cands.push_back(f->second);
        }
        if (cands.empty())
          continue;
        const mkfit::Hit& h = wrapper.hits()[i];
        for (int s : cands) {
          const TrackingParticle& tp = *tps[s];
          const float d = std::hypot(h.x() - tp.vx(), h.y() - tp.vy(), h.z() - tp.vz());
          simHits[s].emplace_back(d, ilay, i);
        }
        int best = cands.front();
        if (cands.size() > 1) {
          ++nHitsShared;
          weight.clear();
          auto ds = links.find(rh->geographicalId());
          if (ds != links.end()) {
            if (pixel) {
              const SiPixelCluster& c = cl.pixelCluster();
              for (int ip = 0; ip < c.size(); ++ip) {
                const SiPixelCluster::Pixel& px = c.pixel(ip);
                addLinks(*ds, PixelChannelIdentifier::pixelToChannel(px.x, px.y), px.adc, simOfTrack, weight);
              }
            } else {
              const Phase2TrackerCluster1D& c = cl.phase2OTCluster();
              for (unsigned int is = 0; is < c.size(); ++is)
                addLinks(*ds, Phase2TrackerDigi::pixelToChannel(c.firstRow() + is, c.column()), 1., simOfTrack, weight);
            }
          }
          double wBest = -1;
          for (int s : cands) {
            auto f = weight.find(s);
            const double w = f != weight.end() ? f->second : 0.;
            if (w > wBest) {
              wBest = w;
              best = s;
            }
          }
        }
        ++nHitsSim;
        L.hit.push_back(i);
        L.sim.push_back(best);
      }
    };
    doHits(iEvent.get(pixelHitsToken_),
           iEvent.get(pixelClusterIndexToken_),
           iEvent.get(pixelLayerToken_),
           pixelLinks,
           true);
    doHits(iEvent.get(stripHitsToken_),
           iEvent.get(stripClusterIndexToken_),
           iEvent.get(stripLayerToken_),
           phase2Links,
           false);

    // The sim tracks, as writeMemoryFile makes them.
    int nSimsHits = 0;
    T.simTracks.reserve(tps.size());
    for (int s = 0; s < (int)tps.size(); ++s) {
      const TrackingParticle& tp = *tps[s];
      const mkfit::SVector3 pos(tp.vx(), tp.vy(), tp.vz());
      const mkfit::SVector3 mom(tp.px(), tp.py(), tp.pz());
      mkfit::SMatrixSym66 err;
      err.At(0, 0) = tp.vx() * tp.vx();
      err.At(1, 1) = tp.vy() * tp.vy();
      err.At(2, 2) = tp.vz() * tp.vz();
      err.At(3, 3) = tp.px() * tp.px();
      err.At(4, 4) = tp.py() * tp.py();
      err.At(5, 5) = tp.pz() * tp.pz();
      mkfit::TrackState state(tp.charge(), pos, mom, err);
      state.convertFromCartesianToCCS();
      mkfit::Track& t = T.simTracks.emplace_back(state, 0.f, s, 0, nullptr);
      if (tp.eventId().bunchCrossing() == 0)
        t.setProdType(tp.eventId().event() == 0 ? mkfit::Track::ProdType::Signal : mkfit::Track::ProdType::InTimePU);
      else
        t.setProdType(mkfit::Track::ProdType::OutOfTimePU);
      auto& hits = simHits[s];
      std::sort(hits.begin(), hits.end());
      t.reserveHits(hits.size());
      for (const auto& [d, lay, idx] : hits)
        t.addHitIdx(idx, lay, 0.f);
      if (!hits.empty())
        ++nSimsHits;
    }

    // The seeds' sim tracks, as the ntuple has them.
    const TrajectorySeedCollection& seeds = iEvent.get(seedsToken_);
    T.seeds.resize(seeds.size());
    int nSeedsBest = 0, nSeedsSim = 0;
    for (int i = 0; i < (int)seeds.size(); ++i) {
      auto& ss = T.seeds[i];
      ss.seed = i;
      const auto& hits = seeds[i].recHits();
      const std::vector<OmniClusterRef> clusters = track_associator::hitsToClusterRefs(hits.begin(), hits.end());
      // TrackingNtuple's findBestMatchingTrackingParticle().
      std::unordered_map<int, std::pair<int, std::size_t>> count;  // sim -> (clusters, innermost cluster)
      for (std::size_t c = 0; c < clusters.size(); ++c) {
        auto range = clusterTP.equal_range(clusters[c]);
        for (auto it = range.first; it != range.second; ++it) {
          auto f = simOfKey.find(it->second.key());
          if (f == simOfKey.end())
            continue;
          auto& e = count.try_emplace(f->second, 0, std::numeric_limits<std::size_t>::max()).first->second;
          ++e.first;
          e.second = std::min(e.second, c);
        }
      }
      int bestCount = 2;
      std::size_t bestInner = std::numeric_limits<std::size_t>::max();
      for (const auto& [s, e] : count) {
        if (e.first > bestCount || (e.first == bestCount && e.second < bestInner)) {
          ss.best_sim = s;
          bestCount = e.first;
          bestInner = e.second;
        }
      }
      if (ss.best_sim >= 0) {
        ss.best_frac = clusters.empty() ? 0.f : float(bestCount) / clusters.size();
        ++nSeedsBest;
      }
    }
    {
      const auto& seedTracks = iEvent.get(seedTracksToken_);
      edm::RefToBaseVector<reco::Track> refs;
      for (std::size_t i = 0; i < seedTracks.size(); ++i)
        refs.push_back(seedTracks.refAt(i));
      const reco::RecoToSimCollection recSim = iEvent.get(associatorToken_).associateRecoToSim(refs, tps);
      for (std::size_t i = 0; i < refs.size(); ++i) {
        const auto key = refs[i]->seedRef().key();
        if (key >= T.seeds.size())
          throw cms::Exception("LogicError")
              << "MkFitTraceTruth: seed track " << i << " points at seed " << key << " of " << T.seeds.size();
        auto found = recSim.find(refs[i]);
        if (found == recSim.end() || found->val.empty())
          continue;
        auto f = simOfKey.find(found->val.front().first.key());
        if (f == simOfKey.end())
          continue;
        T.seeds[key].sim = f->second;
        T.seeds[key].sim_frac = found->val.front().second;
        ++nSeedsSim;
      }
    }

    const int pass = writer->writeTruth(T);

    char buf[512];
    std::snprintf(buf,
                  sizeof(buf),
                  "MKFITTRUTH %u:%u:%llu %s hits=%lld hits_sim=%lld hits_shared=%lld sims=%zu sims_hits=%d seeds=%zu "
                  "seeds_sim=%d seeds_best=%d pass=%d ms=%.1f",
                  T.run,
                  T.lumi,
                  T.event,
                  module_.c_str(),
                  nHits,
                  nHitsSim,
                  nHitsShared,
                  T.simTracks.size(),
                  nSimsHits,
                  T.seeds.size(),
                  nSeedsSim,
                  nSeedsBest,
                  pass,
                  std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
    std::cout << buf << std::endl;
  }

}  // namespace mkfit_fwff

using mkfit_fwff::MkFitTraceTruth;
DEFINE_FWK_MODULE(MkFitTraceTruth);
