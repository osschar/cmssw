// Package:    MkFitFWFF/Trace
// Class:      MkFitTraceTruthNtuple
//
// The truth source "ntuple" of one MkFitProducer's traced event: what the tracking ntuple
// (Validation/RecoTrack/plugins/TrackingNtuple.cc) and mkFit-external's writeMemoryFile
// (RecoTracker/MkFitCMS/standalone/tkNtuple/WriteMemoryFile.cc) give a standalone sample, made in memory from the
// products of the job. It gives the tables to the MkFitTraceWriter service as the truth source given by "source"
// (TraceRNTupleWriter.h: TruthHits_<source>, SimTracks_<source>, TruthSeeds_<source>, SimHitStates_<source>).
//
// The tracking ntuple, as customiseTrackingNtuple and extendedContent configure it (CMSSW_20_1_0_pre2; the code of
// TrackingNtuple.cc and of the association modules is the same in the IB):
// - TrackingParticles: trackingParticles, the ntuple's selection (trackingParticleRefSelector with signalOnly,
//   chargedOnly False, tip and lip 1e5, rapidity within 10, ptMin 0) with intimeOnly False (extendedContent): every
//   bunch crossing of the mixing's TrackingParticles. Their index is the sim track's index and label.
// - Sim hits (fillSimHits): the entries of simHitTPMap (SimHitTPAssociationProducer on the sim hits of the whole
//   crossing, signal first: the playback of the mixing), in its order, whose TrackingParticle is one of the above and
//   whose DetId is in the tracker. A sim hit's index is its position among them. Its global position and momentum at
//   entry (TrackerGeometry, as floats), particle type, process, TrackingParticle and time of flight are kept.
// - Per rec hit (matchCluster(), pixel clusters then outer-tracker clusters, each in cluster index order): for each
//   TrackingParticle that clusterTPMap (ClusterTPAssociation, by track id and event id) gives the cluster, in its order
//   and among the selected ones, its sim hits on the hit's module, in simHitTPMap's order; with keepEleSimHits (as
//   extendedContent sets it) the electron sim hits of a non-electron TrackingParticle only when it has no other sim
//   hit on the module. The rec hit goes on each such sim hit's list of rec hits. A TrackingParticle of the cluster
//   without a sim hit on the module stops the ntuple with an exception; here it is counted (nosimhit).
// - A TrackingParticle's sim hits (sim_simHitIdx): ordered by time of flight (std::sort, as the ntuple sorts them).
// - The first track of a TrackingParticle (sim_trkIdx[0]): the first entry of trackAssociator's associateSimToReco()
//   over tracks (the ntuple's tracks: generalTracks) and the selected TrackingParticles.
//
// writeMemoryFile (no --clean-sim-tracks, as the sample was written), every selection it applies:
// - bestTkIdx(), per rec hit over its sim hits (the list above): a sim hit counts only with a TrackingParticle, with
//   |p| >= 0.05 GeV (p at entry) and |p| >= 0.5 |p| of its TrackingParticle; of those the one with the largest ratio
//   wins (the first on a tie). Split-cluster arbitration: when the rec hit has exactly one sim hit and it won, the
//   label stays only if this rec hit is the first rec hit of the same type on that sim hit's list. The label is the
//   winner's TrackingParticle (TruthHits_<source>); the winning sim hit stays even when the arbitration clears the
//   label (SimHitStates_<source>: Event::simHitStates_, and Event::simTrackStates_ as the state at the sim hit).
//   Computed in floats as writeMemoryFile computes them from the ntuple's float branches, without fused multiply-adds
//   (normNoFma()).
// - Sim tracks: every TrackingParticle, in order; the state at its parent vertex with the placeholder covariance
//   (each variance the square of its value) in CCS; the production type from the event id; chi2 = the number of
//   mkFit layers with a valid pixel or outer-tracker hit on its first track (0 without one); hits: for each sim hit in
//   time order, its rec hits in their order, as (mkFit layer, cluster index). No hit is filtered here; the loader
//   filters the hits of other labels (Event::filterOutMislabeledHitsInSimTracks()), as the standalone does.
// - Seed labels (TruthSeeds_<source> sim): for each TrackingParticle in order whose first track's seed is one of
//   seeds (the producer's input seeds), that seed gets the TrackingParticle's index; a later one overwrites it.
//   writeMemoryFile keeps the initialStep and hltIter0 seeds with at least 3 pixel hits; every input seed is written
//   here, with its hits (hit_lyr, hit_idx) so that seeds can be matched.
// - Not applicable: the strip and glued branches, the cluster charge cut (phase-2 inputs ignore it), the rec tracks
//   (--write-rec-tracks), the hit masks, the beam spot (it is in Inputs; the ntuple writes bsp_z from x0()).
//
// The digi-sim links' sim hit index (PixelDigiSimLink::CFposition()) is read by neither the ntuple nor its association
// modules: they go by track id and event id. Here it serves a check: for each hit's chosen sim hit, whether a link on
// the cluster's channels names it, under the decoding the event's links agree with (linkIndex auto: global, or
// per module as the digitizer of CMSSW_20_0_0_pre1 counted; see MkFitTraceHitOrigin). simHits: the sim hit vectors
// simHitTPMap was made from (the whole crossing, signal first); empty: no check (linked -1).
//
// Prints one line per event:
//   MKFITTRUTHNT <run>:<lumi>:<event> <module> source=<s> tps=<n> tps_bx=<bx:n,...> simhits=<n> hits=<registered>
//     hits_simhit=<with a matched sim hit> labels=<n> lost_p50=<all sim hits below 50 MeV> lost_half=<none at half the
//     particle's p> lost_arb=<cleared by the arbitration> states=<n> nosimhit=<n> sims_hits=<sim tracks with a hit>
//     sims_trk=<with a first track> seeds=<n> seeds_lab=<labelled> link_index=<..> linked=<n>/<checked>
//     linked_other=<n under the other decoding> events_links=<all links> found_global=<n> found_perModule=<n> pass=<p>
//     ms=<time>
// With dumpPrefix, <dumpPrefix>hits.txt (per rec hit with matched sim hits: type 0 pixel or 4 outer, cluster index,
// layer, its sim hit indices) and <dumpPrefix>simhits.txt (per sim hit on a rec hit: index, DetId, x, y, z, its
// TrackingParticle), to compare with the ntuple's pix_simHitIdx, ph2_simHitIdx and simhit_* branches.

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PixelChannelIdentifier.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimGeneral/TrackingAnalysis/interface/SimHitTPAssociationProducer.h"
#include "SimTracker/TrackerHitAssociation/interface/ClusterTPAssociation.h"

#include "RecoTracker/MkFit/interface/MkFitClusterIndexToHit.h"
#include "RecoTracker/MkFit/interface/MkFitGeometry.h"
#include "RecoTracker/MkFit/interface/MkFitOutputWrapper.h"
#include "RecoTracker/MkFitCore/interface/Track.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"

#include "MkFitFWFF/Trace/plugins/MkFitTraceWriter.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mkfit_fwff {

  class MkFitTraceTruthNtuple : public edm::global::EDAnalyzer<> {
  public:
    explicit MkFitTraceTruthNtuple(const edm::ParameterSet& iConfig);

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void analyze(edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

    using SimLinks = edm::DetSetVector<PixelDigiSimLink>;
    using SimHitTPAssociationList = SimHitTPAssociationProducer::SimHitTPAssociationList;
    enum class LinkIndex { kAuto, kGlobal, kPerModule };

    const std::string module_;
    const std::string source_;
    const edm::EDGetTokenT<MkFitOutputWrapper> mkFitOutputToken_;
    const edm::EDGetTokenT<MkFitClusterIndexToHit> pixelClusterIndexToken_;
    const edm::EDGetTokenT<MkFitClusterIndexToHit> stripClusterIndexToken_;
    const edm::EDGetTokenT<std::vector<int>> pixelLayerToken_;
    const edm::EDGetTokenT<std::vector<int>> stripLayerToken_;
    const edm::EDGetTokenT<TrajectorySeedCollection> seedsToken_;
    const edm::EDGetTokenT<edm::View<reco::Track>> tracksToken_;
    const edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> associatorToken_;
    const edm::EDGetTokenT<TrackingParticleRefVector> trackingParticlesToken_;
    const edm::EDGetTokenT<SimHitTPAssociationList> simHitTPMapToken_;
    const edm::EDGetTokenT<ClusterTPAssociation> clusterTPToken_;
    const edm::EDGetTokenT<SimLinks> pixelSimLinkToken_;
    const edm::EDGetTokenT<SimLinks> phase2OTSimLinkToken_;
    std::vector<edm::EDGetTokenT<std::vector<PSimHit>>> simHitTokens_;
    std::vector<unsigned int> simHitTof_;  // PixelDigiSimLink::LowTof or HighTof, from the instance name
    const bool keepEleSimHits_;
    LinkIndex linkIndex_;
    const std::string dumpPrefix_;
    const edm::ESGetToken<MkFitGeometry, TrackerRecoGeometryRecord> mkFitGeomToken_;
    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeomToken_;
  };

  namespace {
    std::uint64_t detKey(unsigned int detId, unsigned int x) { return (std::uint64_t(detId) << 32) | x; }

    // The ntuple's sim hit (fillSimHits()), with what writeMemoryFile reads of it.
    struct NtSimHit {
      float x, y, z, px, py, pz, tof;
      int particle, process, tp;
      unsigned int detId;
      TrackPSimHitRef ref;
      std::vector<int> hitIdx, hitType;  // its rec hits: cluster index and HitType (0 pixel, 4 outer tracker)
    };

    // TrackingNtuple's TPHitIndex and its sort.
    struct TPHitIndex {
      unsigned int tpKey, simHitIdx;
      float tof;
      unsigned int detId;
    };
    bool tpHitIndexListLess(const TPHitIndex& i, const TPHitIndex& j) { return i.tpKey < j.tpKey; }
    bool tpHitIndexListLessSort(const TPHitIndex& i, const TPHitIndex& j) {
      if (i.tpKey == j.tpKey) {
        if (edm::isNotFinite(i.tof) && edm::isNotFinite(j.tof))
          return i.detId < j.detId;
        return i.tof < j.tof;  // works as intended if either one is NaN
      }
      return i.tpKey < j.tpKey;
    }

    constexpr int kHitPixel = 0, kHitPhase2OT = 4;  // TrackingNtuple's HitType

    // |p| from float components as writeMemoryFile computes it, sqrt((x * x + y * y) + z * z) with each operation
    // rounded on its own: the volatile products keep the compiler from fusing them into the additions (an FMA). The
    // sample's writeMemoryFile was built without FMA; with fused products the ratio |p_hit| / |p_particle| of a sim hit
    // at its particle's production point (both near 1) can come out one float ULP apart, which decides between the
    // two members of a pair made at that point.
    float normNoFma(float x, float y, float z) {
      const volatile float xx = x * x, yy = y * y, zz = z * z;
      const volatile float s = xx + yy;
      return std::sqrt(s + zz);
    }
  }  // namespace

  MkFitTraceTruthNtuple::MkFitTraceTruthNtuple(const edm::ParameterSet& iConfig)
      : module_(iConfig.getParameter<std::string>("mkFitProducer")),
        source_(iConfig.getParameter<std::string>("source")),
        mkFitOutputToken_(consumes(edm::InputTag(module_))),
        pixelClusterIndexToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelHits"))),
        stripClusterIndexToken_(consumes(iConfig.getParameter<edm::InputTag>("stripHits"))),
        pixelLayerToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelHits"))),
        stripLayerToken_(consumes(iConfig.getParameter<edm::InputTag>("stripHits"))),
        seedsToken_(consumes(iConfig.getParameter<edm::InputTag>("seeds"))),
        tracksToken_(consumes(iConfig.getParameter<edm::InputTag>("tracks"))),
        associatorToken_(consumes(iConfig.getParameter<edm::InputTag>("trackAssociator"))),
        trackingParticlesToken_(consumes(iConfig.getParameter<edm::InputTag>("trackingParticles"))),
        simHitTPMapToken_(consumes(iConfig.getParameter<edm::InputTag>("simHitTPMap"))),
        clusterTPToken_(consumes(iConfig.getParameter<edm::InputTag>("clusterTPMap"))),
        pixelSimLinkToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelSimLink"))),
        phase2OTSimLinkToken_(consumes(iConfig.getParameter<edm::InputTag>("phase2OTSimLink"))),
        keepEleSimHits_(iConfig.getParameter<bool>("keepEleSimHits")),
        dumpPrefix_(iConfig.getUntrackedParameter<std::string>("dumpPrefix")),
        mkFitGeomToken_(esConsumes()),
        trackerGeomToken_(esConsumes()) {
    if (source_.empty() || source_ == "tp")
      throw cms::Exception("Configuration") << "MkFitTraceTruthNtuple: source must be a name other than 'tp'";
    for (const edm::InputTag& t : iConfig.getParameter<std::vector<edm::InputTag>>("simHits")) {
      simHitTokens_.push_back(consumes<std::vector<PSimHit>>(t));
      simHitTof_.push_back(t.instance().find("HighTof") != std::string::npos ? PixelDigiSimLink::HighTof
                                                                             : PixelDigiSimLink::LowTof);
    }
    const std::string li = iConfig.getParameter<std::string>("linkIndex");
    if (li == "auto")
      linkIndex_ = LinkIndex::kAuto;
    else if (li == "global")
      linkIndex_ = LinkIndex::kGlobal;
    else if (li == "perModule")
      linkIndex_ = LinkIndex::kPerModule;
    else
      throw cms::Exception("Configuration")
          << "MkFitTraceTruthNtuple: linkIndex must be auto, global or perModule, not '" << li << "'";
  }

  void MkFitTraceTruthNtuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("mkFitProducer", "initialStepTrackCandidatesMkFit")
        ->setComment(
            "label of the traced MkFitProducer: the key of the truth, and its output orders this module after it");
    desc.add<std::string>("source", "ntuple")->setComment("the truth source's name in the trace file");
    desc.add<edm::InputTag>("pixelHits", edm::InputTag("mkFitSiPixelHits"));
    desc.add<edm::InputTag>("stripHits", edm::InputTag("mkFitSiPhase2Hits"));
    desc.add<edm::InputTag>("seeds", edm::InputTag("initialStepSeeds"))->setComment("the producer's input seeds");
    desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"))->setComment("the ntuple's tracks");
    desc.add<edm::InputTag>("trackAssociator", edm::InputTag("quickTrackAssociatorByHits"));
    desc.add<edm::InputTag>("trackingParticles", edm::InputTag("trackingParticlesNtuple"))
        ->setComment("TrackingParticleRefVector, the ntuple's selection: the sim tracks, in order");
    desc.add<edm::InputTag>("simHitTPMap", edm::InputTag("simHitTPAssocNtuple"))
        ->setComment("SimHitTPAssociationProducer on the sim hits of the whole crossing, as the ntuple's");
    desc.add<edm::InputTag>("clusterTPMap", edm::InputTag("tpClusterProducer"));
    desc.add<edm::InputTag>("pixelSimLink", edm::InputTag("simSiPixelDigis", "Pixel"));
    desc.add<edm::InputTag>("phase2OTSimLink", edm::InputTag("simSiPixelDigis", "Tracker"));
    desc.add<std::vector<edm::InputTag>>("simHits", {})
        ->setComment(
            "the sim hit vectors of the whole crossing that simHitTPMap was made from, for the link check; empty: no "
            "check");
    desc.add<bool>("keepEleSimHits", true)->setComment("the ntuple's keepEleSimHits (extendedContent sets it)");
    desc.add<std::string>("linkIndex", "auto")
        ->setComment("what a digi-sim link's sim hit index counts: global, perModule or auto (MkFitTraceHitOrigin)");
    desc.addUntracked<std::string>("dumpPrefix", "")
        ->setComment("write <dumpPrefix>hits.txt and simhits.txt for a comparison with the ntuple; empty: none");
    descriptions.add("mkFitTraceTruthNtuple", desc);
  }

  void MkFitTraceTruthNtuple::analyze(edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
    const auto t0 = std::chrono::steady_clock::now();
    edm::Service<MkFitTraceWriter> writer;
    if (!writer.isAvailable())
      throw cms::Exception("Configuration") << "MkFitTraceTruthNtuple needs the MkFitTraceWriter service";

    // Only to run after the producer, whose trace the writer has then taken.
    iEvent.get(mkFitOutputToken_);

    const auto& trackerInfo = iSetup.getData(mkFitGeomToken_).trackerInfo();
    const TrackerGeometry& tracker = iSetup.getData(trackerGeomToken_);
    const TrackingParticleRefVector& tps = iEvent.get(trackingParticlesToken_);
    const SimHitTPAssociationList& simHitsTPAssoc = iEvent.get(simHitTPMapToken_);
    const ClusterTPAssociation& clusterTP = iEvent.get(clusterTPToken_);
    const MkFitClusterIndexToHit& pixelIndex = iEvent.get(pixelClusterIndexToken_);
    const MkFitClusterIndexToHit& stripIndex = iEvent.get(stripClusterIndexToken_);
    const std::vector<int>& pixelLayer = iEvent.get(pixelLayerToken_);
    const std::vector<int>& stripLayer = iEvent.get(stripLayerToken_);

    // ---- The ntuple: TrackingParticles.
    std::unordered_map<unsigned int, int> tpKeyToIndex;
    std::map<int, int> tpsPerBx;
    std::vector<float> tpP(tps.size());
    for (int i = 0; i < (int)tps.size(); ++i) {
      tpKeyToIndex[tps[i].key()] = i;
      ++tpsPerBx[tps[i]->eventId().bunchCrossing()];
      // writeMemoryFile: from the float branches sim_px, sim_py, sim_pz.
      tpP[i] = normNoFma(tps[i]->px(), tps[i]->py(), tps[i]->pz());
    }

    // ---- The ntuple's sim hits (fillSimHits()).
    std::vector<NtSimHit> simHits;
    std::map<std::pair<unsigned int, edm::ProductID>, int> simHitRefKeyToIndex;
    std::vector<TPHitIndex> tpHitList;
    for (const auto& assoc : simHitsTPAssoc) {
      const unsigned int tpKey = assoc.first.key();
      auto found = tpKeyToIndex.find(tpKey);
      if (found == tpKeyToIndex.end())
        continue;
      const PSimHit& simhit = *assoc.second;
      const DetId detId(simhit.detUnitId());
      if (detId.det() != DetId::Tracker)
        continue;
      if (!keepEleSimHits_ && std::abs(simhit.particleType()) == 11 && std::abs(assoc.first->pdgId()) != 11)
        continue;
      const auto key = std::make_pair((unsigned int)assoc.second.key(), assoc.second.id());
      if (!simHitRefKeyToIndex.emplace(key, simHits.size()).second)
        throw cms::Exception("LogicError") << "MkFitTraceTruthNtuple: sim hit " << key.first << " of " << key.second
                                           << " a second time, now with TrackingParticle " << tpKey;
      const auto* det = tracker.idToDetUnit(detId);
      if (!det)
        throw cms::Exception("LogicError") << "MkFitTraceTruthNtuple: no det unit for DetId " << detId.rawId();
      const auto pos = det->surface().toGlobal(simhit.localPosition());
      const auto mom = det->surface().toGlobal(simhit.momentumAtEntry());
      NtSimHit& h = simHits.emplace_back();
      h.x = pos.x();
      h.y = pos.y();
      h.z = pos.z();
      h.px = mom.x();
      h.py = mom.y();
      h.pz = mom.z();
      h.tof = simhit.timeOfFlight();
      h.particle = simhit.particleType();
      h.process = simhit.processType();
      h.tp = found->second;
      h.detId = simhit.detUnitId();
      h.ref = assoc.second;
      tpHitList.push_back(TPHitIndex{tpKey, (unsigned int)(simHits.size() - 1), h.tof, h.detId});
    }

    // ---- The ntuple's matchCluster(), pixel clusters then outer-tracker clusters, in index order.
    struct RecHit {
      int type, idx, layer;
      unsigned int detId;
      const TrackingRecHit* rh;
      std::vector<int> shs;  // the matched sim hits (pix_simHitIdx, ph2_simHitIdx)
    };
    std::vector<RecHit> recHits;
    long long nNoSimHit = 0, nNoLayer = 0;
    auto matchHits = [&](const MkFitClusterIndexToHit& index, const std::vector<int>& layerOf, int type) {
      for (int i = 0; i < (int)index.hits().size(); ++i) {
        const TrackingRecHit* rh = index.hits()[i];
        if (rh == nullptr)
          continue;
        RecHit& r = recHits.emplace_back();
        r.type = type;
        r.idx = i;
        r.layer = i < (int)layerOf.size() ? layerOf[i] : -1;
        nNoLayer += r.layer < 0;
        r.detId = rh->geographicalId().rawId();
        r.rh = rh;
        const OmniClusterRef& cl = static_cast<const BaseTrackerRecHit*>(rh)->firstClusterRef();
        auto range = clusterTP.equal_range(cl);
        for (auto ip = range.first; ip != range.second; ++ip) {
          const TrackingParticleRef& tp = ip->second;
          if (tpKeyToIndex.find(tp.key()) == tpKeyToIndex.end())
            continue;
          const std::pair<TrackingParticleRef, TrackPSimHitRef> dummy(tp, TrackPSimHitRef());
          auto shRange = std::equal_range(simHitsTPAssoc.begin(),
                                          simHitsTPAssoc.end(),
                                          dummy,
                                          SimHitTPAssociationProducer::simHitTPAssociationListGreater);
          const bool eleTp = std::abs(tp->pdgId()) == 11;
          int foundNonElectrons = 0;
          for (auto sp = shRange.first; sp != shRange.second; ++sp)
            if (sp->second->detUnitId() == r.detId && !(std::abs(sp->second->particleType()) == 11 && !eleTp))
              ++foundNonElectrons;
          bool foundSimHit = false, foundElectron = false;
          for (auto sp = shRange.first; sp != shRange.second; ++sp) {
            if (sp->second->detUnitId() != r.detId)
              continue;
            if (std::abs(sp->second->particleType()) == 11 && !eleTp) {
              foundElectron = true;
              if (!keepEleSimHits_)
                continue;
              if (foundNonElectrons > 0)
                continue;
            }
            foundSimHit = true;
            auto f = simHitRefKeyToIndex.find(std::make_pair((unsigned int)sp->second.key(), sp->second.id()));
            if (f == simHitRefKeyToIndex.end())
              throw cms::Exception("LogicError") << "MkFitTraceTruthNtuple: a sim hit of the association is not "
                                                    "among the ntuple's sim hits";
            r.shs.push_back(f->second);
            simHits[f->second].hitIdx.push_back(i);
            simHits[f->second].hitType.push_back(type);
          }
          if (!foundSimHit && !(foundElectron && !keepEleSimHits_))
            ++nNoSimHit;
        }
      }
    };
    matchHits(pixelIndex, pixelLayer, kHitPixel);
    const std::size_t nPixelRecHits = recHits.size();
    matchHits(stripIndex, stripLayer, kHitPhase2OT);

    // The (layer, index) of a rec hit by its type and cluster index, as writeMemoryFile's simHitLayIdx gives them.
    auto layerOfHit = [&](int type, int idx) -> int {
      const std::vector<int>& lo = type == kHitPixel ? pixelLayer : stripLayer;
      return idx >= 0 && idx < (int)lo.size() ? lo[idx] : -1;
    };

    // ---- The TrackingParticles' sim hits in time order, and their first track.
    std::sort(tpHitList.begin(), tpHitList.end(), tpHitIndexListLessSort);
    const edm::Handle<edm::View<reco::Track>> tracksHandle = iEvent.getHandle(tracksToken_);
    edm::RefToBaseVector<reco::Track> trackRefs;
    for (std::size_t i = 0; i < tracksHandle->size(); ++i)
      trackRefs.push_back(tracksHandle->refAt(i));
    const reco::SimToRecoCollection simRecColl = iEvent.get(associatorToken_).associateSimToReco(trackRefs, tps);
    const edm::Handle<TrajectorySeedCollection> seedsHandle = iEvent.getHandle(seedsToken_);
    const TrajectorySeedCollection& seeds = *seedsHandle;

    TraceTruthTables T;
    T.source = source_;
    T.run = iEvent.id().run();
    T.lumi = iEvent.id().luminosityBlock();
    T.event = iEvent.id().event();
    T.module = module_;
    T.layers.resize(trackerInfo.n_layers());
    T.simHitStates.resize(trackerInfo.n_layers());
    for (int l = 0; l < (int)T.layers.size(); ++l)
      T.layers[l].layer = T.simHitStates[l].layer = l;
    T.seeds.resize(seeds.size());

    int nSimsHits = 0, nSimsTrk = 0, nSeedsLab = 0;
    long long nTrkHitsNoLayer = 0;
    T.simTracks.reserve(tps.size());
    for (int s = 0; s < (int)tps.size(); ++s) {
      const TrackingParticle& tp = *tps[s];
      // sim_trkIdx[0] and writeMemoryFile's nlay and seed label.
      int nlay = 0;
      auto foundTracks = simRecColl.find(tps[s]);
      if (foundTracks != simRecColl.end() && !foundTracks->val.empty()) {
        ++nSimsTrk;
        const reco::Track& trk = *foundTracks->val.front().first;
        std::vector<int> hitlay(trackerInfo.n_layers(), 0);
        for (auto ih = trk.recHitsBegin(); ih != trk.recHitsEnd(); ++ih) {
          const TrackingRecHit* hit = *ih;
          if (hit->geographicalId().det() != DetId::Tracker || !hit->isValid())
            continue;
          const OmniClusterRef& cl = static_cast<const BaseTrackerRecHit*>(hit)->firstClusterRef();
          int lay = -1;
          if (cl.isPixel())
            lay = layerOfHit(kHitPixel, cl.cluster_pixel().key());
          else if (cl.isPhase2())
            lay = layerOfHit(kHitPhase2OT, cl.cluster_phase2OT().key());
          if (lay >= 0 && lay < (int)hitlay.size())
            ++hitlay[lay];
          else
            ++nTrkHitsNoLayer;
        }
        for (int c : hitlay)
          nlay += c > 0;
        if (trk.seedRef().id() == seedsHandle.id()) {
          const auto key = trk.seedRef().key();
          if (key >= T.seeds.size())
            throw cms::Exception("LogicError")
                << "MkFitTraceTruthNtuple: a track's seed " << key << " of " << T.seeds.size();
          T.seeds[key].sim = s;
        }
      }
      // writeMemoryFile's sim track, from the ntuple's floats.
      const auto& v = tp.parentVertex()->position();
      const float vx = v.x(), vy = v.y(), vz = v.z();
      const float px = tp.px(), py = tp.py(), pz = tp.pz();
      const mkfit::SVector3 pos(vx, vy, vz);
      const mkfit::SVector3 mom(px, py, pz);
      mkfit::SMatrixSym66 err;
      err.At(0, 0) = vx * vx;
      err.At(1, 1) = vy * vy;
      err.At(2, 2) = vz * vz;
      err.At(3, 3) = px * px;
      err.At(4, 4) = py * py;
      err.At(5, 5) = pz * pz;
      mkfit::TrackState state((int)tp.charge(), pos, mom, err);
      state.convertFromCartesianToCCS();
      mkfit::Track& t = T.simTracks.emplace_back(state, float(nlay), s, 0, nullptr);
      if (tp.eventId().bunchCrossing() == 0)
        t.setProdType(tp.eventId().event() == 0 ? mkfit::Track::ProdType::Signal : mkfit::Track::ProdType::InTimePU);
      else
        t.setProdType(mkfit::Track::ProdType::OutOfTimePU);
      auto range = std::equal_range(
          tpHitList.begin(), tpHitList.end(), TPHitIndex{(unsigned int)tps[s].key(), 0, 0.f, 0}, tpHitIndexListLess);
      for (auto ip = range.first; ip != range.second; ++ip) {
        const NtSimHit& h = simHits[ip->simHitIdx];
        for (std::size_t k = 0; k < h.hitIdx.size(); ++k) {
          const int lay = layerOfHit(h.hitType[k], h.hitIdx[k]);
          if (lay >= 0)
            t.addHitIdx(h.hitIdx[k], lay, 0.f);
        }
      }
      nSimsHits += t.nTotalHits() > 0;
    }

    // ---- The seeds: their hits, and the labels set above.
    for (int i = 0; i < (int)seeds.size(); ++i) {
      auto& ss = T.seeds[i];
      ss.seed = i;
      nSeedsLab += ss.sim >= 0;
      for (const TrackingRecHit& hit : seeds[i].recHits()) {
        if (!hit.isValid())
          continue;
        const OmniClusterRef& cl = static_cast<const BaseTrackerRecHit&>(hit).firstClusterRef();
        int idx = -1, lay = -1;
        if (cl.isPixel()) {
          idx = cl.cluster_pixel().key();
          lay = layerOfHit(kHitPixel, idx);
        } else if (cl.isPhase2()) {
          idx = cl.cluster_phase2OT().key();
          lay = layerOfHit(kHitPhase2OT, idx);
        }
        ss.hit_lyr.push_back(lay);
        ss.hit_idx.push_back(idx);
      }
    }

    // ---- The link check: the sim hit vectors of the whole crossing, per (module, event id) their sim hits in order.
    const bool checkLinks = !simHitTokens_.empty();
    std::vector<edm::Handle<std::vector<PSimHit>>> vecs;
    std::vector<std::unordered_map<std::uint64_t, std::vector<int>>> hitsOfDetEvent(simHitTokens_.size());
    for (std::size_t c = 0; c < simHitTokens_.size(); ++c) {
      vecs.push_back(iEvent.getHandle(simHitTokens_[c]));
      const std::vector<PSimHit>& v = *vecs.back();
      for (int j = 0; j < (int)v.size(); ++j)
        hitsOfDetEvent[c][detKey(v[j].detUnitId(), v[j].eventId().rawId())].push_back(j);
    }
    // The sim hit a link names under a decoding: (vector, index) or (-1, -1).
    auto resolve = [&](unsigned int detId, const PixelDigiSimLink& l, bool perModule) -> std::pair<int, int> {
      const unsigned int pos = l.CFposition();
      const std::uint32_t eid = l.eventId().rawId();
      for (int c = 0; c < (int)vecs.size(); ++c) {
        if (simHitTof_[c] != l.TofBin())
          continue;
        const std::vector<PSimHit>& v = *vecs[c];
        auto hs = hitsOfDetEvent[c].find(detKey(detId, eid));
        if (hs == hitsOfDetEvent[c].end())
          continue;
        int j = -1;
        if (perModule) {
          const long long k = (long long)pos - hs->second.front();
          if (k >= 0 && k < (long long)hs->second.size())
            j = hs->second[k];
        } else if (pos < v.size() && v[pos].detUnitId() == detId && v[pos].eventId().rawId() == eid) {
          j = pos;
        }
        if (j >= 0 && v[j].trackId() == l.SimTrackId())
          return {c, j};
      }
      return {-1, -1};
    };
    const SimLinks& pixelLinks = iEvent.get(pixelSimLinkToken_);
    const SimLinks& phase2Links = iEvent.get(phase2OTSimLinkToken_);
    long long nEvLinks = 0, nFound[2] = {0, 0};
    if (checkLinks)
      for (const SimLinks* links : {&pixelLinks, &phase2Links})
        for (const auto& ds : *links)
          for (const PixelDigiSimLink& l : ds.data) {
            ++nEvLinks;
            for (int pm = 0; pm < 2; ++pm)
              nFound[pm] += resolve(ds.detId(), l, pm == 1).first >= 0;
          }
    const bool perModule =
        linkIndex_ == LinkIndex::kPerModule || (linkIndex_ == LinkIndex::kAuto && nFound[1] > nFound[0]);
    std::map<edm::ProductID, int> vecOfId;
    for (int c = 0; c < (int)vecs.size(); ++c)
      vecOfId[vecs[c].id()] = c;
    // Whether a link on the rec hit's channels names sim hit sh, under each decoding.
    auto linked = [&](const RecHit& r, const NtSimHit& h, bool pm) -> bool {
      auto fv = vecOfId.find(h.ref.id());
      if (fv == vecOfId.end())
        return false;
      const SimLinks& links = r.type == kHitPixel ? pixelLinks : phase2Links;
      auto ds = links.find(r.detId);
      if (ds == links.end())
        return false;
      std::vector<unsigned int> channels;
      const OmniClusterRef& cl = static_cast<const BaseTrackerRecHit*>(r.rh)->firstClusterRef();
      if (r.type == kHitPixel) {
        const SiPixelCluster& c = cl.pixelCluster();
        for (int ip = 0; ip < c.size(); ++ip)
          channels.push_back(PixelChannelIdentifier::pixelToChannel(c.pixel(ip).x, c.pixel(ip).y));
      } else {
        const Phase2TrackerCluster1D& c = cl.phase2OTCluster();
        for (unsigned int is = 0; is < c.size(); ++is)
          channels.push_back(Phase2TrackerDigi::pixelToChannel(c.firstRow() + is, c.column()));
      }
      for (const PixelDigiSimLink& l : ds->data) {
        if (std::find(channels.begin(), channels.end(), l.channel()) == channels.end())
          continue;
        const auto [c, j] = resolve(r.detId, l, pm);
        if (c == fv->second && j == (int)h.ref.key())
          return true;
      }
      return false;
    };

    // ---- writeMemoryFile's bestTkIdx() per rec hit, the labels and the chosen sim hits.
    long long nHitsSimHit = 0, nLabels = 0, nLostP50 = 0, nLostHalf = 0, nLostArb = 0, nStates = 0;
    long long nChecked = 0, nLinked[2] = {0, 0};
    for (std::size_t r = 0; r < recHits.size(); ++r) {
      const RecHit& rh = recHits[r];
      if (rh.layer < 0)
        continue;
      auto& L = T.layers.at(rh.layer);
      L.n_hits = std::max(L.n_hits, rh.idx + 1);
      auto& S = T.simHitStates.at(rh.layer);
      S.n_hits = std::max(S.n_hits, rh.idx + 1);
      if (rh.shs.empty())
        continue;
      ++nHitsSimHit;
      int ibest = -1, shbest = -1;
      float maxfrac = -1;
      bool anyP50 = false;
      for (int sh : rh.shs) {
        const NtSimHit& h = simHits[sh];
        const int tkidx = h.tp;
        if (tkidx < 0)
          continue;
        const float hp = normNoFma(h.px, h.py, h.pz);
        if (hp < 0.05f)
          continue;
        anyP50 = true;
        const float tp = tpP[tkidx];
        if (hp < 0.5 * tp)
          continue;
        if (maxfrac < hp / tp) {
          maxfrac = hp / tp;
          ibest = tkidx;
          shbest = sh;
        }
      }
      if (shbest < 0) {
        ++(anyP50 ? nLostHalf : nLostP50);
      }
      // The split-cluster arbitration.
      if (rh.shs.size() == 1 && ibest >= 0) {
        const NtSimHit& h = simHits[shbest];
        for (std::size_t ih = 0; ih < h.hitType.size(); ++ih) {
          if (h.hitType[ih] == rh.type) {
            if (h.hitIdx[ih] != rh.idx) {
              ibest = -1;
              ++nLostArb;
            }
            break;
          }
        }
      }
      if (ibest >= 0) {
        ++nLabels;
        L.hit.push_back(rh.idx);
        L.sim.push_back(ibest);
      }
      if (shbest >= 0) {
        ++nStates;
        const NtSimHit& h = simHits[shbest];
        S.hit.push_back(rh.idx);
        S.sh_idx.push_back(shbest);
        S.x.push_back(h.x);
        S.y.push_back(h.y);
        S.z.push_back(h.z);
        S.px.push_back(h.px);
        S.py.push_back(h.py);
        S.pz.push_back(h.pz);
        mkfit::TrackState ts(
            (int)tps[h.tp]->charge(), mkfit::SVector3(h.x, h.y, h.z), mkfit::SVector3(h.px, h.py, h.pz), {});
        ts.convertFromCartesianToCCS();
        S.ts_q.push_back(ts.charge);
        S.ts_par.push_back({{ts.parameters[0],
                             ts.parameters[1],
                             ts.parameters[2],
                             ts.parameters[3],
                             ts.parameters[4],
                             ts.parameters[5]}});
        std::int8_t lk = -1;
        if (checkLinks) {
          ++nChecked;
          const bool l0 = linked(rh, h, false), l1 = linked(rh, h, true);
          nLinked[0] += l0;
          nLinked[1] += l1;
          lk = (perModule ? l1 : l0) ? 1 : 0;
        }
        S.linked.push_back(lk);
      }
    }

    if (!dumpPrefix_.empty()) {
      std::ostringstream key;
      key << T.run << "_" << T.lumi << "_" << T.event;
      std::ofstream fh(dumpPrefix_ + key.str() + "-hits.txt"), fs(dumpPrefix_ + key.str() + "-simhits.txt");
      std::vector<char> onHit(simHits.size(), 0);
      for (const RecHit& rh : recHits) {
        if (rh.shs.empty())
          continue;
        fh << rh.type << " " << rh.idx << " " << rh.layer;
        for (int sh : rh.shs) {
          fh << " " << sh;
          onHit[sh] = 1;
        }
        fh << "\n";
      }
      char buf[200];
      for (std::size_t sh = 0; sh < simHits.size(); ++sh) {
        if (!onHit[sh])
          continue;
        const NtSimHit& h = simHits[sh];
        std::snprintf(buf, sizeof(buf), "%zu %u %.9g %.9g %.9g %d\n", sh, h.detId, h.x, h.y, h.z, h.tp);
        fs << buf;
      }
    }

    const int pass = writer->writeTruth(T);

    std::ostringstream bx;
    for (const auto& [b, n] : tpsPerBx)
      bx << (bx.tellp() > 0 ? "," : "") << b << ":" << n;
    char buf[1024];
    std::snprintf(buf,
                  sizeof(buf),
                  "MKFITTRUTHNT %u:%u:%llu %s source=%s tps=%zu tps_bx=%s simhits=%zu hits=%zu (pixel %zu) "
                  "hits_nolayer=%lld hits_simhit=%lld labels=%lld lost_p50=%lld lost_half=%lld lost_arb=%lld "
                  "states=%lld nosimhit=%lld sims_hits=%d sims_trk=%d trk_hits_nolayer=%lld seeds=%zu seeds_lab=%d "
                  "link_index=%s linked=%lld/%lld linked_other=%lld event_links=%lld found_global=%lld "
                  "found_perModule=%lld pass=%d ms=%.1f",
                  T.run,
                  T.lumi,
                  T.event,
                  module_.c_str(),
                  source_.c_str(),
                  tps.size(),
                  bx.str().c_str(),
                  simHits.size(),
                  recHits.size(),
                  nPixelRecHits,
                  nNoLayer,
                  nHitsSimHit,
                  nLabels,
                  nLostP50,
                  nLostHalf,
                  nLostArb,
                  nStates,
                  nNoSimHit,
                  nSimsHits,
                  nSimsTrk,
                  nTrkHitsNoLayer,
                  seeds.size(),
                  nSeedsLab,
                  checkLinks ? (perModule ? "perModule" : "global") : "none",
                  nLinked[perModule ? 1 : 0],
                  nChecked,
                  nLinked[perModule ? 0 : 1],
                  nEvLinks,
                  nFound[0],
                  nFound[1],
                  pass,
                  std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
    std::cout << buf << std::endl;
  }

}  // namespace mkfit_fwff

using mkfit_fwff::MkFitTraceTruthNtuple;
DEFINE_FWK_MODULE(MkFitTraceTruthNtuple);
