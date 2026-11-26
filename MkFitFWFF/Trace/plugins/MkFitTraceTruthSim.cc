// Package:    MkFitFWFF/Trace
// Class:      MkFitTraceTruthSim
//
// The truth source "sim" of one MkFitProducer's traced event, made straight from the simulation's products: the
// SimTracks, SimVertices and PSimHits of the whole crossing (signal and pileup, every bunch crossing) and the digi-sim
// links, with no TrackingParticle in the labels, the sim tracks, their states or the seed labels. It gives the tables
// to the MkFitTraceWriter service as the truth source named by "source" (TraceRNTupleWriter.h: TruthHits_<source>,
// SimTracks_<source>, TruthSeeds_<source>, SimHitStates_<source>, SimTrackMap_<source>), and optionally a second
// source, sourceAncestor, with the ancestor rule below.
//
// The simulation's collections. Geant4 keeps only some particles as SimTracks. A sim hit of a particle that is not
// kept carries the id of a kept ancestor, with its own particle type and creating process (PSimHit::particleType(),
// processType()); TrackingParticles and ClusterTPAssociation go by the id alone, so they give such a deposit (a delta
// ray, a conversion electron, the products of a nuclear interaction or of a decay) to the ancestor. Here:
// - The sim hits are the vectors of the whole crossing (simHits: MkFitTraceCrossingFramePSimHits of the mixing's
//   playback, the signal first, then the pileup events as the mixing accumulated them), the SimTracks and SimVertices
//   likewise (simTracks, simVertices). A SimTrack is (event id, track id); its creating process is the processType()
//   of its SimVertex, whose index (vertIndex()) counts in the SimVertices of its own event.
// - A sim hit is own when a SimTrack with its (event id, track id) exists and has the sim hit's particle type and
//   creating process; otherwise it is secondary (an orphan when no SimTrack has its id).
// - A hit is (layer, index): the index is the cluster's index in its cluster collection (MkFitHitWrapper), as the
//   records keep it. Its charge: per pixel of the cluster the pixel's ADC, per strip of an outer-tracker cluster 1.
//   Each digi-sim link on a channel of the cluster gives its share of the channel's charge (the shares scaled to 1
//   when the noise made their sum exceed 1) to the sim hit it names. The links' sim hit index is decoded as the
//   event's links agree with (linkIndex auto: global, or per module as the digitizer of CMSSW_20_0_0_pre1 counted;
//   MkFitTraceHitOrigin's comment), in the vector of the link's tof bin, the sim hit's event id and track id checked.
//   A link whose sim hit is not found is unmatched charge.
// - The label ("sim"): the SimTrack whose own sim hits put the most charge into the cluster (charged SimTracks only:
//   SimTrack::charge(), or Z for a nucleus, to which CoreSimTrack gives 0; on a tie the one first in the crossing's
//   SimTracks), provided that charge is at least minOwnShare of the cluster's charge. A hit whose cluster holds only
//   secondary charge is not labelled: the particle whose id it carries did not deposit it, and the deposit lies where
//   the secondary went, often off the particle's path or on a module the particle never crossed. A stored secondary
//   (a conversion electron kept as a SimTrack, ...) labels its own hits.
// - The ancestor rule (sourceAncestor, empty: not written): the SimTrack with the most charge of the cluster, own and
//   secondary together (the secondary charge counts for the SimTrack whose id it carries, as TrackingParticles count
//   it), any charge; orphan charge counts for none.
// - The sim tracks, the same list for both sources: every charged SimTrack with at least one own sim hit in the
//   tracker (the simHits vectors), and every SimTrack the ancestor rule labels a hit with; in the order of the
//   crossing's SimTracks. A sim track's label is its index. Its state at production: the position of its SimVertex,
//   the SimTrack's momentum, its charge, as floats, with the placeholder covariance of the other sources (each
//   variance the square of its value), in CCS; the production type from the event id (signal: bunch crossing 0,
//   event 0; in-time pileup: bunch crossing 0; out-of-time pileup otherwise). Its hits: the hits it labels, in the time
//   of flight of their sim hit below, then layer and index. chi2: as the source "ntuple" gives it, the number of
//   mkFit layers with a valid hit on the first track (tracks, by trackAssociator's associateSimToReco) of the
//   TrackingParticle of the tpSources entry named by chi2Source that lists the SimTrack among its g4Tracks; 0 without
//   one.
// - SimHitStates_<source>, per labelled hit: the label's own sim hit with the most charge in the cluster (the
//   ancestor rule: its sim hit with the most charge, own first), its global position (the local mid point) and
//   momentum at entry, the state there (the sim track's charge, CCS); sh_idx counts the sim hits over the simHits
//   vectors in their order; linked 1.
// - TruthSeeds_<source>, per input seed (with its hits): sim and sim_frac by the rule of the ntuple's seed
//   association (quickTrackAssociatorByHits on the seed as a track: the sim track with the most of the seed's hits,
//   when their fraction exceeds 0.75, and for a seed of 3 hits all 3), best_sim and best_frac by the ntuple's
//   see_bestSimTrkIdx (the sim track with the most of the seed's hits, at least 3, the innermost first on a tie); a hit
//   counts for its label only.
// - SimTrackMap_<source>: per sim track its SimTrack (track id, event id, type, creating process, parent's track id,
//   generator index, own sim hits in the tracker) and, per entry of tpSources, the index of the TrackingParticle
//   listing the SimTrack among its g4Tracks in that source's sim tracks (-1 none). MergedTrackTruth merges an
//   electron's bremsstrahlung (and other one-electron vertices of an electron) into one TrackingParticle with several
//   g4Tracks, so several sim tracks can map to one TrackingParticle.
//
// Prints one line per event:
//   MKFITTRUTHSIM <run>:<lumi>:<event> <module> source=<s> ancestor=<s or -> simtracks=<SimTracks of the crossing>
//     simhits=<sim hits> simhits_own=.. simhits_sec=.. simhits_orph=.. sims=<sim tracks> sims_own=<with own sim hits>
//     sims_anc_only=<only by the ancestor rule> sims_prod=<signal>/<in-time>/<out-of-time> sims_hits=<with a "sim" hit>
//     sims_tp=<with a TrackingParticle per tpSources entry> hits=<registered> hits_linked=<with a link>
//     labels=<"sim" labels> labels_anc=<ancestor labels> sec_only=<hits with charge only from secondaries>
//     below_share=<own charge below minOwnShare> neutral_own=<hits where a neutral SimTrack had the most own charge>
//     q_own=.. q_sec=.. q_orph=.. q_unm=.. q_noise=.. (shares of the listed hits' charge) links=<on hits>
//     links_unmatched=<n> link_index=<..> event_links=<n> found_global=<n> found_perModule=<n> seeds=<n>
//     seeds_sim=<n> seeds_best=<n> tp_conflicts=<SimTracks claimed by two TrackingParticles of one source>
//     pass=<p> ms=<time>
// With dumpPrefix, CSV files for the comparison of the sources (<dumpPrefix><run>_<lumi>_<event>_p<pass>-*.csv):
//   hits.csv      per hit with a link: layer, idx, detid, pixel, the cluster's charge and its noise and unmatched
//                 parts, the rec hit's local x, y and global x, y, z, the number of SimTracks with charge, the two
//                 labels (sim track indices), the own share of the "sim" label, the shares of all own, secondary and
//                 orphan charge.
//   hitkeys.csv   per (hit, SimTrack with charge in it): the SimTrack (crossing index, event id, track id, type,
//                 charge, sim track index or -1, TrackingParticle per tpSources entry or -1), its own and secondary
//                 charge, its own sim hit with the most charge in the cluster (local x, y, global x, y, z, |p| at
//                 entry, time of flight; -1 when none) or else its own sim hit on the module nearest the rec hit
//                 (flag own_in 1 in the cluster, 0 on the module, -1 none), its own sim hits on the module, and its
//                 secondary sim hit with the most charge (type, process, |p|).
//   simtracks.csv per sim track: index, SimTrack, state at production (vertex, momentum), charge, type, process,
//                 parent, the TrackingParticles, own sim hits in the tracker, hits labelled in each source, chi2.
//   tps.csv       per TrackingParticle of each tpSources entry: source, index, event id, type, charge, its parent
//                 vertex and momentum, its g4Tracks' crossing indices.

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
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PixelChannelIdentifier.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

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
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mkfit_fwff {

  class MkFitTraceTruthSim : public edm::global::EDAnalyzer<> {
  public:
    explicit MkFitTraceTruthSim(const edm::ParameterSet& iConfig);

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void analyze(edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

    using SimLinks = edm::DetSetVector<PixelDigiSimLink>;
    enum class LinkIndex { kAuto, kGlobal, kPerModule };

    const std::string module_;
    const std::string source_;
    const std::string sourceAncestor_;
    const edm::EDGetTokenT<MkFitOutputWrapper> mkFitOutputToken_;
    const edm::EDGetTokenT<MkFitClusterIndexToHit> pixelClusterIndexToken_;
    const edm::EDGetTokenT<MkFitClusterIndexToHit> stripClusterIndexToken_;
    const edm::EDGetTokenT<std::vector<int>> pixelLayerToken_;
    const edm::EDGetTokenT<std::vector<int>> stripLayerToken_;
    const edm::EDGetTokenT<TrajectorySeedCollection> seedsToken_;
    const edm::EDGetTokenT<SimLinks> pixelSimLinkToken_;
    const edm::EDGetTokenT<SimLinks> phase2OTSimLinkToken_;
    std::vector<edm::EDGetTokenT<std::vector<PSimHit>>> simHitTokens_;
    std::vector<unsigned int> simHitTof_;  // PixelDigiSimLink::LowTof or HighTof, from the instance name
    const edm::EDGetTokenT<edm::SimTrackContainer> simTracksToken_;
    const edm::EDGetTokenT<edm::SimVertexContainer> simVerticesToken_;
    // The TrackingParticles of other sources, for SimTrackMap and chi2.
    std::vector<std::string> tpNames_;
    std::vector<edm::EDGetTokenT<TrackingParticleRefVector>> tpTokens_;
    int chi2Source_ = -1;  // the tpSources entry whose TrackingParticles give chi2, -1 none
    edm::EDGetTokenT<edm::View<reco::Track>> tracksToken_;
    edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> associatorToken_;
    const double minOwnShare_;
    LinkIndex linkIndex_;
    const std::string dumpPrefix_;
    const edm::ESGetToken<MkFitGeometry, TrackerRecoGeometryRecord> mkFitGeomToken_;
    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeomToken_;
  };

  namespace {
    std::uint64_t idKey(std::uint32_t eventId, unsigned int trackId) {
      return (std::uint64_t(eventId) << 32) | trackId;
    }
    std::uint64_t detKey(unsigned int detId, unsigned int x) { return (std::uint64_t(detId) << 32) | x; }

    enum ShClass : std::int8_t { kOwn = 0, kSec = 1, kOrph = 2 };

    constexpr int kHitPixel = 0, kHitPhase2OT = 4;  // TrackingNtuple's HitType

    // A sim hit as (vector, index); packed for maps.
    struct ShRef {
      int c = -1, j = -1;
      bool valid() const { return c >= 0; }
    };

    // A SimTrack's charge in one cluster.
    struct KeyCharge {
      int st = -1;              // SimTrack index in the crossing
      double own = 0, sec = 0;  // its own and secondary charge
      ShRef bestOwn, bestAny;   // its own sim hit with the most charge, and its sim hit with the most charge
      double bestOwnQ = -1, bestAnyQ = -1;
      ShRef bestSec;  // its secondary sim hit with the most charge
      double bestSecQ = -1;
    };

    // A SimTrack's charge: SimTrack::charge() (CoreSimTrack: from the PDG table, 0 for a nucleus), and for a nucleus
    // (PDG code 10LZZZAAAI) its Z, so that a deuteron's or an alpha's own ionisation can label hits.
    int chargeOf(const SimTrack& t) {
      const int a = std::abs(t.type());
      if (a >= 1000000000)
        return (t.type() > 0 ? 1 : -1) * ((a / 10000) % 1000);
      return (int)std::lround(t.charge());
    }

    // Charge per sim hit within one cluster.
    struct ShCharge {
      ShRef r;
      double q = 0;
    };
  }  // namespace

  MkFitTraceTruthSim::MkFitTraceTruthSim(const edm::ParameterSet& iConfig)
      : module_(iConfig.getParameter<std::string>("mkFitProducer")),
        source_(iConfig.getParameter<std::string>("source")),
        sourceAncestor_(iConfig.getParameter<std::string>("sourceAncestor")),
        mkFitOutputToken_(consumes(edm::InputTag(module_))),
        pixelClusterIndexToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelHits"))),
        stripClusterIndexToken_(consumes(iConfig.getParameter<edm::InputTag>("stripHits"))),
        pixelLayerToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelHits"))),
        stripLayerToken_(consumes(iConfig.getParameter<edm::InputTag>("stripHits"))),
        seedsToken_(consumes(iConfig.getParameter<edm::InputTag>("seeds"))),
        pixelSimLinkToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelSimLink"))),
        phase2OTSimLinkToken_(consumes(iConfig.getParameter<edm::InputTag>("phase2OTSimLink"))),
        simTracksToken_(consumes(iConfig.getParameter<edm::InputTag>("simTracks"))),
        simVerticesToken_(consumes(iConfig.getParameter<edm::InputTag>("simVertices"))),
        minOwnShare_(iConfig.getParameter<double>("minOwnShare")),
        dumpPrefix_(iConfig.getUntrackedParameter<std::string>("dumpPrefix")),
        mkFitGeomToken_(esConsumes()),
        trackerGeomToken_(esConsumes()) {
    if (source_.empty() || source_ == "tp" || sourceAncestor_ == "tp" || sourceAncestor_ == source_)
      throw cms::Exception("Configuration")
          << "MkFitTraceTruthSim: source must be a name other than 'tp', and sourceAncestor empty or another name";
    for (const edm::InputTag& t : iConfig.getParameter<std::vector<edm::InputTag>>("simHits")) {
      simHitTokens_.push_back(consumes<std::vector<PSimHit>>(t));
      simHitTof_.push_back(t.instance().find("HighTof") != std::string::npos ? PixelDigiSimLink::HighTof
                                                                             : PixelDigiSimLink::LowTof);
    }
    if (simHitTokens_.empty())
      throw cms::Exception("Configuration") << "MkFitTraceTruthSim: simHits must name the sim hit vectors";
    const std::string chi2From = iConfig.getParameter<std::string>("chi2Source");
    for (const edm::ParameterSet& p : iConfig.getParameter<std::vector<edm::ParameterSet>>("tpSources")) {
      tpNames_.push_back(p.getParameter<std::string>("name"));
      tpTokens_.push_back(consumes<TrackingParticleRefVector>(p.getParameter<edm::InputTag>("trackingParticles")));
      if (tpNames_.back() == chi2From)
        chi2Source_ = tpNames_.size() - 1;
    }
    if (!chi2From.empty()) {
      if (chi2Source_ < 0)
        throw cms::Exception("Configuration")
            << "MkFitTraceTruthSim: chi2Source '" << chi2From << "' is not a tpSources name";
      tracksToken_ = consumes(iConfig.getParameter<edm::InputTag>("tracks"));
      associatorToken_ = consumes(iConfig.getParameter<edm::InputTag>("trackAssociator"));
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
          << "MkFitTraceTruthSim: linkIndex must be auto, global or perModule, not '" << li << "'";
  }

  void MkFitTraceTruthSim::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("mkFitProducer", "initialStepTrackCandidatesMkFit")
        ->setComment(
            "label of the traced MkFitProducer: the key of the truth, and its output orders this module after it");
    desc.add<std::string>("source", "sim")->setComment("the truth source's name in the trace file");
    desc.add<std::string>("sourceAncestor", "")
        ->setComment("the name of a second source with the ancestor rule (secondary charge counts); empty: none");
    desc.add<edm::InputTag>("pixelHits", edm::InputTag("mkFitSiPixelHits"));
    desc.add<edm::InputTag>("stripHits", edm::InputTag("mkFitSiPhase2Hits"));
    desc.add<edm::InputTag>("seeds", edm::InputTag("initialStepSeeds"))->setComment("the producer's input seeds");
    desc.add<edm::InputTag>("pixelSimLink", edm::InputTag("simSiPixelDigis", "Pixel"));
    desc.add<edm::InputTag>("phase2OTSimLink", edm::InputTag("simSiPixelDigis", "Tracker"));
    desc.add<std::vector<edm::InputTag>>("simHits", {})
        ->setComment(
            "the tracker sim hit vectors of the whole crossing, signal first (MkFitTraceCrossingFramePSimHits), one "
            "per collection the digitizer read; the tof bin from the instance name");
    desc.add<edm::InputTag>("simTracks", edm::InputTag("mkFitTracePuSimTracks", "g4SimHits"))
        ->setComment("the SimTracks of the whole crossing, signal first (MkFitTraceCrossingFrameSimTracks)");
    desc.add<edm::InputTag>("simVertices", edm::InputTag("mkFitTracePuSimVertices", "g4SimHits"))
        ->setComment("the SimVertices of the whole crossing, signal first (MkFitTraceCrossingFrameSimVertices)");
    edm::ParameterSetDescription tpDesc;
    tpDesc.add<std::string>("name", "")->setComment("the other source's name: the field tp_<name> of SimTrackMap");
    tpDesc.add<edm::InputTag>("trackingParticles", edm::InputTag(""))
        ->setComment("TrackingParticleRefVector: that source's sim tracks, in order");
    desc.addVPSet("tpSources", tpDesc, {})
        ->setComment("other sources whose sim tracks are TrackingParticles, for SimTrackMap; empty: no map");
    desc.add<std::string>("chi2Source", "")
        ->setComment(
            "the tpSources entry whose TrackingParticles' first track gives chi2 (as the source ntuple); empty: 0");
    desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"))->setComment("the ntuple's tracks, for chi2");
    desc.add<edm::InputTag>("trackAssociator", edm::InputTag("quickTrackAssociatorByHits"));
    desc.add<double>("minOwnShare", 0.)
        ->setComment("the least share of the cluster's charge the label's own charge must have; 0: any own charge");
    desc.add<std::string>("linkIndex", "auto")
        ->setComment("what a digi-sim link's sim hit index counts: global, perModule or auto (MkFitTraceHitOrigin)");
    desc.addUntracked<std::string>("dumpPrefix", "")
        ->setComment("write the CSV files for comparing the sources (hits, hitkeys, simtracks, tps); empty: none");
    descriptions.add("mkFitTraceTruthSim", desc);
  }

  void MkFitTraceTruthSim::analyze(edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
    const auto t0 = std::chrono::steady_clock::now();
    edm::Service<MkFitTraceWriter> writer;
    if (!writer.isAvailable())
      throw cms::Exception("Configuration") << "MkFitTraceTruthSim needs the MkFitTraceWriter service";

    // Only to run after the producer, whose trace the writer has then taken.
    iEvent.get(mkFitOutputToken_);

    const auto& trackerInfo = iSetup.getData(mkFitGeomToken_).trackerInfo();
    const TrackerGeometry& tracker = iSetup.getData(trackerGeomToken_);
    const MkFitClusterIndexToHit& pixelIndex = iEvent.get(pixelClusterIndexToken_);
    const MkFitClusterIndexToHit& stripIndex = iEvent.get(stripClusterIndexToken_);
    const std::vector<int>& pixelLayer = iEvent.get(pixelLayerToken_);
    const std::vector<int>& stripLayer = iEvent.get(stripLayerToken_);
    const bool dumping = !dumpPrefix_.empty();
    const bool withAnc = !sourceAncestor_.empty();

    // ---- SimTracks and SimVertices of the crossing: by (event id, track id); per event its SimVertices in order.
    const edm::SimTrackContainer& simTracks = iEvent.get(simTracksToken_);
    const edm::SimVertexContainer& simVertices = iEvent.get(simVerticesToken_);
    const int nSt = simTracks.size();
    std::unordered_map<std::uint64_t, int> stOfKey;
    stOfKey.reserve(2 * nSt);
    for (int i = 0; i < nSt; ++i)
      stOfKey.emplace(idKey(simTracks[i].eventId().rawId(), simTracks[i].trackId()), i);
    std::unordered_map<std::uint32_t, std::vector<const SimVertex*>> vertsOf;
    for (const SimVertex& v : simVertices)
      vertsOf[v.eventId().rawId()].push_back(&v);
    std::vector<const SimVertex*> stVertex(nSt, nullptr);
    std::vector<int> stProc(nSt, -1), stCharge(nSt, 0);
    long long nVertBad = 0;
    for (int i = 0; i < nSt; ++i) {
      const SimTrack& t = simTracks[i];
      stCharge[i] = chargeOf(t);
      auto f = vertsOf.find(t.eventId().rawId());
      const int vi = t.vertIndex();
      if (f != vertsOf.end() && vi >= 0 && vi < (int)f->second.size()) {
        stVertex[i] = f->second[vi];
        stProc[i] = (int)stVertex[i]->processType();
      } else {
        ++nVertBad;
      }
    }

    // ---- The sim hit vectors: per (module, event id) their sim hits in order; per sim hit its SimTrack and class.
    std::vector<edm::Handle<std::vector<PSimHit>>> vecs;
    std::vector<long long> vecOffset;  // the first sh_idx of each vector
    std::vector<std::unordered_map<std::uint64_t, std::vector<int>>> hitsOfDetEvent(simHitTokens_.size());
    std::vector<std::vector<int>> shSt(simHitTokens_.size());
    std::vector<std::vector<std::int8_t>> shCls(simHitTokens_.size());
    std::vector<int> nOwnSh(nSt, 0);
    // Own sim hits per (module, SimTrack), for the dump's nearest own sim hit.
    std::unordered_map<std::uint64_t, std::vector<ShRef>> ownOnDet;
    long long nSimHits = 0, nShCls[3] = {0, 0, 0};
    for (std::size_t c = 0; c < simHitTokens_.size(); ++c) {
      vecOffset.push_back(nSimHits);
      vecs.push_back(iEvent.getHandle(simHitTokens_[c]));
      const std::vector<PSimHit>& v = *vecs.back();
      nSimHits += v.size();
      shSt[c].assign(v.size(), -1);
      shCls[c].assign(v.size(), kOrph);
      for (int j = 0; j < (int)v.size(); ++j) {
        const PSimHit& h = v[j];
        hitsOfDetEvent[c][detKey(h.detUnitId(), h.eventId().rawId())].push_back(j);
        auto f = stOfKey.find(idKey(h.eventId().rawId(), h.trackId()));
        if (f != stOfKey.end()) {
          const int st = f->second;
          shSt[c][j] = st;
          const bool own = h.particleType() == simTracks[st].type() && (int)h.processType() == stProc[st];
          shCls[c][j] = own ? kOwn : kSec;
          if (own) {
            ++nOwnSh[st];
            if (dumping)
              ownOnDet[detKey(h.detUnitId(), st)].push_back(ShRef{(int)c, j});
          }
        }
        ++nShCls[shCls[c][j]];
      }
    }

    // ---- The links' sim hit index: the sim hit a link names under a decoding, as MkFitTraceTruthNtuple resolves it.
    auto resolve = [&](unsigned int detId, const PixelDigiSimLink& l, bool perModule) -> ShRef {
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
          return ShRef{c, j};
      }
      return ShRef{};
    };
    const SimLinks& pixelLinks = iEvent.get(pixelSimLinkToken_);
    const SimLinks& phase2Links = iEvent.get(phase2OTSimLinkToken_);
    long long nEvLinks = 0, nFound[2] = {0, 0};
    for (const SimLinks* links : {&pixelLinks, &phase2Links})
      for (const auto& ds : *links)
        for (const PixelDigiSimLink& l : ds.data) {
          ++nEvLinks;
          for (int pm = 0; pm < 2; ++pm)
            nFound[pm] += resolve(ds.detId(), l, pm == 1).valid();
        }
    const bool perModule =
        linkIndex_ == LinkIndex::kPerModule || (linkIndex_ == LinkIndex::kAuto && nFound[1] > nFound[0]);

    // ---- Per hit: the charge of each SimTrack, own and secondary, and the labels (as SimTrack indices).
    struct HitInfo {
      int layer, idx, type;
      unsigned int detId;
      const TrackingRecHit* rh;
      double qTot = 0, qNoise = 0, qUnm = 0, qOrph = 0, qOwnAll = 0, qSecAll = 0;
      std::vector<KeyCharge> keys;
      int labK = -1, ancK = -1;  // the labels: their entries in keys
    };
    std::vector<HitInfo> hits;
    long long nHits = 0, nLinks = 0, nLinksUnm = 0, nSecOnly = 0, nBelowShare = 0, nNeutralOwn = 0;
    double sumOwn = 0, sumSec = 0, sumOrph = 0, sumUnm = 0, sumNoise = 0, sumTot = 0;
    std::vector<int> nAncSt(nSt, 0);  // hits each SimTrack labels by the ancestor rule

    auto doHits = [&](const MkFitClusterIndexToHit& clusterIndex,
                      const std::vector<int>& layerOf,
                      const SimLinks& links,
                      int type) {
      unsigned int curDet = 0;
      std::vector<std::pair<unsigned int, const PixelDigiSimLink*>> chanLinks;
      std::vector<ShCharge> shq;
      for (int i = 0; i < (int)layerOf.size(); ++i) {
        const int ilay = layerOf[i];
        if (ilay < 0)
          continue;
        ++nHits;
        const TrackingRecHit* rh = clusterIndex.hits()[i];
        const OmniClusterRef& cl = static_cast<const BaseTrackerRecHit*>(rh)->firstClusterRef();
        const unsigned int detId = rh->geographicalId().rawId();
        if (detId != curDet) {
          curDet = detId;
          chanLinks.clear();
          auto ds = links.find(detId);
          if (ds != links.end())
            for (const PixelDigiSimLink& l : ds->data)
              chanLinks.emplace_back(l.channel(), &l);
          std::sort(chanLinks.begin(), chanLinks.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        }
        std::vector<std::pair<unsigned int, double>> chans;
        if (type == kHitPixel) {
          const SiPixelCluster& c = cl.pixelCluster();
          for (int ip = 0; ip < c.size(); ++ip) {
            const SiPixelCluster::Pixel& px = c.pixel(ip);
            chans.emplace_back(PixelChannelIdentifier::pixelToChannel(px.x, px.y), px.adc);
          }
        } else {
          const Phase2TrackerCluster1D& c = cl.phase2OTCluster();
          for (unsigned int is = 0; is < c.size(); ++is)
            chans.emplace_back(Phase2TrackerDigi::pixelToChannel(c.firstRow() + is, c.column()), 1.);
        }

        HitInfo hi;
        hi.layer = ilay;
        hi.idx = i;
        hi.type = type;
        hi.detId = detId;
        hi.rh = rh;
        bool anyLink = false;
        shq.clear();
        for (const auto& [chan, w] : chans) {
          hi.qTot += w;
          auto range = std::equal_range(chanLinks.begin(),
                                        chanLinks.end(),
                                        std::make_pair(chan, (const PixelDigiSimLink*)nullptr),
                                        [](const auto& a, const auto& b) { return a.first < b.first; });
          double sumf = 0;
          for (auto it = range.first; it != range.second; ++it)
            sumf += it->second->fraction();
          const double scale = sumf > 1 ? 1. / sumf : 1.;
          for (auto it = range.first; it != range.second; ++it) {
            const PixelDigiSimLink& l = *it->second;
            anyLink = true;
            ++nLinks;
            const double q = w * l.fraction() * scale;
            const ShRef r = resolve(detId, l, perModule);
            if (!r.valid()) {
              ++nLinksUnm;
              hi.qUnm += q;
              continue;
            }
            auto f =
                std::find_if(shq.begin(), shq.end(), [&](const ShCharge& s) { return s.r.c == r.c && s.r.j == r.j; });
            if (f == shq.end())
              shq.push_back(ShCharge{r, q});
            else
              f->q += q;
          }
          hi.qNoise += w * std::max(0., 1. - sumf);
        }
        if (!anyLink || hi.qTot <= 0)
          continue;

        // Per SimTrack: its own and secondary charge, its sim hits with the most charge.
        for (const ShCharge& s : shq) {
          const int st = shSt[s.r.c][s.r.j];
          const std::int8_t k = shCls[s.r.c][s.r.j];
          if (st < 0) {
            hi.qOrph += s.q;
            continue;
          }
          auto f = std::find_if(hi.keys.begin(), hi.keys.end(), [&](const KeyCharge& x) { return x.st == st; });
          if (f == hi.keys.end()) {
            hi.keys.emplace_back();
            f = hi.keys.end() - 1;
            f->st = st;
          }
          if (k == kOwn) {
            f->own += s.q;
            hi.qOwnAll += s.q;
            if (s.q > f->bestOwnQ) {
              f->bestOwnQ = s.q;
              f->bestOwn = s.r;
            }
          } else {
            f->sec += s.q;
            hi.qSecAll += s.q;
            if (s.q > f->bestSecQ) {
              f->bestSecQ = s.q;
              f->bestSec = s.r;
            }
          }
          if (s.q > f->bestAnyQ) {
            f->bestAnyQ = s.q;
            f->bestAny = s.r;
          }
        }
        // The keys in the crossing's order, so that a tie goes to the SimTrack first there.
        std::sort(hi.keys.begin(), hi.keys.end(), [](const KeyCharge& a, const KeyCharge& b) { return a.st < b.st; });
        // "sim": the charged SimTrack with the most own charge. The ancestor rule: any SimTrack with the most charge,
        // own and secondary (neutral ones too, as TrackingParticles take them). A strict > keeps the first on a tie.
        double bestOwn = 0, bestOwnNeutral = 0, bestAll = 0;
        for (int k = 0; k < (int)hi.keys.size(); ++k) {
          const KeyCharge& x = hi.keys[k];
          if (x.own > 0) {
            if (stCharge[x.st] == 0) {
              bestOwnNeutral = std::max(bestOwnNeutral, x.own);
            } else if (x.own > bestOwn) {
              bestOwn = x.own;
              hi.labK = k;
            }
          }
          if (x.own + x.sec > bestAll) {
            bestAll = x.own + x.sec;
            hi.ancK = k;
          }
        }
        nNeutralOwn += bestOwnNeutral > bestOwn;
        if (hi.labK >= 0 && bestOwn < minOwnShare_ * hi.qTot) {
          hi.labK = -1;
          ++nBelowShare;
        }
        if (hi.qOwnAll <= 0 && hi.qSecAll + hi.qOrph > 0)
          ++nSecOnly;
        if (hi.ancK >= 0)
          ++nAncSt[hi.keys[hi.ancK].st];
        sumOwn += hi.qOwnAll;
        sumSec += hi.qSecAll;
        sumOrph += hi.qOrph;
        sumUnm += hi.qUnm;
        sumNoise += hi.qNoise;
        sumTot += hi.qTot;
        hits.push_back(std::move(hi));
      }
    };
    doHits(pixelIndex, pixelLayer, pixelLinks, kHitPixel);
    doHits(stripIndex, stripLayer, phase2Links, kHitPhase2OT);

    // ---- The sim tracks: charged SimTracks with own sim hits in the tracker, and every SimTrack the ancestor rule
    // labels a hit with; in the crossing's order.
    std::vector<int> simOfSt(nSt, -1), stOfSim;
    int nSimsOwn = 0, nSimsAncOnly = 0;
    for (int st = 0; st < nSt; ++st) {
      const bool own = nOwnSh[st] > 0 && stCharge[st] != 0;
      const bool anc = withAnc && nAncSt[st] > 0;
      if (!own && !anc)
        continue;
      simOfSt[st] = stOfSim.size();
      stOfSim.push_back(st);
      nSimsOwn += own;
      nSimsAncOnly += !own;
    }
    const int nSims = stOfSim.size();

    // ---- The TrackingParticles of the other sources: per SimTrack its TrackingParticle there.
    std::vector<std::vector<int>> tpOfSt(tpTokens_.size(), std::vector<int>(nSt, -1));
    std::vector<const TrackingParticleRefVector*> tpVecs;
    long long nTpConflicts = 0, nTpG4Missing = 0;
    for (std::size_t s = 0; s < tpTokens_.size(); ++s) {
      const TrackingParticleRefVector& tps = iEvent.get(tpTokens_[s]);
      tpVecs.push_back(&tps);
      for (int i = 0; i < (int)tps.size(); ++i)
        for (const SimTrack& g4 : tps[i]->g4Tracks()) {
          auto f = stOfKey.find(idKey(g4.eventId().rawId(), g4.trackId()));
          if (f == stOfKey.end()) {
            ++nTpG4Missing;
            continue;
          }
          int& t = tpOfSt[s][f->second];
          if (t >= 0 && t != i)
            ++nTpConflicts;
          else
            t = i;
        }
    }

    // ---- chi2 as the source "ntuple": the layers of the first track of the SimTrack's TrackingParticle.
    auto layerOfHit = [&](int type, int idx) -> int {
      const std::vector<int>& lo = type == kHitPixel ? pixelLayer : stripLayer;
      return idx >= 0 && idx < (int)lo.size() ? lo[idx] : -1;
    };
    std::vector<int> tpNlay;
    if (chi2Source_ >= 0) {
      const TrackingParticleRefVector& tps = *tpVecs[chi2Source_];
      tpNlay.assign(tps.size(), 0);
      const edm::Handle<edm::View<reco::Track>> tracksHandle = iEvent.getHandle(tracksToken_);
      edm::RefToBaseVector<reco::Track> trackRefs;
      for (std::size_t i = 0; i < tracksHandle->size(); ++i)
        trackRefs.push_back(tracksHandle->refAt(i));
      const reco::SimToRecoCollection simRecColl = iEvent.get(associatorToken_).associateSimToReco(trackRefs, tps);
      for (int s = 0; s < (int)tps.size(); ++s) {
        auto found = simRecColl.find(tps[s]);
        if (found == simRecColl.end() || found->val.empty())
          continue;
        const reco::Track& trk = *found->val.front().first;
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
        }
        for (int c : hitlay)
          tpNlay[s] += c > 0;
      }
    }

    // ---- The tables of both sources.
    const edm::Handle<TrajectorySeedCollection> seedsHandle = iEvent.getHandle(seedsToken_);
    const TrajectorySeedCollection& seeds = *seedsHandle;
    TraceTruthTables T[2];
    const int nSrc = withAnc ? 2 : 1;
    for (int s = 0; s < nSrc; ++s) {
      TraceTruthTables& t = T[s];
      t.source = s == 0 ? source_ : sourceAncestor_;
      t.run = iEvent.id().run();
      t.lumi = iEvent.id().luminosityBlock();
      t.event = iEvent.id().event();
      t.module = module_;
      t.layers.resize(trackerInfo.n_layers());
      t.simHitStates.resize(trackerInfo.n_layers());
      for (int l = 0; l < (int)t.layers.size(); ++l)
        t.layers[l].layer = t.simHitStates[l].layer = l;
    }
    // n_hits over every registered hit, as the other sources count it.
    for (const auto* lo : {&pixelLayer, &stripLayer})
      for (int i = 0; i < (int)lo->size(); ++i)
        if ((*lo)[i] >= 0)
          for (int s = 0; s < nSrc; ++s) {
            auto& L = T[s].layers.at((*lo)[i]);
            L.n_hits = std::max(L.n_hits, i + 1);
            auto& S = T[s].simHitStates.at((*lo)[i]);
            S.n_hits = std::max(S.n_hits, i + 1);
          }

    // The sim tracks' states (shared), then each source's hits.
    std::vector<mkfit::Track> simTrackBase;
    simTrackBase.reserve(nSims);
    int nProd[3] = {0, 0, 0};
    for (int k = 0; k < nSims; ++k) {
      const SimTrack& st = simTracks[stOfSim[k]];
      const SimVertex* sv = stVertex[stOfSim[k]];
      const float vx = sv ? sv->position().x() : 0.f, vy = sv ? sv->position().y() : 0.f,
                  vz = sv ? sv->position().z() : 0.f;
      const float px = st.momentum().px(), py = st.momentum().py(), pz = st.momentum().pz();
      mkfit::SMatrixSym66 err;
      err.At(0, 0) = vx * vx;
      err.At(1, 1) = vy * vy;
      err.At(2, 2) = vz * vz;
      err.At(3, 3) = px * px;
      err.At(4, 4) = py * py;
      err.At(5, 5) = pz * pz;
      mkfit::TrackState state(stCharge[stOfSim[k]], mkfit::SVector3(vx, vy, vz), mkfit::SVector3(px, py, pz), err);
      state.convertFromCartesianToCCS();
      int nlay = 0;
      if (chi2Source_ >= 0) {
        const int tp = tpOfSt[chi2Source_][stOfSim[k]];
        if (tp >= 0)
          nlay = tpNlay[tp];
      }
      mkfit::Track& t = simTrackBase.emplace_back(state, float(nlay), k, 0, nullptr);
      const EncodedEventId eid = st.eventId();
      if (eid.bunchCrossing() == 0) {
        t.setProdType(eid.event() == 0 ? mkfit::Track::ProdType::Signal : mkfit::Track::ProdType::InTimePU);
        ++nProd[eid.event() == 0 ? 0 : 1];
      } else {
        t.setProdType(mkfit::Track::ProdType::OutOfTimePU);
        ++nProd[2];
      }
    }

    auto shOf = [&](const ShRef& r) -> const PSimHit& { return (*vecs[r.c])[r.j]; };
    int nSimsHits = 0;
    long long nLabels[2] = {0, 0};
    std::vector<std::vector<int>> hitsOfSimPerSrc[2];  // per source, per sim track: its hits (index into hits)
    for (int s = 0; s < nSrc; ++s) {
      TraceTruthTables& t = T[s];
      auto& hos = hitsOfSimPerSrc[s];
      hos.assign(nSims, {});
      for (int h = 0; h < (int)hits.size(); ++h) {
        const HitInfo& hi = hits[h];
        const int k = s == 0 ? hi.labK : hi.ancK;
        if (k < 0)
          continue;
        const KeyCharge& kc = hi.keys[k];
        const int sim = simOfSt[kc.st];
        if (sim < 0)
          throw cms::Exception("LogicError") << "MkFitTraceTruthSim: a label without a sim track";
        ++nLabels[s];
        auto& L = t.layers.at(hi.layer);
        L.hit.push_back(hi.idx);
        L.sim.push_back(sim);
        hos[sim].push_back(h);
        // The sim hit behind the label.
        const ShRef r = kc.bestOwn.valid() ? kc.bestOwn : kc.bestAny;
        const PSimHit& sh = shOf(r);
        const auto* det = tracker.idToDetUnit(DetId(sh.detUnitId()));
        if (!det)
          throw cms::Exception("LogicError") << "MkFitTraceTruthSim: no det unit for DetId " << sh.detUnitId();
        const auto pos = det->surface().toGlobal(sh.localPosition());
        const auto mom = det->surface().toGlobal(sh.momentumAtEntry());
        auto& S = t.simHitStates.at(hi.layer);
        S.hit.push_back(hi.idx);
        S.sh_idx.push_back(vecOffset[r.c] + r.j);
        S.x.push_back(pos.x());
        S.y.push_back(pos.y());
        S.z.push_back(pos.z());
        S.px.push_back(mom.x());
        S.py.push_back(mom.y());
        S.pz.push_back(mom.z());
        mkfit::TrackState ts(stCharge[kc.st],
                             mkfit::SVector3(pos.x(), pos.y(), pos.z()),
                             mkfit::SVector3(mom.x(), mom.y(), mom.z()),
                             {});
        ts.convertFromCartesianToCCS();
        S.ts_q.push_back(ts.charge);
        S.ts_par.push_back({{ts.parameters[0],
                             ts.parameters[1],
                             ts.parameters[2],
                             ts.parameters[3],
                             ts.parameters[4],
                             ts.parameters[5]}});
        S.linked.push_back(1);
      }
      // The sim tracks with their hits, in the time of flight of the sim hit behind each label.
      t.simTracks = simTrackBase;
      for (int k = 0; k < nSims; ++k) {
        auto& hl = hos[k];
        auto tof = [&](int h) {
          const HitInfo& hi = hits[h];
          const KeyCharge& kc = hi.keys[s == 0 ? hi.labK : hi.ancK];
          return shOf(kc.bestOwn.valid() ? kc.bestOwn : kc.bestAny).timeOfFlight();
        };
        std::stable_sort(hl.begin(), hl.end(), [&](int a, int b) {
          const float ta = tof(a), tb = tof(b);
          if (ta != tb)
            return ta < tb;
          if (hits[a].layer != hits[b].layer)
            return hits[a].layer < hits[b].layer;
          return hits[a].idx < hits[b].idx;
        });
        for (int h : hl)
          t.simTracks[k].addHitIdx(hits[h].idx, hits[h].layer, 0.f);
        if (s == 0)
          nSimsHits += !hl.empty();
      }
    }

    // ---- Seeds: their hits, and their labels from the hits' labels.
    std::map<std::pair<int, int>, int> hitOfKey;  // (layer, idx) -> index into hits
    for (int h = 0; h < (int)hits.size(); ++h)
      hitOfKey.emplace(std::make_pair(hits[h].layer, hits[h].idx), h);
    int nSeedsSim = 0, nSeedsBest = 0;
    for (int s = 0; s < nSrc; ++s) {
      TraceTruthTables& t = T[s];
      t.seeds.resize(seeds.size());
      for (int i = 0; i < (int)seeds.size(); ++i) {
        auto& ss = t.seeds[i];
        ss.seed = i;
        std::vector<int> labs;
        int nValid = 0;
        for (const TrackingRecHit& hit : seeds[i].recHits()) {
          if (!hit.isValid())
            continue;
          ++nValid;
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
          int lab = -1;
          auto f = hitOfKey.find({lay, idx});
          if (f != hitOfKey.end()) {
            const HitInfo& hi = hits[f->second];
            const int k = s == 0 ? hi.labK : hi.ancK;
            if (k >= 0)
              lab = simOfSt[hi.keys[k].st];
          }
          labs.push_back(lab);
        }
        // Count per label, with the position of its innermost hit.
        std::map<int, std::pair<int, int>> cnt;  // label -> (count, innermost hit)
        for (int h = 0; h < (int)labs.size(); ++h)
          if (labs[h] >= 0) {
            auto& c = cnt.emplace(labs[h], std::make_pair(0, h)).first->second;
            ++c.first;
          }
        // quickTrackAssociatorByHits: the largest shared fraction above 0.75, and 3 of 3 for a 3-hit seed.
        int bestQ = -1, bestQn = 0;
        for (const auto& [lab, c] : cnt)
          if (c.first > bestQn) {
            bestQn = c.first;
            bestQ = lab;
          }
        if (bestQ >= 0 && nValid > 0) {
          const float frac = float(bestQn) / nValid;
          if (frac > 0.75f && !(nValid == 3 && bestQn < 3)) {
            ss.sim = bestQ;
            ss.sim_frac = frac;
          }
        }
        // see_bestSimTrkIdx: at least 3, the innermost on a tie.
        int bestCount = 2, bestInner = std::numeric_limits<int>::max();
        for (const auto& [lab, c] : cnt)
          if (c.first > bestCount || (c.first == bestCount && c.second < bestInner)) {
            bestCount = c.first;
            bestInner = c.second;
            ss.best_sim = lab;
            ss.best_frac = nValid > 0 ? float(c.first) / nValid : 0.f;
          }
        if (s == 0) {
          nSeedsSim += ss.sim >= 0;
          nSeedsBest += ss.best_sim >= 0;
        }
      }
    }

    // ---- SimTrackMap (both sources).
    std::vector<int> nTpOfSims(tpTokens_.size(), 0);
    for (int s = 0; s < nSrc; ++s) {
      auto& m = T[s].simTrackMap;
      m.present = true;
      m.tp_names = tpNames_;
      m.tp_index.assign(tpNames_.size(), {});
      for (int k = 0; k < nSims; ++k) {
        const int st = stOfSim[k];
        const SimTrack& t = simTracks[st];
        m.g4_id.push_back(t.trackId());
        m.bx.push_back(t.eventId().bunchCrossing());
        m.ev.push_back(t.eventId().event());
        m.pdg.push_back(t.type());
        m.proc.push_back(stProc[st]);
        m.parent_g4.push_back(stVertex[st] ? stVertex[st]->parentIndex() : -1);
        m.gen.push_back(t.genpartIndex());
        m.n_own_simhits.push_back(nOwnSh[st]);
        for (std::size_t n = 0; n < tpNames_.size(); ++n) {
          m.tp_index[n].push_back(tpOfSt[n][st]);
          if (s == 0)
            nTpOfSims[n] += tpOfSt[n][st] >= 0;
        }
      }
    }

    int pass = writer->writeTruth(T[0]);
    if (withAnc)
      writer->writeTruth(T[1]);

    // ---- The CSV files for the comparison of the sources.
    if (dumping) {
      std::ostringstream key;
      key << T[0].run << "_" << T[0].lumi << "_" << T[0].event << "_p" << pass;
      const std::string pfx = dumpPrefix_ + key.str();
      std::ofstream fh(pfx + "-hits.csv"), fk(pfx + "-hitkeys.csv"), fs(pfx + "-simtracks.csv"), ft(pfx + "-tps.csv");
      if (!fh || !fk || !fs || !ft)
        throw cms::Exception("Configuration") << "MkFitTraceTruthSim: cannot write the dump " << pfx << "-*.csv";
      char buf[1024];
      fh << "layer,idx,detid,pix,q_tot,q_noise,q_unm,lx,ly,gx,gy,gz,n_keys,sim,simanc,f_lab_own,f_own,f_sec,f_orph\n";
      std::string tpCols;
      for (const std::string& n : tpNames_)
        tpCols += ",tp_" + n;
      fk << "layer,idx,st,bx,ev,tid,pdg,charge,sim" << tpCols
         << ",q_own,q_sec,own_in,own_lx,own_ly,own_gx,own_gy,own_gz,own_p,own_tof,n_own_mod,sec_pdg,sec_proc,sec_p\n";
      for (const HitInfo& hi : hits) {
        const auto* det = tracker.idToDetUnit(DetId(hi.detId));
        const LocalPoint lp = hi.rh->localPosition();
        const GlobalPoint gp = det->surface().toGlobal(lp);
        const int lab = hi.labK >= 0 ? simOfSt[hi.keys[hi.labK].st] : -1;
        const int anc = hi.ancK >= 0 ? simOfSt[hi.keys[hi.ancK].st] : -1;
        std::snprintf(buf,
                      sizeof(buf),
                      "%d,%d,%u,%d,%.6g,%.6g,%.6g,%.6f,%.6f,%.6f,%.6f,%.6f,%zu,%d,%d,%.5f,%.5f,%.5f,%.5f\n",
                      hi.layer,
                      hi.idx,
                      hi.detId,
                      hi.type == kHitPixel,
                      hi.qTot,
                      hi.qNoise,
                      hi.qUnm,
                      lp.x(),
                      lp.y(),
                      gp.x(),
                      gp.y(),
                      gp.z(),
                      hi.keys.size(),
                      lab,
                      anc,
                      hi.labK >= 0 ? hi.keys[hi.labK].own / hi.qTot : 0.,
                      hi.qOwnAll / hi.qTot,
                      hi.qSecAll / hi.qTot,
                      hi.qOrph / hi.qTot);
        fh << buf;
        for (const KeyCharge& kc : hi.keys) {
          const SimTrack& t = simTracks[kc.st];
          // Its own sim hit with the most charge in the cluster, else its own sim hit on the module nearest the hit.
          ShRef r = kc.bestOwn;
          int ownIn = r.valid() ? 1 : -1;
          int nOwnMod = 0;
          auto od = ownOnDet.find(detKey(hi.detId, kc.st));
          if (od != ownOnDet.end()) {
            nOwnMod = od->second.size();
            if (!r.valid()) {
              double best = std::numeric_limits<double>::max();
              for (const ShRef& o : od->second) {
                const Local3DPoint p = shOf(o).localPosition();
                const double d = std::hypot(p.x() - lp.x(), p.y() - lp.y());
                if (d < best) {
                  best = d;
                  r = o;
                }
              }
              ownIn = 0;
            }
          }
          std::ostringstream os;
          os << hi.layer << "," << hi.idx << "," << kc.st << "," << t.eventId().bunchCrossing() << ","
             << t.eventId().event() << "," << t.trackId() << "," << t.type() << "," << stCharge[kc.st] << ","
             << simOfSt[kc.st];
          for (std::size_t n = 0; n < tpNames_.size(); ++n)
            os << "," << tpOfSt[n][kc.st];
          if (r.valid()) {
            const PSimHit& sh = shOf(r);
            const GlobalPoint sg = det->surface().toGlobal(sh.localPosition());
            std::snprintf(buf,
                          sizeof(buf),
                          ",%.6g,%.6g,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6g,%.6g,%d",
                          kc.own,
                          kc.sec,
                          ownIn,
                          sh.localPosition().x(),
                          sh.localPosition().y(),
                          sg.x(),
                          sg.y(),
                          sg.z(),
                          sh.pabs(),
                          sh.timeOfFlight(),
                          nOwnMod);
          } else {
            std::snprintf(buf, sizeof(buf), ",%.6g,%.6g,-1,0,0,0,0,0,-1,0,%d", kc.own, kc.sec, nOwnMod);
          }
          os << buf;
          if (kc.bestSec.valid()) {
            const PSimHit& sh = shOf(kc.bestSec);
            std::snprintf(buf, sizeof(buf), ",%d,%d,%.6g\n", sh.particleType(), (int)sh.processType(), sh.pabs());
          } else {
            std::snprintf(buf, sizeof(buf), ",0,0,-1\n");
          }
          os << buf;
          fk << os.str();
        }
      }
      fs << "sim,st,bx,ev,tid,pdg,charge,proc,parent,gen,vx,vy,vz,px,py,pz" << tpCols << ",n_own_sh,n_hits_" << source_
         << (withAnc ? ",n_hits_" + sourceAncestor_ : std::string()) << ",chi2\n";
      for (int k = 0; k < nSims; ++k) {
        const int st = stOfSim[k];
        const SimTrack& t = simTracks[st];
        const SimVertex* sv = stVertex[st];
        std::ostringstream os;
        std::snprintf(buf,
                      sizeof(buf),
                      "%d,%d,%d,%d,%u,%d,%d,%d,%d,%d,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g",
                      k,
                      st,
                      t.eventId().bunchCrossing(),
                      t.eventId().event(),
                      t.trackId(),
                      t.type(),
                      stCharge[st],
                      stProc[st],
                      sv ? sv->parentIndex() : -1,
                      t.genpartIndex(),
                      sv ? sv->position().x() : 0.,
                      sv ? sv->position().y() : 0.,
                      sv ? sv->position().z() : 0.,
                      t.momentum().px(),
                      t.momentum().py(),
                      t.momentum().pz());
        os << buf;
        for (std::size_t n = 0; n < tpNames_.size(); ++n)
          os << "," << tpOfSt[n][st];
        os << "," << nOwnSh[st] << "," << hitsOfSimPerSrc[0][k].size();
        if (withAnc)
          os << "," << hitsOfSimPerSrc[1][k].size();
        os << "," << simTrackBase[k].chi2() << "\n";
        fs << os.str();
      }
      ft << "source,idx,bx,ev,pdg,charge,vx,vy,vz,px,py,pz,g4_st\n";
      for (std::size_t n = 0; n < tpVecs.size(); ++n) {
        const TrackingParticleRefVector& tps = *tpVecs[n];
        for (int i = 0; i < (int)tps.size(); ++i) {
          const TrackingParticle& tp = *tps[i];
          const auto& v = tp.parentVertex()->position();
          std::snprintf(buf,
                        sizeof(buf),
                        "%s,%d,%d,%d,%d,%d,%.7g,%.7g,%.7g,%.7g,%.7g,%.7g,",
                        tpNames_[n].c_str(),
                        i,
                        tp.eventId().bunchCrossing(),
                        tp.eventId().event(),
                        tp.pdgId(),
                        (int)tp.charge(),
                        v.x(),
                        v.y(),
                        v.z(),
                        tp.px(),
                        tp.py(),
                        tp.pz());
          std::ostringstream os;
          os << buf;
          bool first = true;
          for (const SimTrack& g4 : tp.g4Tracks()) {
            auto f = stOfKey.find(idKey(g4.eventId().rawId(), g4.trackId()));
            os << (first ? "" : ";") << (f == stOfKey.end() ? -1 : f->second);
            first = false;
          }
          os << "\n";
          ft << os.str();
        }
      }
    }

    std::ostringstream tpCounts;
    for (std::size_t n = 0; n < tpNames_.size(); ++n)
      tpCounts << (n ? "," : "") << tpNames_[n] << ":" << nTpOfSims[n];
    const double tot = sumTot > 0 ? sumTot : 1;
    char buf[2048];
    std::snprintf(buf,
                  sizeof(buf),
                  "MKFITTRUTHSIM %u:%u:%llu %s source=%s ancestor=%s simtracks=%d simvertices=%zu vert_bad=%lld "
                  "simhits=%lld simhits_own=%lld simhits_sec=%lld simhits_orph=%lld sims=%d sims_own=%d "
                  "sims_anc_only=%d sims_prod=%d/%d/%d sims_hits=%d sims_tp=%s hits=%lld hits_linked=%zu labels=%lld "
                  "labels_anc=%lld sec_only=%lld below_share=%lld neutral_own=%lld q_own=%.4f q_sec=%.4f "
                  "q_orph=%.4f q_unm=%.4f q_noise=%.4f links=%lld links_unmatched=%lld link_index=%s "
                  "event_links=%lld found_global=%lld found_perModule=%lld seeds=%zu seeds_sim=%d seeds_best=%d "
                  "tp_conflicts=%lld tp_g4_missing=%lld pass=%d ms=%.1f",
                  T[0].run,
                  T[0].lumi,
                  T[0].event,
                  module_.c_str(),
                  source_.c_str(),
                  withAnc ? sourceAncestor_.c_str() : "-",
                  nSt,
                  simVertices.size(),
                  nVertBad,
                  nSimHits,
                  nShCls[kOwn],
                  nShCls[kSec],
                  nShCls[kOrph],
                  nSims,
                  nSimsOwn,
                  nSimsAncOnly,
                  nProd[0],
                  nProd[1],
                  nProd[2],
                  nSimsHits,
                  tpCounts.str().c_str(),
                  nHits,
                  hits.size(),
                  nLabels[0],
                  nLabels[1],
                  nSecOnly,
                  nBelowShare,
                  nNeutralOwn,
                  sumOwn / tot,
                  sumSec / tot,
                  sumOrph / tot,
                  sumUnm / tot,
                  sumNoise / tot,
                  nLinks,
                  nLinksUnm,
                  perModule ? "perModule" : "global",
                  nEvLinks,
                  nFound[0],
                  nFound[1],
                  seeds.size(),
                  nSeedsSim,
                  nSeedsBest,
                  nTpConflicts,
                  nTpG4Missing,
                  pass,
                  std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
    std::cout << buf << std::endl;
  }

}  // namespace mkfit_fwff

using mkfit_fwff::MkFitTraceTruthSim;
DEFINE_FWK_MODULE(MkFitTraceTruthSim);
