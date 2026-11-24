// Package:    MkFitFWFF/Trace
// Class:      MkFitTraceHitOrigin
//
// The origin of each mkFit hit's charge in one MkFitProducer's traced event, from the simulation's own collections:
// which part of a hit's cluster charge a SimTrack's own particle deposited, which part a secondary that is not a
// SimTrack deposited (a delta ray, the products of an interaction), and which part cannot be told (pileup). It gives
// the table to the MkFitTraceWriter service, which writes it into the trace file next to the records
// (TraceRNTupleWriter.h: TruthHitOrigin). The loader of the standalone joins it to the records, next to the truth of
// MkFitTraceTruth.
//
// Why: Geant4 does not keep every particle as a SimTrack. The tracker's sensitive detector gives a sim hit of a
// particle that is not kept the id of its parent (SimG4CMS/Tracker TkAccumulatingSensitiveDetector::createHit()), and
// its own particle type and creating process (PSimHit::particleType(), processType()). A delta ray that leaves a
// cluster off the track therefore carries the track's id, and a truth that goes by the id alone (TrackingParticles,
// ClusterTPAssociation) gives the cluster the track's label.
//
// - A hit is (layer, index): the index is the cluster's index in its cluster collection (MkFitHitWrapper), as the
//   records keep it. Its charge: per pixel of the cluster the pixel's ADC; per strip of an outer-tracker cluster 1.
// - The phase-2 digitizer (SimTracker/SiPhase2Digitizer) writes one digi-sim link per channel and contributing sim
//   hit (one per charge deposit of it, summed here): the sim hit's track id and event id, its share of the channel's
//   signal (fraction; the total includes the noise, so on a channel the shares can add up to more than 1, and are then
//   scaled to 1) and an index of the sim hit in the collection it was read from (CFposition, with the tof bin). The
//   signal event is accumulated before any pileup, so the signal's indices start at 0 in each collection. What the
//   index counts depends on the release that digitized the sample:
//   global     the sim hit's index in the collection (the per-hit accumulateSimHits() since b901d3215d4, 2026-08-21,
//              first in CMSSW_20_1_0_pre3);
//   perModule  the index of the first sim hit on the module plus the sim hit's rank among the module's sim hits, in
//              collection order (Phase2TrackerDigitizerAlgorithm::accumulateSimHits() up to then counts only the
//              module's hits, e.g. CMSSW_20_0_0_pre1).
//   linkIndex "auto" decodes every signal link of the event both ways and takes the one under which more links name a
//   sim hit on their own module with their own track id (the log gives both counts); "global" and "perModule" force
//   one. A signal link whose sim hit is not found so is unmatched.
// - A link with event id 0 in bunch crossing 0 is signal; the sim hit it names is classified:
//   own        the SimTrack with the sim hit's track id has the sim hit's particle type, and its creating process
//              (the process type of its SimVertex) is the sim hit's process type;
//   secondary  the SimTrack has another particle type or another process: a particle that is not a SimTrack
//              deposited the charge under the id of its parent; or no SimTrack has the sim hit's track id (its parent
//              was not kept either), counted apart as "orphan" in the log.
//   Blind spot of the rule: a secondary of the same type and process as its kept parent (a delta ray of an electron
//   that was itself made by ionisation) counts as own.
// - Any other link is pileup: the job has no pileup sim hits (a RAW file keeps only the signal's), so its charge
//   cannot be classified.
// - Per hit: f_own, f_sec (secondary and orphan), f_pu, f_unm, the fractions of the hit's charge, each link
//   contributing its share times the channel's charge; the rest is charge without a link (noise). cls: the largest of
//   own, secondary and unknown (pileup + unmatched), own winning a tie, then secondary. sim: the sim track (the index
//   into the trace file's SimTracks, the TrackingParticles of trackingParticles in their order, as MkFitTraceTruth
//   writes them) of the signal SimTrack with most of the hit's charge, own and secondary together; -1 without
//   trackingParticles. sec_pdg, sec_proc: the secondary sim hit with most charge in the hit.
// - Listed are the hits with at least one link.
//
// Prints one line per event:
//   MKFITHITORIGIN <run>:<lumi>:<event> <module> hits=<registered mkFit hits> listed=<with a link>
//     cls_own=<n> cls_sec=<n> cls_unk=<n> cls_none=<n> orph_dom=<hits where orphan charge is the largest class share>
//     q_own=.. q_sec=.. q_orph=.. q_pu=.. q_unm=.. q_noise=.. (shares of the listed hits' total charge)
//     links_sig=<on the hits> links_unmatched=<n> links_pu=<n> link_index=<global|perModule>
//     event_links_sig=<all signal links of the event> found_global=<n> found_perModule=<n> simtracks=<n>
//     simhits=<signal sim hits read> simhits_own=<n> simhits_sec=<n> simhits_orph=<n> pass=<the table's pass>
//     ms=<time>
// one line with the most frequent (particle type, process) pairs among the secondary sim hits on hits:
//   MKFITHITORIGIN-SECTYPES <key> <pdg>/<proc>:<sim hits>:<charge share of all secondary charge> ...
// and one line on the unmatched signal links on hits, pixel and outer tracker apart:
//   MKFITHITORIGIN-UNMATCHED <key> pixel: no_simhit_in_tof_bin=.. out_of_range=.. other_module=.. other_track=..
//     track_has_simhit_here=..; outer: ...
// With printSecondary > 0, for every printStride-th hit whose class is secondary, up to printSecondary of them, the
// hit and each sim hit behind its signal charge (MKFITHITORIGIN-HIT, -SIMHIT lines) and every sim hit of the same
// SimTracks on the same module (-ONMODULE lines), with their local positions against the cluster's. With
// printUnmatched > 0, up to that many unmatched signal links with the sim hits around their index (-UNMLINK).

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
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PixelChannelIdentifier.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "RecoTracker/MkFit/interface/MkFitClusterIndexToHit.h"
#include "RecoTracker/MkFit/interface/MkFitGeometry.h"
#include "RecoTracker/MkFit/interface/MkFitOutputWrapper.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"

#include "MkFitFWFF/Trace/plugins/MkFitTraceWriter.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mkfit_fwff {

  class MkFitTraceHitOrigin : public edm::global::EDAnalyzer<> {
  public:
    explicit MkFitTraceHitOrigin(const edm::ParameterSet& iConfig);

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void analyze(edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

    using SimLinks = edm::DetSetVector<PixelDigiSimLink>;
    enum class LinkIndex { kAuto, kGlobal, kPerModule };

    const std::string module_;
    const edm::EDGetTokenT<MkFitOutputWrapper> mkFitOutputToken_;
    const edm::EDGetTokenT<MkFitClusterIndexToHit> pixelClusterIndexToken_;
    const edm::EDGetTokenT<MkFitClusterIndexToHit> stripClusterIndexToken_;
    const edm::EDGetTokenT<std::vector<int>> pixelLayerToken_;
    const edm::EDGetTokenT<std::vector<int>> stripLayerToken_;
    const edm::EDGetTokenT<SimLinks> pixelSimLinkToken_;
    const edm::EDGetTokenT<SimLinks> phase2OTSimLinkToken_;
    std::vector<edm::EDGetTokenT<std::vector<PSimHit>>> simHitTokens_;
    std::vector<std::string> simHitNames_;
    std::vector<unsigned int> simHitTof_;  // PixelDigiSimLink::LowTof or HighTof, from the instance name
    const edm::EDGetTokenT<edm::SimTrackContainer> simTracksToken_;
    const edm::EDGetTokenT<edm::SimVertexContainer> simVerticesToken_;
    const bool useTrackingParticles_;
    edm::EDGetTokenT<TrackingParticleRefVector> trackingParticlesToken_;
    LinkIndex linkIndex_;
    const int printSecondary_;
    const int printStride_;
    const int printUnmatched_;
    const edm::ESGetToken<MkFitGeometry, TrackerRecoGeometryRecord> mkFitGeomToken_;
  };

  namespace {
    std::uint64_t simTrackKey(unsigned int trackId, EncodedEventId eid) {
      return (std::uint64_t(eid.rawId()) << 32) | trackId;
    }
    std::uint64_t detKey(unsigned int detId, unsigned int x) { return (std::uint64_t(detId) << 32) | x; }

    enum SimHitClass { kShOwn = 0, kShSec = 1, kShOrph = 2 };
    const char* const kShClassName[3] = {"own", "secondary", "orphan"};

    bool isSignal(EncodedEventId eid) { return eid.event() == 0 && eid.bunchCrossing() == 0; }

    std::string fmtPoint(const Local3DPoint& p) {
      char buf[96];
      std::snprintf(buf, sizeof(buf), "(%.5f,%.5f,%.5f)", p.x(), p.y(), p.z());
      return buf;
    }
  }  // namespace

  MkFitTraceHitOrigin::MkFitTraceHitOrigin(const edm::ParameterSet& iConfig)
      : module_(iConfig.getParameter<std::string>("mkFitProducer")),
        mkFitOutputToken_(consumes(edm::InputTag(module_))),
        pixelClusterIndexToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelHits"))),
        stripClusterIndexToken_(consumes(iConfig.getParameter<edm::InputTag>("stripHits"))),
        pixelLayerToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelHits"))),
        stripLayerToken_(consumes(iConfig.getParameter<edm::InputTag>("stripHits"))),
        pixelSimLinkToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelSimLink"))),
        phase2OTSimLinkToken_(consumes(iConfig.getParameter<edm::InputTag>("phase2OTSimLink"))),
        simTracksToken_(consumes(iConfig.getParameter<edm::InputTag>("simTracks"))),
        simVerticesToken_(consumes(iConfig.getParameter<edm::InputTag>("simVertices"))),
        useTrackingParticles_(!iConfig.getParameter<edm::InputTag>("trackingParticles").label().empty()),
        printSecondary_(iConfig.getUntrackedParameter<int>("printSecondary")),
        printStride_(std::max(1, iConfig.getUntrackedParameter<int>("printStride"))),
        printUnmatched_(iConfig.getUntrackedParameter<int>("printUnmatched")),
        mkFitGeomToken_(esConsumes()) {
    for (const edm::InputTag& t : iConfig.getParameter<std::vector<edm::InputTag>>("simHits")) {
      simHitTokens_.push_back(consumes<std::vector<PSimHit>>(t));
      simHitNames_.push_back(t.instance());
      simHitTof_.push_back(t.instance().find("HighTof") != std::string::npos ? PixelDigiSimLink::HighTof
                                                                             : PixelDigiSimLink::LowTof);
    }
    if (useTrackingParticles_)
      trackingParticlesToken_ =
          consumes<TrackingParticleRefVector>(iConfig.getParameter<edm::InputTag>("trackingParticles"));
    const std::string li = iConfig.getParameter<std::string>("linkIndex");
    if (li == "auto")
      linkIndex_ = LinkIndex::kAuto;
    else if (li == "global")
      linkIndex_ = LinkIndex::kGlobal;
    else if (li == "perModule")
      linkIndex_ = LinkIndex::kPerModule;
    else
      throw cms::Exception("Configuration")
          << "MkFitTraceHitOrigin: linkIndex must be auto, global or perModule, not '" << li << "'";
  }

  void MkFitTraceHitOrigin::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("mkFitProducer", "initialStepTrackCandidatesMkFit")
        ->setComment(
            "label of the traced MkFitProducer: the key of the table, and its output orders this module after it");
    desc.add<edm::InputTag>("pixelHits", edm::InputTag("mkFitSiPixelHits"));
    desc.add<edm::InputTag>("stripHits", edm::InputTag("mkFitSiPhase2Hits"));
    desc.add<edm::InputTag>("pixelSimLink", edm::InputTag("simSiPixelDigis", "Pixel"));
    desc.add<edm::InputTag>("phase2OTSimLink", edm::InputTag("simSiPixelDigis", "Tracker"));
    desc.add<std::vector<edm::InputTag>>("simHits",
                                         {edm::InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"),
                                          edm::InputTag("g4SimHits", "TrackerHitsPixelBarrelHighTof"),
                                          edm::InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof"),
                                          edm::InputTag("g4SimHits", "TrackerHitsPixelEndcapHighTof")})
        ->setComment("the signal sim hit collections the digitizer read (its ROUList), which the links index");
    desc.add<edm::InputTag>("simTracks", edm::InputTag("g4SimHits"));
    desc.add<edm::InputTag>("simVertices", edm::InputTag("g4SimHits"));
    desc.add<edm::InputTag>("trackingParticles", edm::InputTag(""))
        ->setComment(
            "TrackingParticleRefVector, the sim tracks of MkFitTraceTruth in their order, for the column sim; empty: "
            "sim is -1");
    desc.add<std::string>("linkIndex", "auto")
        ->setComment(
            "what a digi-sim link's sim hit index counts: global, perModule (digitizers before "
            "CMSSW_20_1_0_pre3), or auto (the one more of the event's signal links agree with)");
    desc.addUntracked<int>("printSecondary", 0)->setComment("print up to this many secondary-dominated hits");
    desc.addUntracked<int>("printStride", 1)->setComment("print every printStride-th secondary-dominated hit");
    desc.addUntracked<int>("printUnmatched", 0)
        ->setComment("print up to this many signal links whose sim hit was not found, with the sim hits around");
    descriptions.add("mkFitTraceHitOrigin", desc);
  }

  void MkFitTraceHitOrigin::analyze(edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
    const auto t0 = std::chrono::steady_clock::now();
    edm::Service<MkFitTraceWriter> writer;
    if (!writer.isAvailable())
      throw cms::Exception("Configuration") << "MkFitTraceHitOrigin needs the MkFitTraceWriter service";

    // Only to run after the producer, whose trace the writer has then taken.
    iEvent.get(mkFitOutputToken_);

    const auto& trackerInfo = iSetup.getData(mkFitGeomToken_).trackerInfo();
    const SimLinks& pixelLinks = iEvent.get(pixelSimLinkToken_);
    const SimLinks& phase2Links = iEvent.get(phase2OTSimLinkToken_);
    const edm::SimTrackContainer& simTracks = iEvent.get(simTracksToken_);
    const edm::SimVertexContainer& simVertices = iEvent.get(simVerticesToken_);

    // The signal's SimTracks by id, and the creating process of each (its vertex's process type).
    std::unordered_map<unsigned int, int> stOfId;
    for (int i = 0; i < (int)simTracks.size(); ++i)
      stOfId.emplace(simTracks[i].trackId(), i);
    auto creatorProc = [&](int ist) -> int {
      const int v = simTracks[ist].vertIndex();
      return v >= 0 && v < (int)simVertices.size() ? (int)simVertices[v].processType() : -1;
    };

    // The signal sim hit collections; per (module, collection) the indices of its sim hits in collection order; per
    // (module, track id) its sim hits.
    std::vector<const std::vector<PSimHit>*> coll(simHitTokens_.size());
    std::unordered_map<std::uint64_t, std::vector<int>> hitsOfDetColl;
    std::unordered_map<std::uint64_t, std::vector<std::pair<int, int>>> hitsOfDetTrack;
    long long nSimHits = 0;
    for (int c = 0; c < (int)coll.size(); ++c) {
      coll[c] = &iEvent.get(simHitTokens_[c]);
      nSimHits += coll[c]->size();
      for (int j = 0; j < (int)coll[c]->size(); ++j) {
        const PSimHit& h = (*coll[c])[j];
        hitsOfDetColl[detKey(h.detUnitId(), c)].push_back(j);
        hitsOfDetTrack[detKey(h.detUnitId(), h.trackId())].emplace_back(c, j);
      }
    }
    auto classify = [&](const PSimHit& h, int& ist) -> SimHitClass {
      auto f = stOfId.find(h.trackId());
      if (f == stOfId.end()) {
        ist = -1;
        return kShOrph;
      }
      ist = f->second;
      const SimTrack& t = simTracks[ist];
      return h.particleType() == t.type() && (int)h.processType() == creatorProc(ist) ? kShOwn : kShSec;
    };

    // The sim hit a signal link names, under either decoding of its index: (collection, index), or (-1, -1).
    auto resolve = [&](unsigned int detId, const PixelDigiSimLink& l, bool perModule) -> std::pair<int, int> {
      const unsigned int pos = l.CFposition();
      for (int c = 0; c < (int)coll.size(); ++c) {
        if (simHitTof_[c] != l.TofBin())
          continue;
        auto hs = hitsOfDetColl.find(detKey(detId, c));
        if (hs == hitsOfDetColl.end())
          continue;
        int j = -1;
        if (perModule) {
          const long long k = (long long)pos - hs->second.front();
          if (k >= 0 && k < (long long)hs->second.size())
            j = hs->second[k];
        } else if (pos < coll[c]->size() && (*coll[c])[pos].detUnitId() == detId) {
          j = pos;
        }
        if (j >= 0 && (*coll[c])[j].trackId() == l.SimTrackId())
          return {c, j};
      }
      return {-1, -1};
    };

    // Which decoding: count the signal links of the event each one resolves.
    long long nEvLinksSig = 0, nFound[2] = {0, 0};
    for (const SimLinks* links : {&pixelLinks, &phase2Links})
      for (const auto& ds : *links)
        for (const PixelDigiSimLink& l : ds.data) {
          if (!isSignal(l.eventId()))
            continue;
          ++nEvLinksSig;
          for (int pm = 0; pm < 2; ++pm)
            nFound[pm] += resolve(ds.detId(), l, pm == 1).first >= 0;
        }
    const bool perModule =
        linkIndex_ == LinkIndex::kPerModule || (linkIndex_ == LinkIndex::kAuto && nFound[1] > nFound[0]);

    // TrackingParticles, for the column sim, keyed as MkFitTraceTruth keys them.
    std::unordered_map<std::uint64_t, int> tpOfTrack;
    if (useTrackingParticles_) {
      const TrackingParticleRefVector& tps = iEvent.get(trackingParticlesToken_);
      for (int i = 0; i < (int)tps.size(); ++i) {
        const EncodedEventId eid = tps[i]->eventId();
        for (const SimTrack& g4 : tps[i]->g4Tracks())
          tpOfTrack.emplace(simTrackKey(g4.trackId(), eid), i);
      }
    }

    TraceHitOriginTables T;
    T.run = iEvent.id().run();
    T.lumi = iEvent.id().luminosityBlock();
    T.event = iEvent.id().event();
    T.module = module_;
    T.layers.resize(trackerInfo.n_layers());
    for (int l = 0; l < (int)T.layers.size(); ++l)
      T.layers[l].layer = l;

    long long nHits = 0, nListed = 0, nLinksSig = 0, nLinksUnm = 0, nLinksPu = 0, nOrphDom = 0;
    long long nCls[4] = {0, 0, 0, 0};
    double qOwn = 0, qSec = 0, qOrph = 0, qPu = 0, qUnm = 0, qNoise = 0;
    long long nShClass[3] = {0, 0, 0};  // signal sim hits behind listed hits' charge, by class (counted per hit)
    std::map<std::pair<int, int>, std::pair<long long, double>> secTypes;  // (pdg, proc) -> (sim hits, charge)
    int nPrinted = 0, nUnmPrinted = 0;
    long long nSecSeen = 0;
    // Unmatched signal links on hits, [pixel or outer][reason]: 0 no sim hit on the module in the link's tof bin,
    // 1 the index out of the module's (perModule) or the collection's (global) range, 2 the sim hit at the index on
    // another module (global), 3 another track id; and whether the module has a sim hit of the link's track id.
    long long nUnm[2][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}}, nUnmTrackHere[2] = {0, 0};

    auto doHits = [&](const MkFitClusterIndexToHit& clusterIndex,
                      const std::vector<int>& layerOf,
                      const SimLinks& links,
                      bool pixel) {
      // The links of the current module, sorted by channel.
      unsigned int curDet = 0;
      std::vector<std::pair<unsigned int, const PixelDigiSimLink*>> chanLinks;
      // Per hit: charge per contributing signal sim hit; per (track id, event id), for the TrackingParticle.
      std::map<std::pair<int, int>, double> shCharge;
      std::unordered_map<std::uint64_t, double> keyCharge;

      for (int i = 0; i < (int)layerOf.size(); ++i) {
        const int ilay = layerOf[i];
        if (ilay < 0)
          continue;
        ++nHits;
        auto& L = T.layers.at(ilay);
        L.n_hits = std::max(L.n_hits, i + 1);
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

        // The channels of the cluster and their charge.
        std::vector<std::pair<unsigned int, double>> chans;
        if (pixel) {
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

        double tot = 0, own = 0, sec = 0, orph = 0, pu = 0, unm = 0, noise = 0;
        bool anyLink = false;
        shCharge.clear();
        keyCharge.clear();
        for (const auto& [chan, w] : chans) {
          tot += w;
          auto range = std::equal_range(chanLinks.begin(),
                                        chanLinks.end(),
                                        std::make_pair(chan, (const PixelDigiSimLink*)nullptr),
                                        [](const auto& a, const auto& b) { return a.first < b.first; });
          double sumf = 0;
          for (auto it = range.first; it != range.second; ++it)
            sumf += it->second->fraction();
          // Shares of the channel's signal; scaled to 1 when the noise made the total smaller than the sim hits'.
          const double scale = sumf > 1 ? 1. / sumf : 1.;
          for (auto it = range.first; it != range.second; ++it) {
            const PixelDigiSimLink& l = *it->second;
            anyLink = true;
            const double q = w * l.fraction() * scale;
            const EncodedEventId eid = l.eventId();
            if (!isSignal(eid)) {
              ++nLinksPu;
              pu += q;
              continue;
            }
            ++nLinksSig;
            const auto [cFound, jFound] = resolve(detId, l, perModule);
            if (cFound < 0) {
              ++nLinksUnm;
              unm += q;
              const int po = pixel ? 0 : 1;
              int reason = 0;
              for (int c = 0; c < (int)coll.size(); ++c) {
                if (simHitTof_[c] != l.TofBin())
                  continue;
                auto hs = hitsOfDetColl.find(detKey(detId, c));
                if (hs == hitsOfDetColl.end())
                  continue;
                const unsigned int pos = l.CFposition();
                if (perModule)
                  reason =
                      pos < (unsigned int)hs->second.front() || pos - hs->second.front() >= hs->second.size() ? 1 : 3;
                else
                  reason = pos >= coll[c]->size() ? 1 : (*coll[c])[pos].detUnitId() != detId ? 2 : 3;
                break;
              }
              ++nUnm[po][reason];
              nUnmTrackHere[po] += hitsOfDetTrack.count(detKey(detId, l.SimTrackId())) > 0;
              if (nUnmPrinted < printUnmatched_) {
                ++nUnmPrinted;
                std::ostringstream os;
                os << "MKFITHITORIGIN-UNMLINK det=" << detId << " chan=" << chan << " track=" << l.SimTrackId()
                   << " tof=" << l.TofBin() << " pos=" << l.CFposition() << " fraction=" << l.fraction()
                   << " reason=" << reason << "\n";
                auto th = hitsOfDetTrack.find(detKey(detId, l.SimTrackId()));
                if (th != hitsOfDetTrack.end())
                  for (const auto& [c, j] : th->second)
                    os << "  of the track on the module: " << simHitNames_[c] << "[" << j
                       << "] pdg=" << (*coll[c])[j].particleType() << "\n";
                std::cout << os.str();
              }
              continue;
            }
            const PSimHit& h = (*coll[cFound])[jFound];
            int ist;
            const SimHitClass k = classify(h, ist);
            (k == kShOwn ? own : k == kShSec ? sec : orph) += q;
            shCharge[{cFound, jFound}] += q;
            keyCharge[simTrackKey(l.SimTrackId(), eid)] += q;
          }
          noise += w * std::max(0., 1. - sumf);
        }
        if (!anyLink || tot <= 0)
          continue;

        ++nListed;
        qOwn += own;
        qSec += sec;
        qOrph += orph;
        qPu += pu;
        qUnm += unm;
        qNoise += noise;
        const double fSec = (sec + orph) / tot, fOwn = own / tot, fUnk = (pu + unm) / tot;
        std::uint8_t cls = TraceHitOriginTables::kNone;
        if (fOwn > 0 || fSec > 0 || fUnk > 0) {
          if (fOwn >= fSec && fOwn >= fUnk)
            cls = TraceHitOriginTables::kOwn;
          else if (fSec >= fUnk)
            cls = TraceHitOriginTables::kSecondary;
          else
            cls = TraceHitOriginTables::kUnknown;
        }
        ++nCls[cls];
        if (orph > own && orph > sec && orph > pu + unm)
          ++nOrphDom;

        // The largest secondary sim hit, the signal SimTrack with most charge and its TrackingParticle.
        int secPdg = 0, secProc = 0;
        double secBest = -1;
        for (const auto& [ci, q] : shCharge) {
          const PSimHit& h = (*coll[ci.first])[ci.second];
          int ist;
          const SimHitClass k = classify(h, ist);
          ++nShClass[k];
          if (k == kShOwn)
            continue;
          auto& st = secTypes[{h.particleType(), (int)h.processType()}];
          ++st.first;
          st.second += q;
          if (q > secBest) {
            secBest = q;
            secPdg = h.particleType();
            secProc = h.processType();
          }
        }
        int sim = -1;
        if (useTrackingParticles_) {
          double best = -1;
          std::uint64_t bestKey = 0;
          for (const auto& [k, q] : keyCharge)
            if (q > best || (q == best && k < bestKey)) {
              best = q;
              bestKey = k;
            }
          auto f = tpOfTrack.find(bestKey);
          if (best > 0 && f != tpOfTrack.end())
            sim = f->second;
        }

        L.hit.push_back(i);
        L.sim.push_back(sim);
        L.f_own.push_back(fOwn);
        L.f_sec.push_back(fSec);
        L.f_pu.push_back(pu / tot);
        L.f_unm.push_back(unm / tot);
        L.cls.push_back(cls);
        L.sec_pdg.push_back(secPdg);
        L.sec_proc.push_back(secProc);

        // Evidence for the secondary class: the hit, the sim hits behind it, the sim hits of their SimTracks here.
        if (cls != TraceHitOriginTables::kSecondary)
          continue;
        if (printSecondary_ > 0 && nPrinted < printSecondary_ && nSecSeen % printStride_ == 0) {
          ++nPrinted;
          std::ostringstream os;
          char buf[512];
          const LocalPoint lp = rh->localPosition();
          std::snprintf(buf,
                        sizeof(buf),
                        "MKFITHITORIGIN-HIT %u:%u:%llu layer=%d idx=%d detid=%u %s channels=%zu charge=%.0f f_own=%.3f "
                        "f_sec=%.3f (orphan %.3f) f_pu=%.3f f_unm=%.3f noise=%.3f sim=%d cluster_local=(%.5f,%.5f)",
                        T.run,
                        T.lumi,
                        T.event,
                        ilay,
                        i,
                        detId,
                        pixel ? "pixel" : "outer",
                        chans.size(),
                        tot,
                        fOwn,
                        fSec,
                        orph / tot,
                        pu / tot,
                        unm / tot,
                        noise / tot,
                        sim,
                        lp.x(),
                        lp.y());
          os << buf << "\n";
          std::vector<unsigned int> idsHere;
          for (const auto& [ci, q] : shCharge) {
            const PSimHit& h = (*coll[ci.first])[ci.second];
            int ist;
            const SimHitClass k = classify(h, ist);
            std::snprintf(buf,
                          sizeof(buf),
                          "  MKFITHITORIGIN-SIMHIT %s[%d] class=%s charge_share=%.3f trackId=%u pdg=%d proc=%d "
                          "eloss_keV=%.2f p_GeV=%.5f tof_ns=%.3f entry=%s exit=%s mid_to_cluster_um=%.1f",
                          simHitNames_[ci.first].c_str(),
                          ci.second,
                          kShClassName[k],
                          q / tot,
                          h.trackId(),
                          h.particleType(),
                          (int)h.processType(),
                          h.energyLoss() * 1e6,
                          h.pabs(),
                          h.tof(),
                          fmtPoint(h.entryPoint()).c_str(),
                          fmtPoint(h.exitPoint()).c_str(),
                          1e4 * std::hypot(h.localPosition().x() - lp.x(), h.localPosition().y() - lp.y()));
            os << buf;
            if (ist >= 0) {
              const SimTrack& t = simTracks[ist];
              std::snprintf(buf,
                            sizeof(buf),
                            " | SimTrack pdg=%d creator_proc=%d p_GeV=%.4f pt_GeV=%.4f eta=%.3f",
                            t.type(),
                            creatorProc(ist),
                            t.momentum().P(),
                            t.momentum().Pt(),
                            t.momentum().Eta());
              os << buf;
            } else {
              os << " | no SimTrack with this id";
            }
            os << "\n";
            if (std::find(idsHere.begin(), idsHere.end(), h.trackId()) == idsHere.end())
              idsHere.push_back(h.trackId());
          }
          for (unsigned int id : idsHere) {
            auto th = hitsOfDetTrack.find(detKey(detId, id));
            if (th == hitsOfDetTrack.end())
              continue;
            for (const auto& [c, j] : th->second) {
              const PSimHit& h = (*coll[c])[j];
              int ist2;
              const SimHitClass k = classify(h, ist2);
              std::snprintf(buf,
                            sizeof(buf),
                            "  MKFITHITORIGIN-ONMODULE %s[%d] class=%s trackId=%u pdg=%d proc=%d eloss_keV=%.2f "
                            "p_GeV=%.5f entry=%s exit=%s mid_to_cluster_um=%.1f%s",
                            simHitNames_[c].c_str(),
                            j,
                            kShClassName[k],
                            id,
                            h.particleType(),
                            (int)h.processType(),
                            h.energyLoss() * 1e6,
                            h.pabs(),
                            fmtPoint(h.entryPoint()).c_str(),
                            fmtPoint(h.exitPoint()).c_str(),
                            1e4 * std::hypot(h.localPosition().x() - lp.x(), h.localPosition().y() - lp.y()),
                            shCharge.count({c, j}) ? " (in this hit)" : "");
              os << buf << "\n";
            }
          }
          std::cout << os.str();
        }
        ++nSecSeen;
      }
    };
    doHits(iEvent.get(pixelClusterIndexToken_), iEvent.get(pixelLayerToken_), pixelLinks, true);
    doHits(iEvent.get(stripClusterIndexToken_), iEvent.get(stripLayerToken_), phase2Links, false);

    const int pass = writer->writeHitOrigin(T);

    const double qTot = qOwn + qSec + qOrph + qPu + qUnm + qNoise;
    auto share = [&](double q) { return qTot > 0 ? q / qTot : 0.; };
    char buf[1024];
    std::snprintf(buf,
                  sizeof(buf),
                  "MKFITHITORIGIN %u:%u:%llu %s hits=%lld listed=%lld cls_own=%lld cls_sec=%lld cls_unk=%lld "
                  "cls_none=%lld orph_dom=%lld q_own=%.4f q_sec=%.4f q_orph=%.4f q_pu=%.4f q_unm=%.4f q_noise=%.4f "
                  "links_sig=%lld links_unmatched=%lld links_pu=%lld link_index=%s event_links_sig=%lld "
                  "found_global=%lld found_perModule=%lld simtracks=%zu simhits=%lld simhits_own=%lld "
                  "simhits_sec=%lld simhits_orph=%lld pass=%d ms=%.1f",
                  T.run,
                  T.lumi,
                  T.event,
                  module_.c_str(),
                  nHits,
                  nListed,
                  nCls[1],
                  nCls[2],
                  nCls[3],
                  nCls[0],
                  nOrphDom,
                  share(qOwn),
                  share(qSec),
                  share(qOrph),
                  share(qPu),
                  share(qUnm),
                  share(qNoise),
                  nLinksSig,
                  nLinksUnm,
                  nLinksPu,
                  perModule ? "perModule" : "global",
                  nEvLinksSig,
                  nFound[0],
                  nFound[1],
                  simTracks.size(),
                  nSimHits,
                  nShClass[kShOwn],
                  nShClass[kShSec],
                  nShClass[kShOrph],
                  pass,
                  std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
    std::ostringstream os;
    os << buf << "\n";
    // The most frequent secondary types, by sim hits.
    std::vector<std::pair<std::pair<int, int>, std::pair<long long, double>>> st(secTypes.begin(), secTypes.end());
    std::sort(st.begin(), st.end(), [](const auto& a, const auto& b) { return a.second.first > b.second.first; });
    double qSecAll = 0;
    for (const auto& e : st)
      qSecAll += e.second.second;
    os << "MKFITHITORIGIN-SECTYPES " << T.run << ":" << T.lumi << ":" << T.event << " " << module_;
    for (std::size_t k = 0; k < std::min<std::size_t>(st.size(), 12); ++k) {
      std::snprintf(buf,
                    sizeof(buf),
                    " %d/%d:%lld:%.3f",
                    st[k].first.first,
                    st[k].first.second,
                    st[k].second.first,
                    qSecAll > 0 ? st[k].second.second / qSecAll : 0.);
      os << buf;
    }
    os << "\n";
    os << "MKFITHITORIGIN-UNMATCHED " << T.run << ":" << T.lumi << ":" << T.event << " " << module_;
    for (int po = 0; po < 2; ++po) {
      std::snprintf(buf,
                    sizeof(buf),
                    " %s: no_simhit_in_tof_bin=%lld out_of_range=%lld other_module=%lld other_track=%lld "
                    "track_has_simhit_here=%lld;",
                    po ? "outer" : "pixel",
                    nUnm[po][0],
                    nUnm[po][1],
                    nUnm[po][2],
                    nUnm[po][3],
                    nUnmTrackHere[po]);
      os << buf;
    }
    os << "\n";
    std::cout << os.str() << std::flush;
  }

}  // namespace mkfit_fwff

using mkfit_fwff::MkFitTraceHitOrigin;
DEFINE_FWK_MODULE(MkFitTraceHitOrigin);
