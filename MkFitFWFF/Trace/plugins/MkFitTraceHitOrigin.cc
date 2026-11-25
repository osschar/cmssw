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
// - Any other link is pileup. A RAW file keeps only the signal's sim hits, so without more the pileup's charge cannot
//   be classified. With the pileup's sim hits (pileupSimHits, pileupSimTracks, pileupSimVertices: the vectors
//   MkFitTraceCrossingFrame* copy from a playback of the event's mixing, signal first, then the pileup events in the
//   order the mixing accumulated them), a pileup link is resolved as a signal link, in the vector of its tof bin. The
//   digitizer counted a pileup event's sim hits after the signal's and after those of the pileup events before it, so
//   "global" is the vector's index, and "perModule" the vector index of the module's first sim hit of that pileup event
//   plus the rank among the module's sim hits of that event; the decoding the signal links chose is used. The sim
//   hit is classified as a signal one, with the SimTracks of its own pileup event (the event id: bunch crossing and
//   event number) and their creating processes (a pileup SimTrack's vertex index counts in its own event's
//   SimVertices). A pileup link whose sim hit is not found so, or any pileup link without the pileup's sim hits, is
//   "pileup unmatched".
// - Per hit: f_own, f_sec (secondary and orphan), f_pu, f_unm, the fractions of the hit's charge, each link
//   contributing its share times the channel's charge; the rest is charge without a link (noise). cls: the largest of
//   own, secondary and unknown (pileup + unmatched), own winning a tie, then secondary. sim: the sim track (the index
//   into the trace file's SimTracks, the TrackingParticles of trackingParticles in their order, as MkFitTraceTruth
//   writes them) of the signal SimTrack with most of the hit's charge, own and secondary together; -1 without
//   trackingParticles. sec_pdg, sec_proc: the secondary sim hit with most charge in the hit. These describe the
//   signal's charge only, with or without the pileup's sim hits.
// - The pileup's part: f_pu_own, f_pu_sec (secondary and orphan), f_pu_unm, the pileup's charge f_pu split by the same
//   rule; cls_all, the largest of own (f_own + f_pu_own), secondary (f_sec + f_pu_sec) and unknown (f_unm +
//   f_pu_unm), with cls's tie rule; sim_all and sec_pdg_all, sec_proc_all as sim and sec_pdg, sec_proc over the signal
//   and the pileup (a pileup SimTrack's sim track is the TrackingParticle with its track id and event id, in-time
//   pileup only when trackingParticles is in-time only). Without the pileup's sim hits they reduce to f_pu_unm = f_pu,
//   cls_all = cls, sim_all = sim, sec_*_all = sec_*.
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
// The pileup's lines (also without the pileup's sim hits, pileup=0):
//   MKFITHITORIGIN-PU <key> pileup=<0|1> cls_all_own=<n> cls_all_sec=<n> cls_all_unk=<n> cls_all_none=<n>
//     cls_all_sec_pu=<of cls_all_sec: the pileup's secondary charge the larger part> cls_all_own_pu=<likewise for own>
//     q_pu_own=.. q_pu_sec=.. q_pu_orph=.. q_pu_unm=.. (shares of the listed hits' total charge) links_pu=<on hits>
//     links_pu_unmatched=<n> event_links_pu=<all pileup links of the event> found_pu_global=<n>
//     found_pu_perModule=<n> pu_simhits=<pileup sim hits read> pu_simtracks=<n> pu_simvertices=<n>
//     pu_events=<distinct pileup event ids of the SimTracks> pu_vert_bad=<pileup SimTracks whose vertex index is not in
//     their event's SimVertices> parent_missing_sig=<signal SimVertices whose parent track id names no SimTrack of the
//     event>/<with a parent> parent_missing_pu=<likewise, pileup> sig_prefix_mismatch=<signal entries of the pileup
//     vectors that differ from the signal collections> simhits_pu_own=.. simhits_pu_sec=.. simhits_pu_orph=..
//   MKFITHITORIGIN-PUEVENTS <key> per bunch crossing: <bcr>:<event ids with SimTracks>/<largest event number>/<event
//     ids with tracker sim hits>
//   MKFITHITORIGIN-PLAYBACK <key> <InputTag> min_bunch=.. bunches=.. sources=.. events=<n> per_bcr=<source 0's count
//     per bunch crossing> same_as_first=<1|0: sizes and event ids equal those of the first playbackInfo>
//   MKFITHITORIGIN-PUSECTYPES <key> as -SECTYPES, for the pileup's secondary sim hits on hits
//   MKFITHITORIGIN-PUUNMATCHED <key> pixel: no_simhit_of_event_on_module=.. out_of_range=.. other_module=..
//     other_track=..; outer: ...  (pileup links on hits whose sim hit was not found, with the pileup's sim hits)
// With printPileupSecondary > 0, for every printStride-th hit of class cls_all secondary whose secondary charge is
// mostly the pileup's, up to printPileupSecondary of them, the hit and the pileup sim hits behind it
// (MKFITHITORIGIN-PUHIT, -PUSIMHIT lines).
// With printSecondary > 0, for every printStride-th hit whose class is secondary, up to printSecondary of them, the
// hit and each sim hit behind its signal charge (MKFITHITORIGIN-HIT, -SIMHIT lines) and every sim hit of the same
// SimTracks on the same module (-ONMODULE lines), with their local positions against the cluster's. With
// printUnmatched > 0, up to that many unmatched signal links with the sim hits around their index (-UNMLINK).
//
// The per-hit list (listPrefix not empty): three CSV files, for a comparison with Geant4's own steps (G4Snitch), which
// know the stepping track's own id, the DetId and the global step points, but not the clusters. Every row starts with
// tag (listTag), run, lumi, event, pass (the table's pass). Lengths in cm, in the module's local frame unless named g*.
//   <listPrefix>hits.csv       one row per listed hit (every hit with a link): layer, idx (as in the trace file),
//     detid, subdet, is_pixel, n_chan, charge (ADC sum, or strips); the cluster's channels as edges in the topology's
//     measurement units, rows [row0, row1) and columns [col0, col1) (an outer-tracker cluster has one column), and the
//     same edges in local x and y (lx0, lx1, ly0, ly1, from PixelTopology::localX/localY; the charge was deposited
//     before the drift, so allow for the Lorentz shift); the rec hit's local position and errors (hit_lx, hit_ly,
//     hit_ex, hit_ey) and its global position (hit_gx, hit_gy, hit_gz = Surface::toGlobal of the local one); the
//     module's frame: origin (mod_gx, mod_gy, mod_gz) and its local axes in global coordinates (ux_*, uy_*, uz_*, the
//     rows of Surface::rotation()), so global = origin + lx ux + ly uy + lz uz and lx = ux . (global - origin); the
//     sensor's width, length and thickness (local z within +-thickness/2); cls, mixed (own and secondary charge both
//     above 0), f_own, f_sec (with orphans), f_orph, f_pu, f_unm, f_noise, sim; n_simhits (signal sim hits behind the
//     charge); the sim hit with most charge (dom_*: its track id, particle type, process, class 0 own 1 secondary
//     2 orphan, charge share, and the type of the SimTrack with its id, 0 if none) and the secondary or orphan sim hit
//     with most charge (sec_*, track id -1 if none).
//   <listPrefix>simhits.csv    one row per (listed hit, signal sim hit behind its charge): layer, idx, detid, the sim
//     hit's collection and index, track id, type, process, class, charge share, energy loss (keV), momentum (GeV),
//     time of flight (ns), local entry and exit points, and the type and creating process of the SimTrack with its id.
//   <listPrefix>simtracks.csv  one row per SimTrack of the event: index, track id (Geant4's), type, charge, four-momentum
//     (GeV), its vertex (index, position, time in ns, process type, the parent's track id: SimVertex::parentIndex(),
//     -1 for a primary), the generator particle index, and its signal sim hits in the tracker collections read.
//   The list is of the signal only: f_pu is the whole pileup's charge, also when the pileup's sim hits are given.
// Process types are Geant4 sub-types (PSimHit::processType(), SimVertex::processType()). One MKFITHITLIST line per
// event gives the counts.

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
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PixelChannelIdentifier.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFramePlaybackInfoNew.h"
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
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
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
    // The pileup's sim hits (per entry of simHits), SimTracks and SimVertices, signal first; empty: none.
    std::vector<edm::EDGetTokenT<std::vector<PSimHit>>> puSimHitTokens_;
    edm::EDGetTokenT<edm::SimTrackContainer> puSimTracksToken_;
    edm::EDGetTokenT<edm::SimVertexContainer> puSimVerticesToken_;
    std::vector<edm::InputTag> playbackTags_;
    std::vector<edm::EDGetTokenT<CrossingFramePlaybackInfoNew>> playbackTokens_;
    const int printSecondary_;
    const int printStride_;
    const int printUnmatched_;
    const int printPileupSecondary_;
    const edm::ESGetToken<MkFitGeometry, TrackerRecoGeometryRecord> mkFitGeomToken_;
    // The per-hit list: the files, written by one event at a time.
    const std::string listPrefix_;
    const std::string listTag_;
    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeomToken_;
    mutable std::mutex listMutex_;
    mutable std::ofstream listHits_, listSimHits_, listSimTracks_;
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
        printPileupSecondary_(iConfig.getUntrackedParameter<int>("printPileupSecondary")),
        mkFitGeomToken_(esConsumes()),
        listPrefix_(iConfig.getUntrackedParameter<std::string>("listPrefix")),
        listTag_(iConfig.getUntrackedParameter<std::string>("listTag")),
        trackerGeomToken_(esConsumes()) {
    if (!listPrefix_.empty()) {
      if (listTag_.find_first_of(",\"\n") != std::string::npos)
        throw cms::Exception("Configuration") << "MkFitTraceHitOrigin: listTag must not hold a comma, quote or newline";
      const char* key = "tag,run,lumi,event,pass,";
      auto open = [&](std::ofstream& f, const char* name, const char* columns) {
        f.open(listPrefix_ + name);
        if (!f)
          throw cms::Exception("Configuration") << "MkFitTraceHitOrigin: cannot open " << listPrefix_ << name;
        f << key << columns << "\n";
      };
      open(listHits_,
           "hits.csv",
           "layer,idx,detid,subdet,is_pixel,n_chan,charge,row0,row1,col0,col1,lx0,lx1,ly0,ly1,hit_lx,hit_ly,hit_ex,"
           "hit_ey,hit_gx,hit_gy,hit_gz,mod_gx,mod_gy,mod_gz,ux_x,ux_y,ux_z,uy_x,uy_y,uy_z,uz_x,uz_y,uz_z,mod_width,"
           "mod_length,mod_thick,cls,mixed,f_own,f_sec,f_orph,f_pu,f_unm,f_noise,sim,n_simhits,dom_trackid,dom_pdg,"
           "dom_proc,dom_class,dom_share,dom_st_pdg,sec_trackid,sec_pdg,sec_proc,sec_class,sec_share");
      open(listSimHits_,
           "simhits.csv",
           "layer,idx,detid,coll,simhit,trackid,pdg,proc,class,share,eloss_keV,p_GeV,tof_ns,entry_lx,entry_ly,entry_lz,"
           "exit_lx,exit_ly,exit_lz,st_pdg,st_proc");
      open(listSimTracks_,
           "simtracks.csv",
           "index,trackid,pdg,charge,px,py,pz,e,vert,vx,vy,vz,vt_ns,vert_proc,parent_trackid,genpart,n_tk_simhits");
    }
    for (const edm::InputTag& t : iConfig.getParameter<std::vector<edm::InputTag>>("simHits")) {
      simHitTokens_.push_back(consumes<std::vector<PSimHit>>(t));
      simHitNames_.push_back(t.instance());
      simHitTof_.push_back(t.instance().find("HighTof") != std::string::npos ? PixelDigiSimLink::HighTof
                                                                             : PixelDigiSimLink::LowTof);
    }
    if (useTrackingParticles_)
      trackingParticlesToken_ =
          consumes<TrackingParticleRefVector>(iConfig.getParameter<edm::InputTag>("trackingParticles"));
    const auto puHits = iConfig.getParameter<std::vector<edm::InputTag>>("pileupSimHits");
    const auto puTracks = iConfig.getParameter<edm::InputTag>("pileupSimTracks");
    const auto puVertices = iConfig.getParameter<edm::InputTag>("pileupSimVertices");
    if (!puHits.empty()) {
      if (puHits.size() != simHitTokens_.size() || puTracks.label().empty() || puVertices.label().empty())
        throw cms::Exception("Configuration")
            << "MkFitTraceHitOrigin: pileupSimHits needs one entry per simHits entry, "
               "and pileupSimTracks and pileupSimVertices";
      for (const edm::InputTag& t : puHits)
        puSimHitTokens_.push_back(consumes<std::vector<PSimHit>>(t));
      puSimTracksToken_ = consumes<edm::SimTrackContainer>(puTracks);
      puSimVerticesToken_ = consumes<edm::SimVertexContainer>(puVertices);
    } else if (!puTracks.label().empty() || !puVertices.label().empty()) {
      throw cms::Exception("Configuration")
          << "MkFitTraceHitOrigin: pileupSimTracks and pileupSimVertices need pileupSimHits";
    }
    for (const edm::InputTag& t : iConfig.getParameter<std::vector<edm::InputTag>>("playbackInfo")) {
      playbackTags_.push_back(t);
      playbackTokens_.push_back(consumes<CrossingFramePlaybackInfoNew>(t));
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
    desc.add<std::vector<edm::InputTag>>("pileupSimHits", {})
        ->setComment(
            "per entry of simHits, in its order: the sim hits of the whole crossing, signal first, then the pileup "
            "events in the mixing's order (MkFitTraceCrossingFramePSimHits); empty: the pileup is not classified");
    desc.add<edm::InputTag>("pileupSimTracks", edm::InputTag(""))
        ->setComment("the SimTracks of the whole crossing, signal first (MkFitTraceCrossingFrameSimTracks)");
    desc.add<edm::InputTag>("pileupSimVertices", edm::InputTag(""))
        ->setComment("the SimVertices of the whole crossing, signal first (MkFitTraceCrossingFrameSimVertices)");
    desc.add<std::vector<edm::InputTag>>("playbackInfo", {})
        ->setComment("mixing playback records (CrossingFramePlaybackInfoNew) to print and compare, e.g. mix::HLT");
    desc.addUntracked<int>("printSecondary", 0)->setComment("print up to this many secondary-dominated hits");
    desc.addUntracked<int>("printStride", 1)->setComment("print every printStride-th secondary-dominated hit");
    desc.addUntracked<int>("printUnmatched", 0)
        ->setComment("print up to this many signal links whose sim hit was not found, with the sim hits around");
    desc.addUntracked<int>("printPileupSecondary", 0)
        ->setComment("print up to this many hits whose charge is mostly a pileup secondary's (every printStride-th)");
    desc.addUntracked<std::string>("listPrefix", "")
        ->setComment("write the per-hit list to <listPrefix>hits.csv, simhits.csv and simtracks.csv; empty: none");
    desc.addUntracked<std::string>("listTag", "")->setComment("the first column of every list row, e.g. the sample");
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
    std::unordered_map<unsigned int, int> nSimHitsOfTrack;  // for the list's simtracks.csv
    long long nSimHits = 0;
    for (int c = 0; c < (int)coll.size(); ++c) {
      coll[c] = &iEvent.get(simHitTokens_[c]);
      nSimHits += coll[c]->size();
      for (int j = 0; j < (int)coll[c]->size(); ++j) {
        const PSimHit& h = (*coll[c])[j];
        hitsOfDetColl[detKey(h.detUnitId(), c)].push_back(j);
        hitsOfDetTrack[detKey(h.detUnitId(), h.trackId())].emplace_back(c, j);
        ++nSimHitsOfTrack[h.trackId()];
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

    // The pileup (with pileupSimHits): per collection the vector of the whole crossing, signal first; per (module,
    // event id) its pileup sim hits' indices in that vector, in order; per (track id, event id) the pileup SimTrack's
    // particle type and creating process (from the SimVertices of its own event).
    const bool withPu = !puSimHitTokens_.empty();
    std::vector<const std::vector<PSimHit>*> puColl(coll.size(), nullptr);
    std::vector<std::unordered_map<std::uint64_t, std::vector<int>>> puHitsOfDetEvent(coll.size());
    struct PuTrack {
      int type, proc;
      const SimTrack* st;
    };
    std::unordered_map<std::uint64_t, PuTrack> puTrackOf;
    long long nPuSimHits = 0, nPuSimTracks = 0, nPuSimVertices = 0, nPuVertBad = 0, nSigPrefixMismatch = 0;
    long long nParentSig = 0, nParentMissSig = 0, nParentPu = 0, nParentMissPu = 0;
    // Per bunch crossing: the event numbers with SimTracks, with tracker sim hits.
    std::map<int, std::set<int>> puEvTracks, puEvHits;
    for (int c = 0; c < (int)coll.size() && withPu; ++c) {
      const std::vector<PSimHit>& v = iEvent.get(puSimHitTokens_[c]);
      puColl[c] = &v;
      const std::size_t nSig = coll[c]->size();
      if (v.size() < nSig)
        throw cms::Exception("LogicError") << "MkFitTraceHitOrigin: pileupSimHits[" << c << "] has " << v.size()
                                           << " sim hits, fewer than the signal's " << nSig;
      for (std::size_t j = 0; j < nSig; ++j) {
        const PSimHit &a = v[j], &b = (*coll[c])[j];
        nSigPrefixMismatch += a.detUnitId() != b.detUnitId() || a.trackId() != b.trackId() ||
                              a.eventId().rawId() != b.eventId().rawId() || a.energyLoss() != b.energyLoss();
      }
      for (int j = nSig; j < (int)v.size(); ++j) {
        const PSimHit& h = v[j];
        puHitsOfDetEvent[c][detKey(h.detUnitId(), h.eventId().rawId())].push_back(j);
        puEvHits[h.eventId().bunchCrossing()].insert(h.eventId().event());
        ++nPuSimHits;
      }
    }
    if (withPu) {
      const edm::SimTrackContainer& tr = iEvent.get(puSimTracksToken_);
      const edm::SimVertexContainer& vx = iEvent.get(puSimVerticesToken_);
      if (tr.size() < simTracks.size() || vx.size() < simVertices.size())
        throw cms::Exception("LogicError") << "MkFitTraceHitOrigin: pileupSimTracks or pileupSimVertices shorter than "
                                              "the signal's";
      for (std::size_t j = 0; j < simTracks.size(); ++j)
        nSigPrefixMismatch += tr[j].trackId() != simTracks[j].trackId() || tr[j].type() != simTracks[j].type();
      for (std::size_t j = 0; j < simVertices.size(); ++j)
        nSigPrefixMismatch += vx[j].processType() != simVertices[j].processType();
      // The SimVertices of each pileup event, in order: a pileup SimTrack's vertex index counts in them.
      std::unordered_map<std::uint32_t, std::vector<const SimVertex*>> vertsOf;
      for (std::size_t j = simVertices.size(); j < vx.size(); ++j) {
        vertsOf[vx[j].eventId().rawId()].push_back(&vx[j]);
        ++nPuSimVertices;
      }
      for (std::size_t j = simTracks.size(); j < tr.size(); ++j) {
        const SimTrack& t = tr[j];
        const EncodedEventId eid = t.eventId();
        const auto& vs = vertsOf[eid.rawId()];
        const int vi = t.vertIndex();
        int proc = -1;
        if (vi >= 0 && vi < (int)vs.size())
          proc = (int)vs[vi]->processType();
        else
          ++nPuVertBad;
        puTrackOf.emplace(simTrackKey(t.trackId(), eid), PuTrack{t.type(), proc, &t});
        puEvTracks[eid.bunchCrossing()].insert(eid.event());
        ++nPuSimTracks;
      }
      // Is the grouping by event right: does a vertex's parent track id name a SimTrack of the same event?
      for (const SimVertex& v : simVertices)
        if (v.parentIndex() >= 0) {
          ++nParentSig;
          nParentMissSig += stOfId.count(v.parentIndex()) == 0;
        }
      for (std::size_t j = simVertices.size(); j < vx.size(); ++j)
        if (vx[j].parentIndex() >= 0) {
          ++nParentPu;
          nParentMissPu += puTrackOf.count(simTrackKey(vx[j].parentIndex(), vx[j].eventId())) == 0;
        }
    }
    auto classifyPu = [&](const PSimHit& h, const PuTrack*& pt) -> SimHitClass {
      auto f = puTrackOf.find(simTrackKey(h.trackId(), h.eventId()));
      if (f == puTrackOf.end()) {
        pt = nullptr;
        return kShOrph;
      }
      pt = &f->second;
      return h.particleType() == pt->type && (int)h.processType() == pt->proc ? kShOwn : kShSec;
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

    // The pileup sim hit a pileup link names, under either decoding of its index: (collection, index into the
    // collection's vector of the whole crossing), or (-1, -1). Without the pileup's sim hits always (-1, -1).
    auto resolvePu = [&](unsigned int detId, const PixelDigiSimLink& l, bool perModule) -> std::pair<int, int> {
      if (!withPu)
        return {-1, -1};
      const unsigned int pos = l.CFposition();
      const std::uint32_t eid = l.eventId().rawId();
      for (int c = 0; c < (int)coll.size(); ++c) {
        if (simHitTof_[c] != l.TofBin())
          continue;
        auto hs = puHitsOfDetEvent[c].find(detKey(detId, eid));
        if (hs == puHitsOfDetEvent[c].end())
          continue;
        const std::vector<PSimHit>& v = *puColl[c];
        int j = -1;
        if (perModule) {
          const long long k = (long long)pos - hs->second.front();
          if (k >= 0 && k < (long long)hs->second.size())
            j = hs->second[k];
        } else if (pos >= coll[c]->size() && pos < v.size() && v[pos].detUnitId() == detId &&
                   v[pos].eventId().rawId() == eid) {
          j = pos;
        }
        if (j >= 0 && v[j].trackId() == l.SimTrackId())
          return {c, j};
      }
      return {-1, -1};
    };

    // Which decoding: count the signal links of the event each one resolves (and the pileup links, for the log).
    long long nEvLinksSig = 0, nFound[2] = {0, 0}, nEvLinksPu = 0, nFoundPu[2] = {0, 0};
    for (const SimLinks* links : {&pixelLinks, &phase2Links})
      for (const auto& ds : *links)
        for (const PixelDigiSimLink& l : ds.data) {
          if (!isSignal(l.eventId())) {
            ++nEvLinksPu;
            for (int pm = 0; pm < 2 && withPu; ++pm)
              nFoundPu[pm] += resolvePu(ds.detId(), l, pm == 1).first >= 0;
            continue;
          }
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
    // The pileup's part: unmatched pileup links on hits (with the pileup's sim hits; reasons as nUnm, 0 meaning no sim
    // hit of the link's event on the module in its tof bin), cls_all, charge, sim hits by class, secondary types.
    long long nLinksPuUnm = 0, nPuUnm[2][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}};
    long long nClsAll[4] = {0, 0, 0, 0}, nClsAllSecPu = 0, nClsAllOwnPu = 0;
    double qPuOwn = 0, qPuSec = 0, qPuOrph = 0, qPuUnm = 0;
    long long nShClassPu[3] = {0, 0, 0};
    std::map<std::pair<int, int>, std::pair<long long, double>> puSecTypes;
    int nPuPrinted = 0;
    long long nPuSecSeen = 0;

    // The per-hit list: rows without their key (tag, run, lumi, event, pass), written at the end of the event.
    const bool listing = !listPrefix_.empty();
    const TrackerGeometry* tkGeom = listing ? &iSetup.getData(trackerGeomToken_) : nullptr;
    std::vector<std::string> rowsHits, rowsSimHits, rowsSimTracks;
    long long nListMixed = 0, nListOwnOnly = 0, nListSecOnly = 0, nListNeither = 0;
    auto stPdg = [&](unsigned int trackId) -> int {
      auto f = stOfId.find(trackId);
      return f == stOfId.end() ? 0 : simTracks[f->second].type();
    };

    auto doHits = [&](const MkFitClusterIndexToHit& clusterIndex,
                      const std::vector<int>& layerOf,
                      const SimLinks& links,
                      bool pixel) {
      // The links of the current module, sorted by channel.
      unsigned int curDet = 0;
      std::vector<std::pair<unsigned int, const PixelDigiSimLink*>> chanLinks;
      // Per hit: charge per contributing signal sim hit; per (track id, event id), for the TrackingParticle.
      std::map<std::pair<int, int>, double> shCharge;
      // The same for the pileup's sim hits (collection, index into its vector of the whole crossing), and the charge
      // per (track id, event id) over the signal and the pileup, for sim_all.
      std::map<std::pair<int, int>, double> puShCharge;
      std::unordered_map<std::uint64_t, double> keyChargeAll;
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
        double puOwn = 0, puSec = 0, puOrph = 0, puUnm = 0;
        bool anyLink = false;
        shCharge.clear();
        keyCharge.clear();
        puShCharge.clear();
        keyChargeAll.clear();
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
              const auto [cPu, jPu] = resolvePu(detId, l, perModule);
              if (cPu < 0) {
                puUnm += q;
                if (withPu) {
                  ++nLinksPuUnm;
                  int reason = 0;
                  for (int c = 0; c < (int)coll.size(); ++c) {
                    if (simHitTof_[c] != l.TofBin())
                      continue;
                    auto hs = puHitsOfDetEvent[c].find(detKey(detId, eid.rawId()));
                    if (hs == puHitsOfDetEvent[c].end())
                      continue;
                    const unsigned int pos = l.CFposition();
                    const std::vector<PSimHit>& v = *puColl[c];
                    if (perModule)
                      reason = pos < (unsigned int)hs->second.front() || pos - hs->second.front() >= hs->second.size()
                                   ? 1
                                   : 3;
                    else
                      reason = pos >= v.size() || pos < coll[c]->size()                                 ? 1
                               : v[pos].detUnitId() != detId || v[pos].eventId().rawId() != eid.rawId() ? 2
                                                                                                        : 3;
                    break;
                  }
                  ++nPuUnm[pixel ? 0 : 1][reason];
                }
                continue;
              }
              const PSimHit& h = (*puColl[cPu])[jPu];
              const PuTrack* pt;
              const SimHitClass k = classifyPu(h, pt);
              (k == kShOwn ? puOwn : k == kShSec ? puSec : puOrph) += q;
              puShCharge[{cPu, jPu}] += q;
              keyChargeAll[simTrackKey(l.SimTrackId(), eid)] += q;
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
            keyChargeAll[simTrackKey(l.SimTrackId(), eid)] += q;
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

        // The pileup's part and the class over all the charge. Without the pileup's sim hits puUnm is pu, summed in
        // the same order, and these equal the signal's fields (cls_all = cls, ...).
        qPuOwn += puOwn;
        qPuSec += puSec;
        qPuOrph += puOrph;
        qPuUnm += puUnm;
        const double fOwnAll = (own + puOwn) / tot, fSecAll = ((sec + orph) + (puSec + puOrph)) / tot,
                     fUnkAll = (unm + puUnm) / tot;
        std::uint8_t clsAll = TraceHitOriginTables::kNone;
        if (fOwnAll > 0 || fSecAll > 0 || fUnkAll > 0) {
          if (fOwnAll >= fSecAll && fOwnAll >= fUnkAll)
            clsAll = TraceHitOriginTables::kOwn;
          else if (fSecAll >= fUnkAll)
            clsAll = TraceHitOriginTables::kSecondary;
          else
            clsAll = TraceHitOriginTables::kUnknown;
        }
        ++nClsAll[clsAll];
        const bool secMostlyPu = puSec + puOrph > sec + orph;
        nClsAllSecPu += clsAll == TraceHitOriginTables::kSecondary && secMostlyPu;
        nClsAllOwnPu += clsAll == TraceHitOriginTables::kOwn && puOwn > own;
        // The largest secondary sim hit over the signal (first, as above) and the pileup.
        int secPdgAll = secPdg, secProcAll = secProc;
        for (const auto& [ci, q] : puShCharge) {
          const PSimHit& h = (*puColl[ci.first])[ci.second];
          const PuTrack* pt;
          const SimHitClass k = classifyPu(h, pt);
          ++nShClassPu[k];
          if (k == kShOwn)
            continue;
          auto& st = puSecTypes[{h.particleType(), (int)h.processType()}];
          ++st.first;
          st.second += q;
          if (q > secBest) {
            secBest = q;
            secPdgAll = h.particleType();
            secProcAll = h.processType();
          }
        }
        int simAll = -1;
        if (useTrackingParticles_) {
          double best = -1;
          std::uint64_t bestKey = 0;
          for (const auto& [k, q] : keyChargeAll)
            if (q > best || (q == best && k < bestKey)) {
              best = q;
              bestKey = k;
            }
          auto f = tpOfTrack.find(bestKey);
          if (best > 0 && f != tpOfTrack.end())
            simAll = f->second;
        }
        L.f_pu_own.push_back(puOwn / tot);
        L.f_pu_sec.push_back((puSec + puOrph) / tot);
        L.f_pu_unm.push_back(puUnm / tot);
        L.cls_all.push_back(clsAll);
        L.sim_all.push_back(simAll);
        L.sec_pdg_all.push_back(secPdgAll);
        L.sec_proc_all.push_back(secProcAll);

        // Evidence for the pileup's secondary class: the hit and the pileup sim hits behind its charge.
        if (clsAll == TraceHitOriginTables::kSecondary && secMostlyPu) {
          if (printPileupSecondary_ > 0 && nPuPrinted < printPileupSecondary_ && nPuSecSeen % printStride_ == 0) {
            ++nPuPrinted;
            std::ostringstream os;
            char buf[512];
            const LocalPoint lp = rh->localPosition();
            std::snprintf(buf,
                          sizeof(buf),
                          "MKFITHITORIGIN-PUHIT %u:%u:%llu layer=%d idx=%d detid=%u %s channels=%zu charge=%.0f "
                          "f_own=%.3f f_sec=%.3f f_pu_own=%.3f f_pu_sec=%.3f f_pu_unm=%.3f f_unm=%.3f sim_all=%d "
                          "cluster_local=(%.5f,%.5f)",
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
                          puOwn / tot,
                          (puSec + puOrph) / tot,
                          puUnm / tot,
                          unm / tot,
                          simAll,
                          lp.x(),
                          lp.y());
            os << buf << "\n";
            for (const auto& [ci, q] : puShCharge) {
              const PSimHit& h = (*puColl[ci.first])[ci.second];
              const PuTrack* pt;
              const SimHitClass k = classifyPu(h, pt);
              std::snprintf(buf,
                            sizeof(buf),
                            "  MKFITHITORIGIN-PUSIMHIT %s[%d] class=%s charge_share=%.3f bcr=%d event=%d trackId=%u "
                            "pdg=%d proc=%d eloss_keV=%.2f p_GeV=%.5f tof_ns=%.3f mid_to_cluster_um=%.1f",
                            simHitNames_[ci.first].c_str(),
                            ci.second,
                            kShClassName[k],
                            q / tot,
                            h.eventId().bunchCrossing(),
                            h.eventId().event(),
                            h.trackId(),
                            h.particleType(),
                            (int)h.processType(),
                            h.energyLoss() * 1e6,
                            h.pabs(),
                            h.tof(),
                            1e4 * std::hypot(h.localPosition().x() - lp.x(), h.localPosition().y() - lp.y()));
              os << buf;
              if (pt) {
                std::snprintf(buf,
                              sizeof(buf),
                              " | SimTrack pdg=%d creator_proc=%d p_GeV=%.4f pt_GeV=%.4f eta=%.3f",
                              pt->type,
                              pt->proc,
                              pt->st->momentum().P(),
                              pt->st->momentum().Pt(),
                              pt->st->momentum().Eta());
                os << buf;
              } else {
                os << " | no SimTrack with this id in its event";
              }
              os << "\n";
            }
            std::cout << os.str();
          }
          ++nPuSecSeen;
        }

        if (listing) {
          const bool hasOwn = own > 0, hasSec = sec + orph > 0;
          nListMixed += hasOwn && hasSec;
          nListOwnOnly += hasOwn && !hasSec;
          nListSecOnly += hasSec && !hasOwn;
          nListNeither += !hasOwn && !hasSec;
          // The sim hit with most charge, and the secondary or orphan one with most.
          std::pair<int, int> dom{-1, -1}, domSec{-1, -1};
          double domQ = -1, domSecQ = -1;
          for (const auto& [ci, q] : shCharge) {
            int ist;
            const SimHitClass k = classify((*coll[ci.first])[ci.second], ist);
            if (q > domQ) {
              domQ = q;
              dom = ci;
            }
            if (k != kShOwn && q > domSecQ) {
              domSecQ = q;
              domSec = ci;
            }
          }
          // The cluster's channel edges, in measurement units and in the local frame; the module's frame.
          const auto* det = dynamic_cast<const PixelGeomDetUnit*>(tkGeom->idToDetUnit(DetId(detId)));
          if (det == nullptr)
            throw cms::Exception("LogicError")
                << "MkFitTraceHitOrigin: DetId " << detId << " is no pixel-type det unit";
          const PixelTopology& topo = det->specificTopology();
          int r0, r1, c0, c1;
          if (pixel) {
            const SiPixelCluster& c = cl.pixelCluster();
            r0 = c.minPixelRow();
            r1 = c.maxPixelRow() + 1;
            c0 = c.minPixelCol();
            c1 = c.maxPixelCol() + 1;
          } else {
            const Phase2TrackerCluster1D& c = cl.phase2OTCluster();
            r0 = c.firstRow();
            r1 = c.firstRow() + c.size();
            c0 = c.column();
            c1 = c.column() + 1;
          }
          const float x0 = topo.localX(r0), x1 = topo.localX(r1), y0 = topo.localY(c0), y1 = topo.localY(c1);
          const Surface& s = det->surface();
          const GlobalPoint& o = s.position();
          const Surface::RotationType& R = s.rotation();
          const LocalPoint lp = rh->localPosition();
          const LocalError le = rh->localPositionError();
          const GlobalPoint gp = s.toGlobal(lp);
          auto clsOf = [&](const std::pair<int, int>& ci) {
            int ist;
            return ci.first < 0 ? -1 : (int)classify((*coll[ci.first])[ci.second], ist);
          };
          auto shOf = [&](const std::pair<int, int>& ci) -> const PSimHit* {
            return ci.first < 0 ? nullptr : &(*coll[ci.first])[ci.second];
          };
          const PSimHit* hd = shOf(dom);
          const PSimHit* hs = shOf(domSec);
          char row[1536];
          std::snprintf(row,
                        sizeof(row),
                        "%d,%d,%u,%d,%d,%zu,%.0f,%d,%d,%d,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,"
                        "%.5f,%.5f,%.5f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%.5f,%.5f,%.5f,"
                        "%d,%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%d,%zu,%lld,%d,%d,%d,%.5f,%d,%lld,%d,%d,%d,%.5f\n",
                        ilay,
                        i,
                        detId,
                        DetId(detId).subdetId(),
                        pixel ? 1 : 0,
                        chans.size(),
                        tot,
                        r0,
                        r1,
                        c0,
                        c1,
                        std::min(x0, x1),
                        std::max(x0, x1),
                        std::min(y0, y1),
                        std::max(y0, y1),
                        lp.x(),
                        lp.y(),
                        std::sqrt(le.xx()),
                        std::sqrt(le.yy()),
                        gp.x(),
                        gp.y(),
                        gp.z(),
                        o.x(),
                        o.y(),
                        o.z(),
                        R.xx(),
                        R.xy(),
                        R.xz(),
                        R.yx(),
                        R.yy(),
                        R.yz(),
                        R.zx(),
                        R.zy(),
                        R.zz(),
                        s.bounds().width(),
                        s.bounds().length(),
                        s.bounds().thickness(),
                        (int)cls,
                        hasOwn && hasSec ? 1 : 0,
                        fOwn,
                        fSec,
                        orph / tot,
                        pu / tot,
                        unm / tot,
                        noise / tot,
                        sim,
                        shCharge.size(),
                        hd ? (long long)hd->trackId() : -1LL,
                        hd ? hd->particleType() : 0,
                        hd ? (int)hd->processType() : -1,
                        clsOf(dom),
                        hd ? domQ / tot : 0.,
                        hd ? stPdg(hd->trackId()) : 0,
                        hs ? (long long)hs->trackId() : -1LL,
                        hs ? hs->particleType() : 0,
                        hs ? (int)hs->processType() : -1,
                        clsOf(domSec),
                        hs ? domSecQ / tot : 0.);
          rowsHits.emplace_back(row);
          for (const auto& [ci, q] : shCharge) {
            const PSimHit& h = (*coll[ci.first])[ci.second];
            int ist;
            const SimHitClass k = classify(h, ist);
            const Local3DPoint en = h.entryPoint(), ex = h.exitPoint();
            std::snprintf(row,
                          sizeof(row),
                          "%d,%d,%u,%s,%d,%u,%d,%d,%d,%.5f,%.3f,%.6f,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%d,%d\n",
                          ilay,
                          i,
                          detId,
                          simHitNames_[ci.first].c_str(),
                          ci.second,
                          h.trackId(),
                          h.particleType(),
                          (int)h.processType(),
                          (int)k,
                          q / tot,
                          h.energyLoss() * 1e6,
                          h.pabs(),
                          h.tof(),
                          en.x(),
                          en.y(),
                          en.z(),
                          ex.x(),
                          ex.y(),
                          ex.z(),
                          ist >= 0 ? simTracks[ist].type() : 0,
                          ist >= 0 ? creatorProc(ist) : -1);
            rowsSimHits.emplace_back(row);
          }
        }

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

    if (listing) {
      char row[512];
      for (int ist = 0; ist < (int)simTracks.size(); ++ist) {
        const SimTrack& t = simTracks[ist];
        const int v = t.vertIndex();
        const SimVertex* sv = v >= 0 && v < (int)simVertices.size() ? &simVertices[v] : nullptr;
        auto nsh = nSimHitsOfTrack.find(t.trackId());
        std::snprintf(row,
                      sizeof(row),
                      "%d,%u,%d,%g,%.6f,%.6f,%.6f,%.6f,%d,%.5f,%.5f,%.5f,%.5f,%d,%d,%d,%d\n",
                      ist,
                      t.trackId(),
                      t.type(),
                      t.charge(),
                      t.momentum().px(),
                      t.momentum().py(),
                      t.momentum().pz(),
                      t.momentum().e(),
                      v,
                      sv ? sv->position().x() : 0.,
                      sv ? sv->position().y() : 0.,
                      sv ? sv->position().z() : 0.,
                      sv ? sv->position().t() * 1e9 : 0.,
                      sv ? (int)sv->processType() : -1,
                      sv ? sv->parentIndex() : -1,
                      t.genpartIndex(),
                      nsh == nSimHitsOfTrack.end() ? 0 : nsh->second);
        rowsSimTracks.emplace_back(row);
      }
      char key[256];
      std::snprintf(key, sizeof(key), "%s,%u,%u,%llu,%d,", listTag_.c_str(), T.run, T.lumi, T.event, pass);
      {
        std::lock_guard<std::mutex> lock(listMutex_);
        for (const std::string& r : rowsHits)
          listHits_ << key << r;
        for (const std::string& r : rowsSimHits)
          listSimHits_ << key << r;
        for (const std::string& r : rowsSimTracks)
          listSimTracks_ << key << r;
        listHits_.flush();
        listSimHits_.flush();
        listSimTracks_.flush();
      }
      std::snprintf(row,
                    sizeof(row),
                    "MKFITHITLIST %u:%u:%llu %s tag=%s pass=%d hits=%zu simhits=%zu simtracks=%zu own_only=%lld "
                    "sec_only=%lld mixed=%lld neither=%lld\n",
                    T.run,
                    T.lumi,
                    T.event,
                    module_.c_str(),
                    listTag_.c_str(),
                    pass,
                    rowsHits.size(),
                    rowsSimHits.size(),
                    rowsSimTracks.size(),
                    nListOwnOnly,
                    nListSecOnly,
                    nListMixed,
                    nListNeither);
      std::cout << row;
    }

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

    // The pileup's lines.
    const std::string key =
        std::to_string(T.run) + ":" + std::to_string(T.lumi) + ":" + std::to_string(T.event) + " " + module_;
    std::snprintf(
        buf,
        sizeof(buf),
        "MKFITHITORIGIN-PU %s pileup=%d cls_all_own=%lld cls_all_sec=%lld cls_all_unk=%lld cls_all_none=%lld "
        "cls_all_sec_pu=%lld cls_all_own_pu=%lld q_pu_own=%.4f q_pu_sec=%.4f q_pu_orph=%.4f q_pu_unm=%.4f "
        "links_pu=%lld links_pu_unmatched=%lld event_links_pu=%lld found_pu_global=%lld "
        "found_pu_perModule=%lld pu_simhits=%lld pu_simtracks=%lld pu_simvertices=%lld pu_events=%zu "
        "pu_vert_bad=%lld parent_missing_sig=%lld/%lld parent_missing_pu=%lld/%lld sig_prefix_mismatch=%lld "
        "simhits_pu_own=%lld simhits_pu_sec=%lld simhits_pu_orph=%lld",
        key.c_str(),
        withPu ? 1 : 0,
        nClsAll[1],
        nClsAll[2],
        nClsAll[3],
        nClsAll[0],
        nClsAllSecPu,
        nClsAllOwnPu,
        share(qPuOwn),
        share(qPuSec),
        share(qPuOrph),
        share(qPuUnm),
        nLinksPu,
        nLinksPuUnm,
        nEvLinksPu,
        nFoundPu[0],
        nFoundPu[1],
        nPuSimHits,
        nPuSimTracks,
        nPuSimVertices,
        [&] {
          std::size_t n = 0;
          for (const auto& [b, evs] : puEvTracks)
            n += evs.size();
          return n;
        }(),
        nPuVertBad,
        nParentMissSig,
        nParentSig,
        nParentMissPu,
        nParentPu,
        nSigPrefixMismatch,
        nShClassPu[kShOwn],
        nShClassPu[kShSec],
        nShClassPu[kShOrph]);
    os << buf << "\n";
    // Pileup events per bunch crossing: event ids with SimTracks, the largest event number, event ids with sim hits.
    os << "MKFITHITORIGIN-PUEVENTS " << key;
    {
      std::set<int> bcrs;
      for (const auto& [b, evs] : puEvTracks)
        bcrs.insert(b);
      for (const auto& [b, evs] : puEvHits)
        bcrs.insert(b);
      for (int b : bcrs) {
        const auto& et = puEvTracks[b];
        os << " " << b << ":" << et.size() << "/" << (et.empty() ? 0 : *et.rbegin()) << "/" << puEvHits[b].size();
      }
    }
    os << "\n";
    // The playback records: per bunch crossing the number of pileup events of source 0, and whether a record equals
    // the first one given (sizes and event ids).
    const CrossingFramePlaybackInfoNew* pbFirst = nullptr;
    for (std::size_t k = 0; k < playbackTokens_.size(); ++k) {
      const CrossingFramePlaybackInfoNew& pb = iEvent.get(playbackTokens_[k]);
      os << "MKFITHITORIGIN-PLAYBACK " << key << " " << playbackTags_[k].encode() << " min_bunch=" << pb.minBunch_
         << " bunches=" << pb.nBcrossings_ << " sources=" << pb.maxNbSources_ << " events=" << pb.eventInfo_.size()
         << " per_bcr=";
      for (int b = 0; b < pb.nBcrossings_; ++b)
        os << (b ? "," : "") << (pb.minBunch_ + b) << ":" << pb.getNumberOfEvents(pb.minBunch_ + b, 0);
      if (pbFirst == nullptr) {
        pbFirst = &pb;
        os << " same_as_first=-";
      } else {
        bool same = pb.sizes_ == pbFirst->sizes_ && pb.eventInfo_.size() == pbFirst->eventInfo_.size();
        for (std::size_t e = 0; same && e < pb.eventInfo_.size(); ++e)
          same = pb.eventInfo_[e].eventID() == pbFirst->eventInfo_[e].eventID();
        os << " same_as_first=" << (same ? 1 : 0);
      }
      os << "\n";
    }
    // The pileup's secondary types on hits, and its unmatched links.
    st.assign(puSecTypes.begin(), puSecTypes.end());
    std::sort(st.begin(), st.end(), [](const auto& a, const auto& b) { return a.second.first > b.second.first; });
    qSecAll = 0;
    for (const auto& e : st)
      qSecAll += e.second.second;
    os << "MKFITHITORIGIN-PUSECTYPES " << key;
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
    os << "MKFITHITORIGIN-PUUNMATCHED " << key;
    for (int po = 0; po < 2; ++po) {
      std::snprintf(buf,
                    sizeof(buf),
                    " %s: no_simhit_of_event_on_module=%lld out_of_range=%lld other_module=%lld other_track=%lld;",
                    po ? "outer" : "pixel",
                    nPuUnm[po][0],
                    nPuUnm[po][1],
                    nPuUnm[po][2],
                    nPuUnm[po][3]);
      os << buf;
    }
    os << "\n";
    std::cout << os.str() << std::flush;
  }

}  // namespace mkfit_fwff

using mkfit_fwff::MkFitTraceHitOrigin;
DEFINE_FWK_MODULE(MkFitTraceHitOrigin);
