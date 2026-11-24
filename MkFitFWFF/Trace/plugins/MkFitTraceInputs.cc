// Package:    MkFitFWFF/Trace
// Class:      MkFitTraceInputs
//
// mkFit's inputs of one MkFitProducer's traced event, so that a trace file can be compared with another run without
// reconstructing them: the hits, the beam spot and the dead bins of the producer's EventOfHits, and the identity of
// mkFit's TrackerInfo (the geometry, with the material). It gives them to the MkFitTraceWriter service, which writes
// them into the trace file next to the records (TraceRNTupleWriter.h: Hits, Inputs, Geometry). Unlike
// MkFitTraceTruth it needs no simulation.
//
// - Hits: every hit of the two MkFitHitWrappers that has a layer, the Hit as mkFit has it (position, covariance,
//   detIDinLayer, charge per cm, spans) and its module's DetId. A hit's index is the one the records use, the index
//   of its cluster in the cluster collection (convertHits.h); the module checks that for every hit
//   (cluster_idx_mismatch below).
// - The beam spot, as MkFitEventOfHitsProducer set it into the EventOfHits (mkfit::BeamSpot).
// - The dead regions as mkFit uses them: the EventOfHits keeps only the bins they make dead
//   (LayerOfHits::isBinDead()), not the regions MkFitEventOfHitsProducer gave it, so the bins are written.
// - The geometry: the md5 of the file TrackerInfo::write_bin_file(geometryVersion) writes for the job's TrackerInfo,
//   computed again when TrackerRecoGeometryRecord changes; the file itself once per md5 (TraceRNTupleWriter's
//   Geometry). With geometryVersion as the dumper (RecoTracker/MkFit/test/DumpMkFitGeometry) was given, the md5 is the
//   md5 of the dumper's file.
//
// Prints one line per event:
//   MKFITINPUTS <run>:<lumi>:<event> <module> hits=<hits written> layers=<layers with hits>
//     cluster_idx_mismatch=<hits whose cluster index is not their index> eoh_mismatch=<layers whose EventOfHits hit
//     count or hit vector is not the wrapper's> n_hits_pixel=<n> n_hits_outer=<n> bs=<x>,<y>,<z>,<sigmaZ>,<width x>,
//     <width y>,<dxdz>,<dydz> dead_bins=<n> dead_layers=<layers with one> geom_md5=<md5> geom_version=<stamp>
//     geom_bytes=<n> geom_written=<0|1> pass=<the inputs' pass> ms=<time>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Digest.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"

#include "RecoTracker/MkFit/interface/MkFitClusterIndexToHit.h"
#include "RecoTracker/MkFit/interface/MkFitEventOfHits.h"
#include "RecoTracker/MkFit/interface/MkFitGeometry.h"
#include "RecoTracker/MkFit/interface/MkFitHitWrapper.h"
#include "RecoTracker/MkFit/interface/MkFitOutputWrapper.h"
#include "RecoTracker/MkFitCore/interface/Hit.h"
#include "RecoTracker/MkFitCore/interface/HitStructures.h"
#include "RecoTracker/MkFitCore/interface/TrackerInfo.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"

#include "MkFitFWFF/Trace/plugins/MkFitTraceWriter.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <string>
#include <unistd.h>
#include <vector>

namespace mkfit_fwff {

  class MkFitTraceInputs : public edm::global::EDAnalyzer<> {
  public:
    explicit MkFitTraceInputs(const edm::ParameterSet& iConfig);

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void analyze(edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const override;

    // The job's TrackerInfo as a geometry file: its bytes and their md5.
    struct Geom {
      std::string md5;
      std::vector<std::uint8_t> bin;
    };
    std::shared_ptr<const Geom> geometry(const edm::EventSetup& iSetup) const;

    const std::string module_;
    const std::string geometryVersion_;
    const edm::EDGetTokenT<MkFitOutputWrapper> mkFitOutputToken_;
    const edm::EDGetTokenT<MkFitHitWrapper> pixelHitsToken_;
    const edm::EDGetTokenT<MkFitHitWrapper> stripHitsToken_;
    const edm::EDGetTokenT<MkFitClusterIndexToHit> pixelClusterIndexToken_;
    const edm::EDGetTokenT<MkFitClusterIndexToHit> stripClusterIndexToken_;
    const edm::EDGetTokenT<std::vector<int>> pixelLayerToken_;
    const edm::EDGetTokenT<std::vector<int>> stripLayerToken_;
    const edm::EDGetTokenT<MkFitEventOfHits> eventOfHitsToken_;
    const edm::ESGetToken<MkFitGeometry, TrackerRecoGeometryRecord> mkFitGeomToken_;

    mutable std::mutex geomMutex_;
    mutable unsigned long long geomCacheId_ = 0;
    mutable std::shared_ptr<const Geom> geom_;
  };

  MkFitTraceInputs::MkFitTraceInputs(const edm::ParameterSet& iConfig)
      : module_(iConfig.getParameter<std::string>("mkFitProducer")),
        geometryVersion_(iConfig.getParameter<std::string>("geometryVersion")),
        mkFitOutputToken_(consumes(edm::InputTag(module_))),
        pixelHitsToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelHits"))),
        stripHitsToken_(consumes(iConfig.getParameter<edm::InputTag>("stripHits"))),
        pixelClusterIndexToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelHits"))),
        stripClusterIndexToken_(consumes(iConfig.getParameter<edm::InputTag>("stripHits"))),
        pixelLayerToken_(consumes(iConfig.getParameter<edm::InputTag>("pixelHits"))),
        stripLayerToken_(consumes(iConfig.getParameter<edm::InputTag>("stripHits"))),
        eventOfHitsToken_(consumes(iConfig.getParameter<edm::InputTag>("eventOfHits"))),
        mkFitGeomToken_(esConsumes()) {}

  void MkFitTraceInputs::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("mkFitProducer", "initialStepTrackCandidatesMkFit")
        ->setComment(
            "label of the traced MkFitProducer: the key of the inputs, and its output orders this module after it");
    desc.add<edm::InputTag>("pixelHits", edm::InputTag("mkFitSiPixelHits"));
    desc.add<edm::InputTag>("stripHits", edm::InputTag("mkFitSiPhase2Hits"));
    desc.add<edm::InputTag>("eventOfHits", edm::InputTag("mkFitEventOfHits"))
        ->setComment("the producer's EventOfHits: the beam spot and the dead bins");
    desc.add<std::string>("geometryVersion", "")
        ->setComment(
            "the stamp written into the geometry file whose md5 is taken; as DumpMkFitGeometry's geometryVersion");
    descriptions.add("mkFitTraceInputs", desc);
  }

  std::shared_ptr<const MkFitTraceInputs::Geom> MkFitTraceInputs::geometry(const edm::EventSetup& iSetup) const {
    const unsigned long long id = iSetup.get<TrackerRecoGeometryRecord>().cacheIdentifier();
    std::lock_guard<std::mutex> lock(geomMutex_);
    if (geom_ && id == geomCacheId_)
      return geom_;

    // TrackerInfo writes only into a file: a temporary one, read back and removed.
    std::string tmpl = (std::filesystem::temp_directory_path() / "mkfit-trace-geom-XXXXXX").string();
    const int fd = mkstemp(tmpl.data());
    if (fd < 0)
      throw cms::Exception("FileOpenError") << "MkFitTraceInputs: cannot create a temporary file " << tmpl;
    ::close(fd);
    auto g = std::make_shared<Geom>();
    try {
      iSetup.getData(mkFitGeomToken_).trackerInfo().write_bin_file(tmpl, geometryVersion_);
      std::ifstream in(tmpl, std::ios::binary);
      g->bin.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    } catch (...) {
      std::filesystem::remove(tmpl);
      throw;
    }
    std::filesystem::remove(tmpl);
    if (g->bin.empty())
      throw cms::Exception("FileReadError") << "MkFitTraceInputs: the geometry file " << tmpl << " is empty";
    cms::Digest dg;
    dg.append(reinterpret_cast<const char*>(g->bin.data()), g->bin.size());
    g->md5 = dg.digest().toString();
    geom_ = std::move(g);
    geomCacheId_ = id;
    return geom_;
  }

  void MkFitTraceInputs::analyze(edm::StreamID, const edm::Event& iEvent, const edm::EventSetup& iSetup) const {
    const auto t0 = std::chrono::steady_clock::now();
    edm::Service<MkFitTraceWriter> writer;
    if (!writer.isAvailable())
      throw cms::Exception("Configuration") << "MkFitTraceInputs needs the MkFitTraceWriter service";

    // Only to run after the producer.
    iEvent.get(mkFitOutputToken_);

    const auto& trackerInfo = iSetup.getData(mkFitGeomToken_).trackerInfo();
    const mkfit::EventOfHits& eoh = iEvent.get(eventOfHitsToken_).get();
    const MkFitHitWrapper& pixelHits = iEvent.get(pixelHitsToken_);
    const MkFitHitWrapper& stripHits = iEvent.get(stripHitsToken_);

    TraceInputTables T;
    T.run = iEvent.id().run();
    T.lumi = iEvent.id().luminosityBlock();
    T.event = iEvent.id().event();
    T.module = module_;
    const int nLayers = trackerInfo.n_layers();
    T.layers.resize(nLayers);
    for (int l = 0; l < nLayers; ++l)
      T.layers[l].layer = l;

    long long nHits = 0, nClusterMismatch = 0;
    auto doHits = [&](const MkFitHitWrapper& wrapper,
                      const MkFitClusterIndexToHit& clusterIndex,
                      const std::vector<int>& layerOf) {
      for (int i = 0; i < (int)layerOf.size(); ++i) {
        const int ilay = layerOf[i];
        if (ilay < 0)
          continue;
        ++nHits;
        auto& L = T.layers.at(ilay);
        L.n_hits = std::max(L.n_hits, i + 1);
        const TrackingRecHit* rh = clusterIndex.hits()[i];
        if (static_cast<const BaseTrackerRecHit*>(rh)->firstClusterRef().index() != (unsigned int)i)
          ++nClusterMismatch;
        const mkfit::Hit& h = wrapper.hits()[i];
        const float* e = h.errArray();
        L.idx.push_back(i);
        L.x.push_back(h.x());
        L.y.push_back(h.y());
        L.z.push_back(h.z());
        L.err.push_back({{e[0], e[1], e[2], e[3], e[4], e[5]}});
        L.detid.push_back(rh->geographicalId().rawId());
        L.det_in_layer.push_back(h.detIDinLayer());
        L.charge_pcm.push_back(h.chargePerCM());
        L.span_rows.push_back(h.spanRows());
        L.span_cols.push_back(h.spanCols());
      }
    };
    doHits(pixelHits, iEvent.get(pixelClusterIndexToken_), iEvent.get(pixelLayerToken_));
    doHits(stripHits, iEvent.get(stripClusterIndexToken_), iEvent.get(stripLayerToken_));
    T.n_hits_pixel = pixelHits.hits().size();
    T.n_hits_outer = stripHits.hits().size();

    // The EventOfHits must hold these hits: per layer the same count, and its hit vector the wrapper's.
    int nLayersWithHits = 0, nEohMismatch = 0;
    for (int l = 0; l < nLayers; ++l) {
      const mkfit::LayerOfHits& loh = eoh[l];
      const auto& L = T.layers[l];
      nLayersWithHits += !L.idx.empty();
      const mkfit::HitVec& hv = loh.is_pixel() ? pixelHits.hits() : stripHits.hits();
      if (loh.nHits() != L.idx.size() || (!L.idx.empty() && loh.hitArray() != hv.data()))
        ++nEohMismatch;
    }

    const mkfit::BeamSpot& bs = eoh.refBeamSpot();
    T.bs_x = bs.x;
    T.bs_y = bs.y;
    T.bs_z = bs.z;
    T.bs_sigma_z = bs.sigmaZ;
    T.bs_width_x = bs.beamWidthX;
    T.bs_width_y = bs.beamWidthY;
    T.bs_dxdz = bs.dxdz;
    T.bs_dydz = bs.dydz;

    using LoH = mkfit::LayerOfHits;
    T.n_phi_bins = LoH::axis_phi_t::c_N_end;
    int nDeadLayers = 0;
    for (int l = 0; l < nLayers; ++l) {
      const LoH& loh = eoh[l];
      const std::size_t before = T.dead_layer.size();
      for (unsigned int q = 0; q < loh.qNBins(); ++q)
        for (unsigned int p = 0; p < LoH::axis_phi_t::c_N_end; ++p)
          if (loh.isBinDead(LoH::bin_index_t(p), LoH::bin_index_t(q))) {
            T.dead_layer.push_back(l);
            T.dead_q_bin.push_back(q);
            T.dead_phi_bin.push_back(p);
          }
      nDeadLayers += T.dead_layer.size() > before;
    }

    const std::shared_ptr<const Geom> geom = geometry(iSetup);
    T.geom_md5 = geom->md5;
    T.geom_version = geometryVersion_;
    T.geom_bin = &geom->bin;

    bool geomWritten = false;
    const int pass = writer->writeInputs(T, &geomWritten);

    char buf[1024];
    std::snprintf(buf,
                  sizeof(buf),
                  "MKFITINPUTS %u:%u:%llu %s hits=%lld layers=%d cluster_idx_mismatch=%lld eoh_mismatch=%d "
                  "n_hits_pixel=%d n_hits_outer=%d bs=%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g dead_bins=%zu "
                  "dead_layers=%d geom_md5=%s geom_version=%s geom_bytes=%zu geom_written=%d pass=%d ms=%.1f",
                  T.run,
                  T.lumi,
                  T.event,
                  module_.c_str(),
                  nHits,
                  nLayersWithHits,
                  nClusterMismatch,
                  nEohMismatch,
                  T.n_hits_pixel,
                  T.n_hits_outer,
                  T.bs_x,
                  T.bs_y,
                  T.bs_z,
                  T.bs_sigma_z,
                  T.bs_width_x,
                  T.bs_width_y,
                  T.bs_dxdz,
                  T.bs_dydz,
                  T.dead_layer.size(),
                  nDeadLayers,
                  T.geom_md5.c_str(),
                  T.geom_version.c_str(),
                  geom->bin.size(),
                  geomWritten ? 1 : 0,
                  pass,
                  std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
    std::cout << buf << std::endl;
  }

}  // namespace mkfit_fwff

using mkfit_fwff::MkFitTraceInputs;
DEFINE_FWK_MODULE(MkFitTraceInputs);
