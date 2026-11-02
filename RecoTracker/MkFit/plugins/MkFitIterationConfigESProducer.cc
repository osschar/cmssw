#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"

#include "RecoTracker/MkFit/interface/MkFitGeometry.h"

// mkFit includes
#include "RecoTracker/MkFitCore/interface/IterationConfig.h"
#include "RecoTracker/MkFitCore/interface/V2p2Params.h"

#include <charconv>
#include <cstdlib>
#include <string>
#include <vector>

namespace {
  // The python PSet "v2p2" and its two score PSets. The defaults are those of
  // mkfit::V2p2Params, so an unset parameter leaves the C++ default.

  // A float default as the double with the same shortest decimal form, so that the
  // generated cfi shows 9.21 rather than 9.210000038146973; it converts back exactly.
  double dflt(float f) {
    char buf[32];
    const auto res = std::to_chars(buf, buf + sizeof(buf), f);
    return std::strtod(std::string(buf, res.ptr).c_str(), nullptr);
  }

  std::vector<double> hitEffGroupDefaults(const mkfit::V2p2ScoreParams &p) {
    std::vector<double> v;
    for (float f : p.hit_eff_grp)
      v.push_back(dflt(f));
    return v;
  }

  void fillScoreDescription(edm::ParameterSetDescription &d, const mkfit::V2p2ScoreParams &p) {
    d.add<double>("hitBonus", dflt(p.hit_bonus))->setComment("linear score: per hit taken");
    d.add<double>("overlapBonus", dflt(p.overlap_bonus))
        ->setComment("linear score: per hit beyond the first in a layer");
    d.add<double>("chi2Weight", dflt(p.chi2_weight))->setComment("linear score: per unit of chi2");
    d.add<double>("missPenalty", dflt(p.miss_penalty))->setComment("a real hole");
    d.add<double>("edgePenalty", dflt(p.edge_penalty))->setComment("crossing at the layer boundary");
    d.add<double>("gapPenalty", dflt(p.gap_penalty))->setComment("crossing in an inactive region");
    d.add<double>("stopPenalty", dflt(p.stop_penalty))->setComment("out of hole budget");
    d.add<double>("hitEff", dflt(p.hit_eff))->setComment("likelihood score: per-layer hit efficiency");
    d.add<std::vector<double>>("hitEffByGroup", hitEffGroupDefaults(p))
        ->setComment(
            "likelihood score: hit efficiency per layer group (pixel barrel, pixel endcap, strip barrel below layer "
            "10, strip barrel from layer 10, strip endcap); negative takes hitEff");
  }

  mkfit::V2p2ScoreParams readScore(const edm::ParameterSet &ps, const char *name) {
    mkfit::V2p2ScoreParams p;
    p.hit_bonus = ps.getParameter<double>("hitBonus");
    p.overlap_bonus = ps.getParameter<double>("overlapBonus");
    p.chi2_weight = ps.getParameter<double>("chi2Weight");
    p.miss_penalty = ps.getParameter<double>("missPenalty");
    p.edge_penalty = ps.getParameter<double>("edgePenalty");
    p.gap_penalty = ps.getParameter<double>("gapPenalty");
    p.stop_penalty = ps.getParameter<double>("stopPenalty");
    p.hit_eff = ps.getParameter<double>("hitEff");
    const auto grp = ps.getParameter<std::vector<double>>("hitEffByGroup");
    if (grp.size() != mkfit::V2p2ScoreParams::n_hit_eff_grp)
      throw cms::Exception("Configuration") << "v2p2." << name << ".hitEffByGroup needs "
                                            << mkfit::V2p2ScoreParams::n_hit_eff_grp << " values, got " << grp.size();
    for (int i = 0; i < mkfit::V2p2ScoreParams::n_hit_eff_grp; ++i)
      p.hit_eff_grp[i] = grp[i];
    return p;
  }

  mkfit::V2p2Params readV2p2(const edm::ParameterSet &ps) {
    mkfit::V2p2Params p;
    p.score_outward = readScore(ps.getParameter<edm::ParameterSet>("scoreOutward"), "scoreOutward");
    p.score_inward = readScore(ps.getParameter<edm::ParameterSet>("scoreInward"), "scoreInward");
    p.hit_chi2_cut = ps.getParameter<double>("hitChi2Cut");
    p.max_presel_hits = ps.getParameter<int>("maxPreselHits");
    p.max_sec_depth = ps.getParameter<int>("maxSecDepth");
    if (p.max_presel_hits < 1 || p.max_presel_hits > mkfit::V2p2Params::max_presel_hits_limit)
      throw cms::Exception("Configuration")
          << "v2p2.maxPreselHits must be in [1, " << mkfit::V2p2Params::max_presel_hits_limit << "], got "
          << p.max_presel_hits;
    if (p.max_sec_depth < 1 || p.max_sec_depth > mkfit::V2p2Params::max_sec_depth_limit)
      throw cms::Exception("Configuration") << "v2p2.maxSecDepth must be in [1, "
                                            << mkfit::V2p2Params::max_sec_depth_limit << "], got " << p.max_sec_depth;
    return p;
  }
}  // namespace

class MkFitIterationConfigESProducer : public edm::ESProducer {
public:
  MkFitIterationConfigESProducer(const edm::ParameterSet &iConfig);

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  std::unique_ptr<mkfit::IterationConfig> produce(const TrackerRecoGeometryRecord &iRecord);

private:
  const edm::ESGetToken<MkFitGeometry, TrackerRecoGeometryRecord> geomToken_;
  const std::string configFile_;
  const float minPtCut_;
  const unsigned int maxClusterSize_;
  const int backwardSearchMinPixelLayers_;
  const float backwardSearchPromptMaxD0_;
  const mkfit::V2p2Params v2p2_;
};

MkFitIterationConfigESProducer::MkFitIterationConfigESProducer(const edm::ParameterSet &iConfig)
    : geomToken_{setWhatProduced(this, iConfig.getParameter<std::string>("ComponentName")).consumes()},
      configFile_{iConfig.getParameter<edm::FileInPath>("config").fullPath()},
      minPtCut_{(float)iConfig.getParameter<double>("minPt")},
      maxClusterSize_{iConfig.getParameter<unsigned int>("maxClusterSize")},
      backwardSearchMinPixelLayers_{iConfig.getParameter<int>("backwardSearchMinPixelLayers")},
      backwardSearchPromptMaxD0_{(float)iConfig.getParameter<double>("backwardSearchPromptMaxD0")},
      v2p2_{readV2p2(iConfig.getParameter<edm::ParameterSet>("v2p2"))} {}

void MkFitIterationConfigESProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("ComponentName", "")->setComment("Product label");
  desc.add<edm::FileInPath>("config", edm::FileInPath())
      ->setComment("Path to the JSON file for the mkFit configuration parameters");
  desc.add<double>("minPt", 0.0)->setComment("min pT cut applied during track building");
  desc.add<unsigned int>("maxClusterSize", 8)->setComment("Max cluster size of SiStrip hits");
  desc.add<int>("backwardSearchMinPixelLayers", 0)
      ->setComment("If > 0, keep a backward-search extension only if its new hits lie on this many pixel layers");
  desc.add<double>("backwardSearchPromptMaxD0", 0.0)
      ->setComment("Candidates with |d0| to the beam spot below this (cm) need new hits on only one pixel layer");

  const mkfit::V2p2Params v2p2Defaults;
  edm::ParameterSetDescription scoreOutward, scoreInward, v2p2;
  fillScoreDescription(scoreOutward, v2p2Defaults.score_outward);
  fillScoreDescription(scoreInward, v2p2Defaults.score_inward);
  v2p2.add<edm::ParameterSetDescription>("scoreOutward", scoreOutward)->setComment("layer-step score, outward search");
  v2p2.add<edm::ParameterSetDescription>("scoreInward", scoreInward)->setComment("layer-step score, inward search");
  v2p2.add<double>("hitChi2Cut", dflt(v2p2Defaults.hit_chi2_cut))
      ->setComment("Kalman hit acceptance: chi2 against the scaled track and the hit covariance");
  v2p2.add<int>("maxPreselHits", v2p2Defaults.max_presel_hits)
      ->setComment("hits kept per candidate and sub-layer after the pre-selection");
  v2p2.add<int>("maxSecDepth", v2p2Defaults.max_sec_depth)->setComment("most hits one path may take in a layer");
  desc.add<edm::ParameterSetDescription>("v2p2", v2p2)->setComment("MkFinderV2p2 tunings (buildingRoutine 'v2p2')");

  descriptions.addWithDefaultLabel(desc);
}

std::unique_ptr<mkfit::IterationConfig> MkFitIterationConfigESProducer::produce(
    const TrackerRecoGeometryRecord &iRecord) {
  mkfit::ConfigJson cj;
  auto it_conf = cj.load_File(configFile_);
  it_conf->m_params.minPtCut = minPtCut_;
  it_conf->m_backward_params.minPtCut = minPtCut_;
  it_conf->m_params.maxClusterSize = maxClusterSize_;
  it_conf->m_backward_params.maxClusterSize = maxClusterSize_;
  it_conf->m_backward_search_min_pixel_layers = backwardSearchMinPixelLayers_;
  it_conf->m_backward_search_prompt_max_d0 = backwardSearchPromptMaxD0_;
  it_conf->m_v2p2 = v2p2_;
  it_conf->setupStandardFunctionsFromNames();
  return it_conf;
}

DEFINE_FWK_EVENTSETUP_MODULE(MkFitIterationConfigESProducer);
