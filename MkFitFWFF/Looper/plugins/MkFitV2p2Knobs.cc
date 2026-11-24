// -*- C++ -*-
//
// Package:    MkFitFWFF/Looper
// Class:      MkFitV2p2Knobs
//
// A service that sets mkFit's V2p2 knobs by name when the job starts, before any module runs: a fresh job with knobs
// other than the defaults, e.g. to compare with a looper loop that changed them. Needs a build with
// -DMKFIT_V2P2_MUTABLE; in the normal build a non-empty 'knobs' stops the job. See FWFFKnobs.h for the PSet. Prints
// one line per knob, starting with MKFITKNOB.

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "FWFFKnobs.h"

#include <iostream>
#include <sstream>

class MkFitV2p2Knobs {
public:
  MkFitV2p2Knobs(edm::ParameterSet const& iConfig, edm::ActivityRegistry&) {
    std::vector<fwff::KnobChange> changes;
    try {
      changes = fwff::prepareKnobs(iConfig.getParameter<edm::ParameterSet>("knobs"));
    } catch (cms::Exception& e) {
      throw cms::Exception("Configuration") << "MkFitV2p2Knobs: " << e.explainSelf();
    }
    std::ostringstream os;
    for (auto const& l : fwff::applyKnobs(changes))
      os << "MKFITKNOB service set " << l << "\n";
    os << "MKFITKNOB service mutable=" << (fwff::v2p2KnobsMutable() ? 1 : 0) << " set=" << changes.size() << "\n";
    std::cout << os.str() << std::flush;
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription knobs;
    knobs.setAllowAnything();
    edm::ParameterSetDescription desc;
    desc.add<edm::ParameterSetDescription>("knobs", knobs)
        ->setComment(
            "one PSet per knob group with the knobs to set, e.g. Window = cms.PSet(dphi_trk_fac = cms.double(1.5)); "
            "bool knobs take cms.bool, int knobs cms.int32, float knobs cms.double");
    descriptions.add("MkFitV2p2Knobs", desc);
  }
};

DEFINE_FWK_SERVICE(MkFitV2p2Knobs);
