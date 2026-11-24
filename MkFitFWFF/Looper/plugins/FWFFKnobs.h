#ifndef MkFitFWFF_Looper_FWFFKnobs_h
#define MkFitFWFF_Looper_FWFFKnobs_h

// Setting mkFit's V2p2 knobs (RecoTracker/MkFitCore/src/V2p2Config.h) by name from a PSet, for FWFFMkFitLooper and
// the service MkFitV2p2Knobs. Only a build with -DMKFIT_V2P2_MUTABLE has settable knobs; in the normal build they are
// constexpr and a non-empty knobs PSet is an error.
//
// The PSet has one PSet per knob group, holding the knobs to set, e.g.
//   knobs = cms.PSet(Window = cms.PSet(dphi_trk_fac = cms.double(1.5)), Policy = cms.PSet(use_wsr = cms.bool(False)))
// A bool knob takes a cms.bool, an int knob a cms.int32, a float knob a cms.double (tracked or untracked). An unknown
// group or knob, or another type, is an error. Prepare all changes first (checks only), then apply them.

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#if defined(MKFIT_V2P2_MUTABLE)
#include "RecoTracker/MkFitCore/interface/V2p2Knobs.h"
#endif

#include <string>
#include <vector>

namespace fwff {

  constexpr bool v2p2KnobsMutable() {
#if defined(MKFIT_V2P2_MUTABLE)
    return true;
#else
    return false;
#endif
  }

  struct KnobChange {
    std::string name;  // "<group>.<knob>"
#if defined(MKFIT_V2P2_MUTABLE)
    mkfit::Config::V2p2::Knob const* knob = nullptr;
    mkfit::Config::V2p2::Knob::Value value;
#endif
  };

  // Checks every knob of 'knobs' against the table; throws cms::Exception on the first problem. Changes nothing.
  inline std::vector<KnobChange> prepareKnobs(edm::ParameterSet const& knobs) {
    std::vector<KnobChange> out;
    auto const groups = knobs.getParameterNames();
    if (groups.empty())
      return out;
#if !defined(MKFIT_V2P2_MUTABLE)
    throw cms::Exception("FWFFKnobs") << "this build's V2p2 knobs are constexpr; setting them needs a build with "
                                         "USER_CXXFLAGS=-DMKFIT_V2P2_MUTABLE";
#else
    using mkfit::Config::V2p2::Knob;
    for (auto const& g : groups) {
      edm::ParameterSet const* gp = nullptr;
      if (knobs.existsAs<edm::ParameterSet>(g, true))
        gp = &knobs.getParameterSet(g);
      else if (knobs.existsAs<edm::ParameterSet>(g, false))
        gp = &knobs.getUntrackedParameterSet(g);
      else
        throw cms::Exception("FWFFKnobs")
            << "knobs: '" << g << "' must be a cms.PSet of the knobs of the group '" << g << "'";
      for (auto const& k : gp->getParameterNames()) {
        std::string const name = g + "." + k;
        Knob const* knob = mkfit::Config::V2p2::find_knob(name);
        if (knob == nullptr)
          throw cms::Exception("FWFFKnobs") << "knobs: no V2p2 knob '" << name << "'";
        KnobChange c{name, knob, Knob::Value{}};
        auto take = [&](auto tag, char const* cmsType) {
          using T = decltype(tag);
          if (gp->existsAs<T>(k, true))
            return T(gp->getParameter<T>(k));
          if (gp->existsAs<T>(k, false))
            return T(gp->getUntrackedParameter<T>(k));
          throw cms::Exception("FWFFKnobs") << "knobs: " << name << " is a " << knob->type_name()
                                            << " knob and takes a " << cmsType << ", got " << *gp->retrieveUnknown(k);
        };
        if (gp->retrieveUnknown(k) == nullptr)
          throw cms::Exception("FWFFKnobs") << "knobs: " << name << " must be a plain value, not a PSet or VPSet";
        switch (knob->address.index()) {
          case 0:
            c.value = take(bool(), "cms.bool");
            break;
          case 1:
            c.value = int(take(int(), "cms.int32"));
            break;
          default:
            c.value = float(take(double(), "cms.double"));
            break;
        }
        out.push_back(std::move(c));
      }
    }
    return out;
#endif
  }

  // Applies prepared changes; returns one line per knob: "<name> <type>: <old> -> <new>".
  inline std::vector<std::string> applyKnobs(std::vector<KnobChange> const& changes) {
    std::vector<std::string> lines;
#if defined(MKFIT_V2P2_MUTABLE)
    using namespace mkfit::Config::V2p2;
    for (auto const& c : changes) {
      auto const old = c.knob->get();
      c.knob->set(c.value);
      lines.push_back(c.name + " " + c.knob->type_name() + ": " + knob_value_string(old) + " -> " +
                      knob_value_string(c.value));
    }
#endif
    return lines;
  }

  // The current value of every knob, "<name>=<value>", in table order; empty in the normal build.
  inline std::vector<std::string> knobValues() {
    std::vector<std::string> out;
#if defined(MKFIT_V2P2_MUTABLE)
    using namespace mkfit::Config::V2p2;
    for (auto const& k : knobs())
      out.push_back(k.name + "=" + knob_value_string(k.get()));
#endif
    return out;
  }

}  // namespace fwff

#endif
