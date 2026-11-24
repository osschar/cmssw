// -*- C++ -*-
//
// Package:     FWCore/Framework
// Class  :     edm::ModuleRegistry
//
// Implementation:
//     [Notes on implementation]
//
// Original Author:  Chris Jones
//         Created:  Fri, 23 Aug 2013 16:39:58 GMT
//

// system include files

// user include files
#include "FWCore/Framework/interface/ModuleRegistry.h"
#include "FWCore/Framework/interface/maker/ModuleMakerPluginFactory.h"
#include "FWCore/Framework/interface/resolveMaker.h"
#include "FWCore/Framework/src/ModuleHolderFactory.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/SignalSentry.h"
#include "FWCore/Utilities/interface/thread_safety_macros.h"

namespace {
  // Makers of replacement modules whose configured type is resolved at run time ("X@alpaka" becomes, e.g.,
  // "alpaka_serial_sync::X"). ModuleHolderFactory caches such a maker under the resolved type only, and its
  // makeReplacementModule looks up the configured type. Replacements happen between loops, on one thread.
  edm::ModuleHolderFactory::MakerMap& replacementMakers() {
    CMS_SA_ALLOW static edm::ModuleHolderFactory::MakerMap s_makers;
    return s_makers;
  }
}  // namespace

namespace edm {
  std::shared_ptr<maker::ModuleHolder> ModuleRegistry::getModule(
      MakeModuleParams const& p,
      std::string const& moduleLabel,
      signalslot::Signal<void(ModuleDescription const&)>& iPre,
      signalslot::Signal<void(ModuleDescription const&)>& iPost) {
    auto modItr = labelToModule_.find(moduleLabel);
    if (modItr == labelToModule_.end()) {
      auto modPtr = ModuleHolderFactory::get()->makeModule(p, typeResolverMaker_, iPre, iPost);

      if (maxModuleID_ < modPtr->moduleDescription().id()) {
        maxModuleID_ = modPtr->moduleDescription().id();
      }
      // Transfer ownership of worker to the registry
      labelToModule_[moduleLabel] = modPtr;
      return modPtr;
    }
    return get_underlying_safe(modItr->second);
  }

  std::shared_ptr<maker::ModuleHolder> ModuleRegistry::getExistingModule(std::string const& moduleLabel) {
    auto modItr = labelToModule_.find(moduleLabel);
    if (modItr == labelToModule_.end()) {
      return {};
    }
    return get_underlying_safe(modItr->second);
  }

  maker::ModuleHolder* ModuleRegistry::replaceModule(std::string const& iModuleLabel,
                                                     edm::ParameterSet const& iPSet,
                                                     edm::PreallocationConfiguration const& iPrealloc) {
    auto modItr = labelToModule_.find(iModuleLabel);
    if (modItr == labelToModule_.end()) {
      return nullptr;
    }

    ModuleDescription const& oldDescription = modItr->second->moduleDescription();
    std::string const& modtype = iPSet.getParameter<std::string>("@module_type");

    auto modPtr = ModuleHolderFactory::get()->makeReplacementModule(iPSet);
    if (not modPtr) {
      // The factory has no maker under the configured type, so the type is resolved at run time. Resolve it as the
      // construction of the module did (ModuleHolderFactory::findMaker), with the same resolver maker, and give the
      // maker the configuration as it is. The original configuration of the module must resolve to the same maker:
      // another C++ type, e.g. another Alpaka backend, would put products of other types than the ones registered
      // for this label.
      auto resolvedName = [this](ParameterSet const& pset) -> std::string {
        std::string const type = pset.getParameter<std::string>("@module_type");
        if (typeResolverMaker_ == nullptr) {
          return type;
        }
        auto resolver = typeResolverMaker_->makeResolver(pset);
        return resolver->resolveType(type, resolver->kInitialIndex).first;
      };
      auto const* maker =
          detail::resolveMaker<ModuleMakerPluginFactory>(modtype, typeResolverMaker_, iPSet, replacementMakers());
      ParameterSet const* oldPSet = pset::Registry::instance()->getMapped(oldDescription.parameterSetID());
      if (oldPSet == nullptr) {
        throw Exception(errors::LogicError)
            << "ModuleRegistry::replaceModule: the configuration of the module '" << iModuleLabel
            << "' is not in the ParameterSet registry, so its resolved type cannot be checked.";
      }
      auto const* oldMaker = detail::resolveMaker<ModuleMakerPluginFactory>(
          oldPSet->getParameter<std::string>("@module_type"), typeResolverMaker_, *oldPSet, replacementMakers());
      if (maker == nullptr or maker != oldMaker) {
        throw Exception(errors::Configuration)
            << "ModuleRegistry::replaceModule: the module '" << iModuleLabel << "' of type '" << modtype
            << "' resolves to '" << resolvedName(iPSet) << "' with the new configuration, but was made as '"
            << resolvedName(*oldPSet) << "'. A replacement must keep the module's C++ type (for an Alpaka module, "
            << "its backend).";
      }
      modPtr = maker->makeReplacementModule(iPSet);
    }
    if (not modPtr) {
      throw Exception(errors::LogicError) << "ModuleRegistry::replaceModule: no module was made for the label '"
                                          << iModuleLabel << "' of type '" << modtype << "'.";
    }
    modPtr->finishModuleInitialization(oldDescription, iPrealloc, nullptr);

    if (maxModuleID_ < modPtr->moduleDescription().id()) {
      maxModuleID_ = modPtr->moduleDescription().id();
    }

    // Transfer ownership of worker to the registry
    modItr->second = modPtr;
    return modItr->second.get();
  }

  void ModuleRegistry::deleteModule(std::string const& iModuleLabel,
                                    signalslot::Signal<void(ModuleDescription const&)>& iPre,
                                    signalslot::Signal<void(ModuleDescription const&)>& iPost) {
    auto modItr = labelToModule_.find(iModuleLabel);
    if (modItr == labelToModule_.end()) {
      throw cms::Exception("LogicError")
          << "Trying to delete module " << iModuleLabel
          << " but it does not exist in the ModuleRegistry. Please contact framework developers.";
    }
    // The pre signal must be given a reference to the ModuleDescription of the module before the module is deleted (the
    // signal is promised to be given the ModuleDescription in a stable location). After the module has been deleted
    // that stable location is gone, so the best we can do is to give a copy of the ModuleDescription.
    auto md = modItr->second->moduleDescription();
    auto guard = signalslot::make_sentry([&iPost, &md]() { iPost.emit(md); });
    iPre.emit(modItr->second->moduleDescription());
    labelToModule_.erase(modItr);
    guard.succeeded();
  }
}  // namespace edm
