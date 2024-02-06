#include "SimG4Core/HelpfulWatchers/src/BeginOfTrackCounter.h"
#include "SimG4Core/HelpfulWatchers/src/SimTracer.h"
#include "SimG4Core/Watcher/interface/SimWatcherFactory.h"
// Adding a Watcher to collect G4step statistics:
#include "SimG4Core/HelpfulWatchers/interface/MonopoleSteppingAction.h"
#include "SimG4Core/HelpfulWatchers/src/G4StepStatistics.h"
#include "SimG4Core/HelpfulWatchers/interface/G4Snitch.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"

using namespace simwatcher;

DEFINE_SIMWATCHER(SimTracer);
DEFINE_SIMWATCHER(BeginOfTrackCounter);

// Adding a Watcher to collect G4step statistics:
DEFINE_SIMWATCHER(G4StepStatistics);

// Adding a Watcher to take care of steps of a monopole:
DEFINE_SIMWATCHER(MonopoleSteppingAction);

// Watcher for extracting detailed information for debugging, reco tunning, and visualization.
DEFINE_SIMWATCHER(G4Snitch);
