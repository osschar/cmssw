#ifndef RecoTracker_MkFit_MkFitIterationConfigRecord_h
#define RecoTracker_MkFit_MkFitIterationConfigRecord_h

#include "FWCore/Framework/interface/DependentRecordImplementation.h"
#include "FWCore/Utilities/interface/mplVector.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"

// The record of mkFit's iteration configurations (mkfit::IterationConfig, made by MkFitIterationConfigESProducer from
// a JSON file and its own parameters). It depends on TrackerRecoGeometryRecord, from which the producer consumes
// MkFitGeometry. A record of its own can be reset alone: a new iteration configuration then does not remake the
// tracker reco geometry, MkFitGeometry and LST's geometry and module data, which are in TrackerRecoGeometryRecord.
class MkFitIterationConfigRecord
    : public edm::eventsetup::DependentRecordImplementation<MkFitIterationConfigRecord,
                                                            edm::mpl::Vector<TrackerRecoGeometryRecord>> {};

#endif
