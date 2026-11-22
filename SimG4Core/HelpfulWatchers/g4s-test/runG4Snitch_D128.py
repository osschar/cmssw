# MT/Claude 2026-11-20: CMSSW_20_1_0_pre5, D128. The cmsDriver base below, with the gun and the G4Snitch watcher
# of runG4Snitch.py (the D110 configuration of 15.1). The beamspot HGCALCloseBy was renamed CloseBy.
# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: CE_E_Front_120um_cfi -s GEN,SIM -n 10 --conditions auto:phase2_realistic_T35 --beamspot CloseBy --datatier GEN-SIM --eventcontent FEVTDEBUG --geometry ExtendedRun4D128 --era Phase2C26I13M9 --relval 9000,100 --fileout file:step1.root --python_filename base_D128_GENSIM.py --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C26I13M9_cff import Phase2C26I13M9

process = cms.Process('SIM',Phase2C26I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D128Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D128_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedCloseBy_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

import os
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(int(os.environ.get("G4SNITCH_NEVT", "20"))),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.obsolete.untracked.bool,
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('CE_E_Front_120um_cfi nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step1.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

# Generate at HGCal boundary or closer to origin.
HGC = True

process.generator = cms.EDProducer("CloseByParticleGunProducer",
    AddAntiParticle = cms.bool(False),
    PGunParameters = cms.PSet(
        ControlledByEta = cms.bool(False),
        Delta = cms.double(5),

        # Single particle
        # VarMin = cms.double(199.99),
        # VarMax = cms.double(200.01),
        # PartID = cms.vint32(2212),
        # NParticles = cms.int32(1),

        # Mix50
        # VarMin = cms.double(0.5),
        # VarMax = cms.double(20.00),
        # PartID = cms.vint32(11, -11, 12, -12, 22, 211, -211, 130, 310, 321, -321, 2212),
        # NParticles = cms.int32(500),

        # 2-particles
		VarMin = cms.double(20.00),
        VarMax = cms.double(50.00),
        PartID = cms.vint32(11, -11, 22, 111, 211, -211, 2112, 2212),
        NParticles = cms.int32(2),

        MaxVarSpread = cms.bool(False),
        Overlapping = cms.bool(True),
        Pointing = cms.bool(True),
        RandomShoot = cms.bool(False),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(1.7 if HGC else 0),
        MaxEta = cms.double(2.7 if HGC else 2.4),
        RMin = cms.double(54.99  if HGC else 0.01),
        RMax = cms.double(55.01  if HGC else 10.0),
        ZMin = cms.double(320.99 if HGC else -7.5),
        ZMax = cms.double(321.01 if HGC else  7.5),
        TMax = cms.double(0.05),
        TMin = cms.double(0.0),
        UseDeltaT = cms.bool(False),
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('random particles in phi and r windows')
)

# G4SNITCH_GUN=vertex replaces the gun at the HGCal face with particles from the origin (the CloseBy vertex smearing
# is zero): one mu- and one pi+ per event, pT 10-50 GeV, eta 1.7-2.7, through the tracker and the timing layer.
import os
if os.environ.get("G4SNITCH_GUN", "") == "vertex":
    process.generator = cms.EDProducer("FlatRandomPtGunProducer",
        PGunParameters = cms.PSet(
            MinPt = cms.double(10.0),
            MaxPt = cms.double(50.0),
            PartID = cms.vint32(13, 211),
            MinEta = cms.double(1.7),
            MaxEta = cms.double(2.7),
            MinPhi = cms.double(-3.14159265359),
            MaxPhi = cms.double(3.14159265359),
        ),
        Verbosity = cms.untracked.int32(0),
        psethack = cms.string('mu- and pi+ from the origin, pt 10-50, eta 1.7-2.7'),
        AddAntiParticle = cms.bool(False),
        firstRun = cms.untracked.uint32(1)
    )
# G4SNITCH_GUN=cms sends particles into most of CMS: mu-, e-, gamma, pi+ and n from the origin, pT 5-20 GeV, eta -5
# to 5, each with its anti-particle back to back (10 per event): trackers, MTD, ECAL barrel, HGCal, HCAL barrel and
# outer, HF, and the muon chambers (pT 5 GeV still reaches the barrel muon stations).
# G4SNITCH_GUN=single: one particle per event from the origin, G4SNITCH_PDG (default 211) at a fixed G4SNITCH_PT
# (GeV, default 10), eta uniform in G4SNITCH_ETAMIN..G4SNITCH_ETAMAX (default 1.6..2.9), phi uniform.
elif os.environ.get("G4SNITCH_GUN", "") == "single":
    pt = float(os.environ.get("G4SNITCH_PT", "10"))
    process.generator = cms.EDProducer("FlatRandomPtGunProducer",
        PGunParameters = cms.PSet(
            MinPt = cms.double(pt - 0.001),
            MaxPt = cms.double(pt + 0.001),
            PartID = cms.vint32(int(os.environ.get("G4SNITCH_PDG", "211"))),
            MinEta = cms.double(float(os.environ.get("G4SNITCH_ETAMIN", "1.6"))),
            MaxEta = cms.double(float(os.environ.get("G4SNITCH_ETAMAX", "2.9"))),
            MinPhi = cms.double(-3.14159265359),
            MaxPhi = cms.double(3.14159265359),
        ),
        Verbosity = cms.untracked.int32(0),
        psethack = cms.string('one particle from the origin'),
        AddAntiParticle = cms.bool(False),
        firstRun = cms.untracked.uint32(1)
    )
elif os.environ.get("G4SNITCH_GUN", "") == "cms":
    process.generator = cms.EDProducer("FlatRandomPtGunProducer",
        PGunParameters = cms.PSet(
            MinPt = cms.double(5.0),
            MaxPt = cms.double(20.0),
            PartID = cms.vint32(13, 11, 22, 211, 2112),
            MinEta = cms.double(-5.0),
            MaxEta = cms.double(5.0),
            MinPhi = cms.double(-3.14159265359),
            MaxPhi = cms.double(3.14159265359),
        ),
        Verbosity = cms.untracked.int32(0),
        psethack = cms.string('mu-, e-, gamma, pi+, n and anti-particles from the origin, pt 5-20, eta -5..5'),
        AddAntiParticle = cms.bool(True),
        firstRun = cms.untracked.uint32(1)
    )


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.generator)



# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion


### Added by MT -- TFile and SimWatcher

# process.TFileService = cms.Service("TFileService", 
#     fileName = cms.string("G4Snitch.root"),
#     closeFileFast = cms.untracked.bool(True)
#     )

# The G4Snitch watcher, configured from the environment (see SimG4Core/HelpfulWatchers/python/G4Snitch_cfi.py).
from SimG4Core.HelpfulWatchers.G4Snitch_cfi import g4snitch_watchers
process.g4SimHits.Watchers = g4snitch_watchers()
