# rerunG4Snitch_D128.py: re-simulate events of a D128 GEN-SIM output (step1.root of runG4Snitch_D128.py, or any
# sample with the random-engine states saved by randomEngineStateProducer and the generatorSmeared HepMC product)
# with the G4Snitch watcher on. The engine states are restored, so Geant4 repeats the original simulation exactly.
#   G4SNITCH_INPUT=path/step1.root      the sample (required)
#   G4SNITCH_EVENTS=1:5[,1:9]           run:event ranges to re-simulate (default all)
#   G4SNITCH_LOOPS=N                    process them N times, with a looper (IterateNTimesLooper)
#   and the watcher's G4SNITCH_* settings (SimG4Core/HelpfulWatchers/python/G4Snitch_cfi.py), e.g. G4SNITCH_SELECT_ID.
# Output: rerun.root with only this process's g4SimHits products, and G4Snitch.root.
import os
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C26I13M9_cff import Phase2C26I13M9

process = cms.Process('SNITCH', Phase2C26I13M9)
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D128Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D128_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:' + os.environ['G4SNITCH_INPUT']))
if os.environ.get("G4SNITCH_EVENTS"):
    process.source.eventsToProcess = cms.untracked.VEventRange(*os.environ["G4SNITCH_EVENTS"].split(","))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
# G4SNITCH_NTHREADS=N: N framework threads (default 1); the snitch then writes one file per thread
process.options = cms.untracked.PSet(numberOfThreads = cms.untracked.uint32(int(os.environ.get('G4SNITCH_NTHREADS', '1'))),
                                     numberOfStreams = cms.untracked.uint32(0))

# Restore the engine states the original job saved at the start of each event.
process.RandomNumberGeneratorService.restoreStateLabel = cms.untracked.string('randomEngineStateProducer')

# G4SNITCH_LOOPS=N: process the input N times with the framework's IterateNTimesLooper (a looper runs one stream;
# each loop is a new run for Geant4 and the snitch).
if int(os.environ.get('G4SNITCH_LOOPS', '0')) > 0:
    process.looper = cms.Looper("IterateNTimesLooper", nTimes = cms.uint32(int(os.environ['G4SNITCH_LOOPS'])))

from SimG4Core.HelpfulWatchers.G4Snitch_cfi import g4snitch_watchers
process.g4SimHits.Watchers = g4snitch_watchers()

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:rerun.root'),
    outputCommands = cms.untracked.vstring('drop *', 'keep *_g4SimHits_*_SNITCH'))
process.simulation_step = cms.Path(process.psim)
process.out_step = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.simulation_step, process.out_step)
