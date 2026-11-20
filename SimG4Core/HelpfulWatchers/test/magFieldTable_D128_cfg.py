# The magnetic field of the D128 simulation (era Phase2C26I13M9, auto:phase2_realistic_T35, as runG4Snitch_D128.py)
# sampled on cylindrical grids into magFieldTable_D128.root (MagFieldTableDumper). One event, no input.
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C26I13M9_cff import Phase2C26I13M9

process = cms.Process("FIELDTABLE", Phase2C26I13M9)
process.load("Configuration.Geometry.GeometryExtendedRun4D128Reco_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
GT = "auto:phase2_realistic_T35"
process.GlobalTag = GlobalTag(process.GlobalTag, GT, "")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

def grid(rMin, rMax, dr, zMin, zMax, dz, phiDeg):
    return cms.PSet(rMin = cms.double(rMin), rMax = cms.double(rMax), dr = cms.double(dr),
                    zMin = cms.double(zMin), zMax = cms.double(zMax), dz = cms.double(dz),
                    phiDeg = cms.vdouble(phiDeg))

process.dump = cms.EDAnalyzer("MagFieldTableDumper",
    fileName = cms.string("magFieldTable_D128.root"),
    label = cms.string("CMSSW_20_1_0_pre5, era Phase2C26I13M9, GlobalTag " + GT +
                       " (" + process.GlobalTag.globaltag.value() + "), Configuration.StandardSequences.MagneticField_cff"),
    grids = cms.VPSet(
        # 0: fine, inside the coil (tracker, MTD, HGCal)
        grid(0., 300., 1., -700., 700., 1., [0.]),
        # 1: coarse, out to the yoke
        grid(0., 800., 5., -1100., 1100., 5., [0.]),
        # 2: phi symmetry check, the yoke repeats every 30 degrees
        grid(0., 800., 10., -1100., 1100., 10., [0., 7.5, 15., 22.5, 30.]),
    )
)
process.p = cms.Path(process.dump)
