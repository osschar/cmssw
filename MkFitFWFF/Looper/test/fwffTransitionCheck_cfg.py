# fwffTransitionCheck_cfg.py: ModuleChanger and the begin/end transitions, framework only (no conditions, no input file).
# EmptySource with 2 events per loop; FWFFMkFitLooper rebuilds the two checker modules of
# MkFitFWFF/Looper/test/FWFFTransitionCheckers.cc. Each checker instance prints FWFFCHECK lines and, when destroyed,
# its counts and OK or UNMATCHED.
#   mode=scripted : 6 edits, each checker rebuilt 3 times (value 1, 2, 3), 7 loops
#   mode=fail     : interactive (controlDir=ctl); the driver sends a rebuild whose constructor throws, then
#                   rerun, then a good rebuild, then stop, e.g. with scripts/fwffLooperCtl.py:
#                     fwffLooperCtl.py ctl wait
#                     fwffLooperCtl.py ctl edit --name bad --module globalChecker --set "value=cms.int32(7)" \
#                                               --set "throwInConstructor=cms.bool(True)"   (refused, exit 1)
#                     fwffLooperCtl.py ctl rerun
#                     fwffLooperCtl.py ctl edit --name good --module globalChecker --set "value=cms.int32(8)"
#                     fwffLooperCtl.py ctl stop
# Expected with this branch's framework: every instance that ran is destroyed with beginJob=endJob=beginStream=
# endStream=1 and OK; the instance whose constructor threw, with all four 0.
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

o = VarParsing.VarParsing()
o.register('mode', 'scripted', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, '')
o.register('summary', 1, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.int, '')
o.parseArguments()

process = cms.Process('FWFFCHECK')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source('EmptySource')
process.maxEvents.input = 2
process.options.wantSummary = bool(o.summary)

process.globalChecker = cms.EDProducer('FWFFGlobalStreamCacheChecker', value = cms.int32(0),
                                       throwInConstructor = cms.bool(False))
process.streamChecker = cms.EDProducer('FWFFStreamChecker', value = cms.int32(0),
                                       throwInConstructor = cms.bool(False))
process.p = cms.Path(process.globalChecker + process.streamChecker)

from MkFitFWFF.Looper.customizeFWFFLooper import fwffEdit
if o.mode == 'scripted':
    edits = []
    for v in (1, 2, 3):
        edits.append(fwffEdit('g%d' % v, 'globalChecker', value = cms.int32(v)))
        edits.append(fwffEdit('s%d' % v, 'streamChecker', value = cms.int32(v)))
    process.looper = cms.Looper('FWFFMkFitLooper', edits = cms.VPSet(*edits))
elif o.mode == 'fail':
    process.looper = cms.Looper('FWFFMkFitLooper', edits = cms.VPSet(), controlDir = cms.string('ctl'),
                                idleTimeout = cms.double(120.))
else:
    raise ValueError('mode ' + o.mode)
