"""Attach FWFFMkFitLooper and MkFitOutputDigest to an existing process.

The looper runs the whole input once per loop: loop 0 as configured, then one loop per entry of 'edits'. Give it an
input of one event (a pick-events skim, or PoolSource.eventsToProcess) so that every loop re-runs that event.

Example, the offline phase-2 initialStep with mkFit (trackingMkFitInitialStep), at the end of a cmsDriver config:

    import FWCore.ParameterSet.Config as cms
    from MkFitFWFF.Looper.customizeFWFFLooper import customizeFWFFLooper, fwffEdit
    process = customizeFWFFLooper(process,
        edits = [
            fwffEdit('same'),                                            # loop 1: nothing changed, no rebuild
            fwffEdit('rebuilt', 'initialStepTrackCandidatesMkFit'),        # loop 2: module rebuilt from its own PSet
            fwffEdit('noSeedCleaning', 'initialStepTrackCandidatesMkFit',  # loop 3: one switch flipped
                     seedCleaning = cms.bool(False)),
            fwffEdit('jsonRoundTrip', jsonPatch = '{}'),                   # loop 4: JSON rewritten unchanged, record reset
            fwffEdit('maxCands3', jsonPatch = '{"m_params": {"maxCandsPerSeed": 3}}'),  # loop 5: JSON edit
        ],
        jsonPath = 'RecoTracker/MkFit/data/mkfit-phase2-initialStep.json',
        candidates = 'initialStepTrackCandidates',
        tracks = 'initialStepTracks',
        tag = 'ttbar-ev123',
        event = '1:5:123')

Module edits ('overrides'): cms-typed values only, with the type and trackiness of the module's parameter
(cms.untracked.bool for 'limitConcurrency', for example); nested PSets are merged, VPSets replaced. An unknown name
or a type mismatch stops the job at the end of the loop before. Edits accumulate: loop 4 above still has
seedCleaning = False.

JSON edits ('jsonPatch'): mkFit's ConfigJsonPatcher syntax (RecoTracker/MkFitCore/src/IterationConfig.cc,
ConfigJsonPatcher::replace): a JSON object whose keys are member names or JSON pointers ("m_params/maxCandsPerSeed"),
with "[i..j]" keys for array ranges, e.g. '{"m_layer_configs": [{"[0..3]": {"m_select_max_dphi": 0.1}}]}'. Unknown
keys stop the job. The patches accumulate too. Fields that MkFitIterationConfigESProducer sets from its own PSet after
reading the file (minPtCut, maxClusterSize, the backward search, the low-pT duplicate cleaning, the v2p2 PSet) are
overwritten there, so patching them in the JSON has no effect. A non-empty jsonPatch ('{}' included) has the
framework reset TrackerRecoGeometryRecord and its dependent records before the next loop.

The JSON file: the looper rewrites the file at jsonPath, so it must be a scratch copy that
MkFitIterationConfigESProducer reads instead of the original. edm::FileInPath accepts only files under $CMSSW_BASE,
the release or the data area, and the looper refuses $CMSSW_BASE/src, the release, the data area and /cvmfs. So,
before cmsRun, copy the original afresh for every job and put its top first in the search path:

    S=$CMSSW_BASE/fwff-json/job1
    mkdir -p $S/RecoTracker/MkFit/data
    cp $CMSSW_BASE/src/RecoTracker/MkFit/data/mkfit-phase2-initialStep.json $S/RecoTracker/MkFit/data/
    export CMSSW_SEARCH_PATH=$S:$CMSSW_SEARCH_PATH

and pass the same relative path as the ESProducer's 'config' as jsonPath. The looper prints the resolved path and a
hash of the file on its 'construct' line.

Output: lines starting with FWFFLOOP (the looper's transitions, with t_ms since construction and dt_ms since the
previous FWFFLOOP line; the applied changes on the lines after changeModule.begin and jsonPatch.begin) and with
MKFITDIGEST (one per event: counts and order-independent hashes of the candidates and tracks, and the loop number).
Both go to std::cout.

The looper forces one stream, one concurrent lumi and one concurrent run; the number of threads stays. Each loop
reopens the input and every output: a PoolOutputModule writes out.root, out001.root, ...; drop outputs that are not
needed.

A module whose type the framework resolves at run time ('X@alpaka', e.g. LST's modules) can be rebuilt with the
framework of this branch: its ModuleRegistry::replaceModule resolves the type as the module's construction did. An
edit that changes the resolved type, e.g. the Alpaka backend through 'alpaka', is refused with an exception (an error
status in the interactive mode). The IB's framework crashes on such a module.

Interactive mode (controlDir set, edits empty): the looper is driven from outside, as a viewer would drive it.

    process = customizeFWFFLooper(process, controlDir = '/some/dir', jsonPath = '...', candidates = ..., tracks = ...)

After each loop the looper writes <controlDir>/status.json (written to status.json.tmp, then renamed) and waits for
the file <controlDir>/command, looking every pollMs milliseconds (a sleep between looks, no busy loop). It moves the
command to <controlDir>/commands/cmd-NNNN.py and reads it:
  - 'rerun' or 'stop' as the whole content (plain words);
  - otherwise Python, read with FWCore/PythonParameterSet's readPSetsFrom (a fresh interpreter per command), which
    must define exactly one of the PSets 'edit', 'rerun' or 'stop' at top level. 'edit' has the fields of an entry of
    'edits', each optional (name, moduleLabel, overrides, jsonPatch); fwffEdit() from this module makes one:

        import FWCore.ParameterSet.Config as cms
        from MkFitFWFF.Looper.customizeFWFFLooper import fwffEdit
        edit = fwffEdit('ce', 'hltInitialStepTrackCandidatesMkFit', buildingRoutine = cms.string('cloneEngine'))

An edit is applied as a whole or not at all. A command that cannot be read or applied (a Python error, an unknown
parameter, a type mismatch, a bad JSON patch, a module that cannot be rebuilt) leaves the job as it was: the looper
writes the error into status.json with state 'error' and waits for the next command. A client writes the command to
a temporary name in controlDir and renames it to 'command', then waits for the next status.json whose state is not
'running'. scripts/fwffLooperCtl.py does that and prints the result; it is also the reference for a viewer.

status.json fields:
  seq       increases with every write;  state  starting | running | waiting | error | stopped | ended
  waiting   true when the looper waits for a command (states waiting and error)
  loop      the loop just finished (waiting, error, stopped) or about to start (running)
  pid, t_ms (looper clock), epoch (wall clock, s), commands (number taken), jsonPath
  command   the last command taken: id, file, kind (edit, rerun, stop, invalid), name, text
  error     the reason a command was refused (state error), or why the looper stopped without a stop command
  applied   every command applied so far, in order: loop, command, kind, name, moduleLabel, jsonPatch, changes
  lastLoop  the loop just finished: loop, edit, events, digest (its MKFITDIGEST lines), times_ms with
            wait (status 'waiting' or 'error' to the command seen), apply (command seen to applied),
            next (applied to the loop's start), runLumi, event, end (event done to the end of the loop),
            loop (applied to the end of the loop), result (command seen to the end of the loop)
SIGUSR2 (cmsRun's graceful-stop signal) or idleTimeout end the wait with state 'stopped'.
"""

import FWCore.ParameterSet.Config as cms


def fwffEdit(name, moduleLabel='', jsonPatch='', **overrides):
    """One entry of FWFFMkFitLooper.edits: the edit applied at the end of a loop, before the next.

    name        printed on the FWFFLOOP lines of the loop this edit starts
    moduleLabel ED module to rebuild through ModuleChanger ('' for none)
    jsonPatch   ConfigJsonPatcher patch of the iteration JSON at the looper's jsonPath ('' for none)
    overrides   parameters of the module to replace, as cms-typed values
    """
    if overrides and not moduleLabel:
        raise ValueError("fwffEdit('%s'): overrides need a moduleLabel" % name)
    return cms.PSet(
        name = cms.string(name),
        moduleLabel = cms.string(moduleLabel),
        overrides = cms.PSet(**overrides),
        jsonPatch = cms.string(jsonPatch),
    )


def customizeFWFFLooper(process, edits=(), jsonPath='', candidates='initialStepTrackCandidates', tracks='', tag='',
                        event=None, controlDir='', pollMs=None, idleTimeout=None):
    """Adds the looper as process.looper and the digest on its own EndPath.

    edits      list of fwffEdit(...) PSets, one per loop after the first (scripted mode)
    jsonPath   iteration JSON the looper patches (see the module docstring); '' when no edit has a jsonPatch
    controlDir directory for the interactive mode (then edits must be empty); '' for the scripted mode
    pollMs     interactive mode: interval between looks for the command file (default 100)
    idleTimeout interactive mode: stop after this many seconds without a command (default 0: wait forever)
    candidates TrackCandidateCollection to digest ('' for none)
    tracks     reco::TrackCollection to digest ('' for none)
    tag        string printed on every MKFITDIGEST line
    event      'run:lumi:event' to keep only that event of a PoolSource (sets eventsToProcess); None leaves the
               source as it is
    """
    from MkFitFWFF.Looper.fwffMkFitLooper_cfi import fwffMkFitLooper
    from MkFitFWFF.Looper.mkFitOutputDigest_cfi import mkFitOutputDigest

    if controlDir and len(edits):
        raise ValueError('customizeFWFFLooper: edits and controlDir exclude each other')
    process.looper = fwffMkFitLooper.clone(
        edits = cms.VPSet(*edits),
        jsonPath = jsonPath,
        controlDir = controlDir,
    )
    if pollMs is not None:
        process.looper.pollMs = pollMs
    if idleTimeout is not None:
        process.looper.idleTimeout = float(idleTimeout)
    process.mkFitOutputDigest = mkFitOutputDigest.clone(
        candidates = candidates,
        tracks = tracks,
        tag = tag,
    )
    process.mkFitOutputDigestEndPath = cms.EndPath(process.mkFitOutputDigest)
    if process.schedule_() is not None:
        process.schedule_().append(process.mkFitOutputDigestEndPath)

    if event is not None:
        process.source.eventsToProcess = cms.untracked.VEventRange(event)

    # A looper forces these to 1 anyway; setting them avoids the framework's warning.
    process.options.numberOfStreams = 1
    process.options.numberOfConcurrentLuminosityBlocks = 1
    process.options.numberOfConcurrentRuns = 1
    return process
