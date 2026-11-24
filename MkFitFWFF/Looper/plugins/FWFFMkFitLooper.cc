// -*- C++ -*-
//
// Package:    MkFitFWFF/Looper
// Class:      FWFFMkFitLooper
//
// An EDLooper that runs the whole input once per loop and changes the configuration between loops, for re-running
// one event with edited mkFit configuration.
//
// Loop 0 runs the job as configured. At the end of each loop the looper takes the next edit and starts the next
// loop. An edit can
//  - rebuild an ED module (usually MkFitProducer) from its current PSet with 'overrides' copied over it, through
//    ModuleChanger (moduleLabel non-empty; empty overrides rebuild the module unchanged);
//  - patch the mkFit iteration JSON at 'jsonPath' with 'jsonPatch' (mkFit's ConfigJsonPatcher syntax) and have the
//    framework reset TrackerRecoGeometryRecord before the next loop, so that MkFitIterationConfigESProducer reads the
//    file again;
//  - do neither (re-run unchanged).
// Module and JSON edits are cumulative: each applies to the state left by the previous ones. An edit is applied as a
// whole or not at all: both parts are prepared and checked before either is committed.
//
// Two modes:
//  - scripted (controlDir empty): entry k of 'edits' is applied at the end of loop k; when the entries are used up,
//    the looper stops. A bad edit stops the job.
//  - interactive (controlDir set): after each loop the looper writes <controlDir>/status.json and waits for a file
//    <controlDir>/command, polling every pollMs. The command is moved to <controlDir>/commands/cmd-NNNN.py and read:
//    the plain words 'rerun' or 'stop', or Python read with FWCore/PythonParameterSet's readPSetsFrom that defines
//    exactly one of the PSets 'edit' (the fields of an entry of 'edits', e.g. made with fwffEdit()), 'rerun' or
//    'stop'. A command that cannot be read or applied is reported in status.json (state 'error') and the looper waits
//    again; the job goes on. See python/customizeFWFFLooper.py for the protocol and the status fields.
//
// Every transition prints one line starting with FWFFLOOP, with the loop number, the time since the looper was
// constructed and the time since the previous FWFFLOOP line (steady_clock, ms). The changes themselves are printed
// on the lines that follow the changeModule and jsonPatch lines.

#include "FWCore/Framework/interface/EDLooper.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetupRecordKey.h"
#include "FWCore/Framework/interface/LooperFactory.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/ModuleChanger.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/ScheduleInfo.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/Entry.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSetEntry.h"
#include "FWCore/ParameterSet/interface/VParameterSetEntry.h"
#include "FWCore/PythonParameterSet/interface/MakePyBind11ParameterSets.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/FileInPath.h"
#include "FWCore/Utilities/interface/UnixSignalHandlers.h"

#include "RecoTracker/MkFitCore/interface/IterationConfig.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"

#include "FWFFLoopCounter.h"

#include "nlohmann/json.hpp"

#include <unistd.h>

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace {
  using Clock = std::chrono::steady_clock;

  double msBetween(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
  }

  struct Edit {
    std::string name;
    std::string moduleLabel;
    edm::ParameterSet overrides;
    std::string jsonPatch;
  };

  // FNV-1a, to identify the contents of the JSON file in the log.
  uint64_t fnv1a(std::string const& s) {
    uint64_t h = 0xcbf29ce484222325ULL;
    for (unsigned char c : s) {
      h ^= c;
      h *= 0x100000001b3ULL;
    }
    return h;
  }

  std::string hex16(uint64_t v) {
    std::ostringstream os;
    os << std::hex << std::setw(16) << std::setfill('0') << v;
    return os.str();
  }

  std::string readFile(std::string const& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in)
      throw cms::Exception("FWFFMkFitLooper") << "cannot open " << path << " for reading";
    std::ostringstream os;
    os << in.rdbuf();
    return os.str();
  }

  void writeFile(std::string const& path, std::string const& text) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out << text;
    out.close();
    if (!out)
      throw cms::Exception("FWFFMkFitLooper") << "cannot write " << path;
  }

  std::string trim(std::string const& s) {
    auto const b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos)
      return std::string();
    auto const e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
  }

  bool pathIsUnder(std::filesystem::path const& p, std::filesystem::path const& top) {
    std::string const ps = p.string();
    std::string ts = top.string();
    while (ts.size() > 1 && ts.back() == '/')
      ts.pop_back();
    if (ps == ts)
      return true;
    return ps.size() > ts.size() && ps.compare(0, ts.size(), ts) == 0 && (ts == "/" || ps[ts.size()] == '/');
  }

  // Resolves jsonPath (absolute, or relative as an edm::FileInPath, which is how
  // MkFitIterationConfigESProducer resolves its 'config') and refuses a file in the source area, the release, the
  // data area or cvmfs. FileInPath accepts only files under CMSSW_BASE, CMSSW_RELEASE_BASE or CMSSW_DATA_PATH, so the
  // scratch copy goes under CMSSW_BASE but outside its src/, in a directory put first in CMSSW_SEARCH_PATH.
  std::string resolveJsonPath(std::string const& jsonPath) {
    namespace fs = std::filesystem;
    std::string const full = fs::path(jsonPath).is_absolute() ? jsonPath : edm::FileInPath(jsonPath).fullPath();
    fs::path const lexical = fs::absolute(full).lexically_normal();
    fs::path const canonical = fs::weakly_canonical(lexical);

    std::vector<fs::path> forbidden{"/cvmfs"};
    auto addEnv = [&forbidden](char const* var, char const* sub) {
      char const* v = std::getenv(var);
      if (v != nullptr && *v != '\0')
        forbidden.emplace_back(fs::path(v) / sub);
    };
    addEnv("CMSSW_BASE", "src");
    addEnv("CMSSW_RELEASE_BASE", "");
    addEnv("CMSSW_FULL_RELEASE_BASE", "");
    addEnv("CMSSW_DATA_PATH", "");
    for (auto const& f : forbidden) {
      fs::path const fl = fs::absolute(f).lexically_normal();
      fs::path const fc = fs::weakly_canonical(fl);
      if (pathIsUnder(lexical, fl) || pathIsUnder(canonical, fc) || pathIsUnder(lexical, fc) ||
          pathIsUnder(canonical, fl)) {
        throw cms::Exception("FWFFMkFitLooper")
            << "jsonPath '" << jsonPath << "' resolves to " << lexical.string() << " (" << canonical.string()
            << "), which is under " << fl.string()
            << ". The looper rewrites this file and writes only to a scratch copy: put the copy under "
               "$CMSSW_BASE/<scratch>/<same relative path>, outside $CMSSW_BASE/src, and prepend "
               "$CMSSW_BASE/<scratch> to CMSSW_SEARCH_PATH.";
      }
    }
    return lexical.string();
  }

  // Copies every parameter of 'over' into 'target'. Each must exist in 'target' with the same type and trackiness.
  // Nested PSets are merged recursively; plain parameters and VPSets are replaced. Each change is written to 'log'.
  void applyOverrides(edm::ParameterSet& target,
                      edm::ParameterSet const& over,
                      std::string const& where,
                      std::ostream& log) {
    for (auto const& name : over.getParameterNames()) {
      std::string const path = where + "." + name;
      if (auto const* oe = over.retrieveUnknown(name)) {
        auto const* te = target.retrieveUnknown(name);
        if (te == nullptr) {
          throw cms::Exception("FWFFMkFitLooper")
              << "override " << path << ": the module has no plain parameter of this name";
        }
        if (te->typeCode() != oe->typeCode() || te->isTracked() != oe->isTracked()) {
          throw cms::Exception("FWFFMkFitLooper")
              << "override " << path << ": type or trackiness differs; module has " << *te << ", override has " << *oe;
        }
        log << "FWFFLOOP   set " << path << ": " << *te << " -> " << *oe << "\n";
        target.insert(true, name, *oe);
      } else if (auto const* op = over.retrieveUnknownParameterSet(name)) {
        auto const* tp = target.retrieveUnknownParameterSet(name);
        if (tp == nullptr) {
          throw cms::Exception("FWFFMkFitLooper") << "override " << path << ": the module has no PSet of this name";
        }
        if (tp->isTracked() != op->isTracked()) {
          throw cms::Exception("FWFFMkFitLooper") << "override " << path << ": trackiness of the PSet differs";
        }
        edm::ParameterSet sub;
        sub.copyForModify(tp->pset());
        applyOverrides(sub, op->pset(), path, log);
        target.insertParameterSet(true, name, edm::ParameterSetEntry(sub, tp->isTracked()));
      } else if (auto const* ov = over.retrieveUnknownVParameterSet(name)) {
        auto const* tv = target.retrieveUnknownVParameterSet(name);
        if (tv == nullptr) {
          throw cms::Exception("FWFFMkFitLooper") << "override " << path << ": the module has no VPSet of this name";
        }
        if (tv->isTracked() != ov->isTracked()) {
          throw cms::Exception("FWFFMkFitLooper") << "override " << path << ": trackiness of the VPSet differs";
        }
        log << "FWFFLOOP   set " << path << ": VPSet of " << tv->size() << " replaced by VPSet of " << ov->size()
            << "\n";
        target.insertVParameterSet(true, name, *ov);
      }
    }
  }

  // The 'edit' PSet of an interactive command: the fields of an entry of 'edits', each optional, nothing else.
  Edit editFromCommandPSet(edm::ParameterSet const& p, std::string const& defaultName) {
    static std::set<std::string> const allowed{"name", "moduleLabel", "overrides", "jsonPatch"};
    for (auto const& n : p.getParameterNames()) {
      if (allowed.count(n) == 0) {
        throw cms::Exception("FWFFMkFitLooper")
            << "the edit has an unknown parameter '" << n << "'; allowed: name, moduleLabel, overrides, jsonPatch";
      }
    }
    auto str = [&p](char const* n, std::string const& def) -> std::string {
      if (!p.exists(n))
        return def;
      if (p.existsAs<std::string>(n, true))
        return p.getParameter<std::string>(n);
      if (p.existsAs<std::string>(n, false))
        return p.getUntrackedParameter<std::string>(n);
      throw cms::Exception("FWFFMkFitLooper") << "the edit's '" << n << "' must be a cms.string";
    };
    Edit e;
    e.name = str("name", defaultName);
    e.moduleLabel = str("moduleLabel", "");
    e.jsonPatch = str("jsonPatch", "");
    if (p.exists("overrides")) {
      if (p.existsAs<edm::ParameterSet>("overrides", true))
        e.overrides = p.getParameter<edm::ParameterSet>("overrides");
      else if (p.existsAs<edm::ParameterSet>("overrides", false))
        e.overrides = p.getUntrackedParameter<edm::ParameterSet>("overrides");
      else
        throw cms::Exception("FWFFMkFitLooper") << "the edit's 'overrides' must be a cms.PSet";
    }
    if (e.moduleLabel.empty() && !e.overrides.getParameterNames().empty())
      throw cms::Exception("FWFFMkFitLooper") << "the edit has overrides but no moduleLabel";
    return e;
  }

  std::string errorText(std::exception_ptr ep) {
    try {
      std::rethrow_exception(ep);
    } catch (cms::Exception const& e) {
      return e.explainSelf();
    } catch (std::exception const& e) {
      return std::string("std::exception: ") + e.what();
    } catch (...) {
      return "unknown exception";
    }
  }
}  // namespace

class FWFFMkFitLooper : public edm::EDLooper {
public:
  explicit FWFFMkFitLooper(edm::ParameterSet const& iConfig);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  std::set<edm::eventsetup::EventSetupRecordKey> modifyingRecords() const override;

private:
  enum class CommandKind { kEdit, kRerun, kStop };
  struct Command {
    CommandKind kind = CommandKind::kRerun;
    Edit edit;
  };

  void startingNewLoop(unsigned int iLoop) override;
  Status duringLoop(edm::Event const& iEvent, edm::EventSetup const&) override;
  Status endOfLoop(edm::EventSetup const&, unsigned int iLoop) override;
  void beginRun(edm::Run const& iRun, edm::EventSetup const&) override;
  void endRun(edm::Run const& iRun, edm::EventSetup const&) override;
  void beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) override;
  void endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) override;
  void endOfJob() override;

  void mark(char const* what, std::string const& extra = std::string()) const;

  // Applies an edit as a whole or not at all; throws when it cannot be applied. Returns the changes as text.
  std::vector<std::string> applyEdit(Edit const& edit);
  edm::ParameterSet prepareModuleEdit(Edit const& edit, std::vector<std::string>& changes);
  void commitModuleEdit(Edit const& edit, edm::ParameterSet const& next);

  // Interactive mode.
  Status interactiveEndOfLoop(unsigned int iLoop);
  std::string waitForCommand(std::string& reason);
  Command readCommand(std::string const& file, std::string const& text, unsigned int id) const;
  void fillLastLoop(unsigned int iLoop);
  void writeStatus(std::string const& state, int loop, std::string const& error = std::string());

  std::vector<Edit> edits_;
  std::string jsonPath_;
  edm::eventsetup::EventSetupRecordKey const recordKey_;

  // The current PSet of each edited module: ScheduleInfo keeps returning the original after a change.
  std::map<std::string, edm::ParameterSet> currentPSets_;
  // The current iteration JSON, as loaded from jsonPath_ at construction and patched since.
  std::unique_ptr<mkfit::ConfigJsonPatcher> jsonState_;

  int loop_ = -1;
  unsigned int nEvents_ = 0;
  bool resetNext_ = false;

  Clock::time_point const t0_;
  mutable Clock::time_point tLast_;

  // Interactive mode: the control directory, the poll interval, the idle limit, and the state for status.json.
  std::string controlDir_;
  unsigned int pollMs_ = 100;
  double idleTimeout_ = 0;
  unsigned long statusSeq_ = 0;
  unsigned int nCommands_ = 0;
  std::string loopEditName_ = "-";
  nlohmann::json applied_ = nlohmann::json::array();
  nlohmann::json lastCommand_;
  nlohmann::json lastLoop_;
  Clock::time_point tLoopStart_, tBeginLumi_, tEventDone_, tEndOfLoopBegin_;
  Clock::time_point tStatusWaiting_, tCommandSeen_, tApplied_;
  bool haveCommandTimes_ = false;
  double lastWaitMs_ = 0;
};

FWFFMkFitLooper::FWFFMkFitLooper(edm::ParameterSet const& iConfig)
    : recordKey_{edm::eventsetup::EventSetupRecordKey::makeKey<TrackerRecoGeometryRecord>()},
      t0_{Clock::now()},
      tLast_{t0_},
      controlDir_{iConfig.getParameter<std::string>("controlDir")},
      pollMs_{iConfig.getParameter<unsigned int>("pollMs")},
      idleTimeout_{iConfig.getParameter<double>("idleTimeout")} {
  bool anyJson = false;
  for (auto const& p : iConfig.getParameter<std::vector<edm::ParameterSet>>("edits")) {
    Edit e{p.getParameter<std::string>("name"),
           p.getParameter<std::string>("moduleLabel"),
           p.getParameter<edm::ParameterSet>("overrides"),
           p.getParameter<std::string>("jsonPatch")};
    if (e.moduleLabel.empty() && !e.overrides.getParameterNames().empty()) {
      throw cms::Exception("Configuration")
          << "FWFFMkFitLooper: edit " << edits_.size() << " ('" << e.name << "') has overrides but no moduleLabel";
    }
    anyJson = anyJson || !e.jsonPatch.empty();
    edits_.push_back(std::move(e));
  }
  bool const interactive = !controlDir_.empty();
  if (interactive && !edits_.empty()) {
    throw cms::Exception("Configuration")
        << "FWFFMkFitLooper: 'edits' and 'controlDir' exclude each other (scripted or interactive mode)";
  }
  if (interactive && pollMs_ == 0)
    pollMs_ = 1;

  std::string const jsonPath = iConfig.getParameter<std::string>("jsonPath");
  if (anyJson && jsonPath.empty()) {
    throw cms::Exception("Configuration") << "FWFFMkFitLooper: an edit has a jsonPatch but jsonPath is empty";
  }
  // In the interactive mode a JSON edit can come at any time: load the file whenever jsonPath is given.
  if (interactive && !jsonPath.empty())
    anyJson = true;

  std::ostringstream os;
  os << "edits=" << edits_.size();
  if (anyJson) {
    jsonPath_ = resolveJsonPath(jsonPath);
    // Load the file as MkFitIterationConfigESProducer does, then keep it in the patcher as JSON.
    std::string const text = readFile(jsonPath_);
    mkfit::ConfigJson cj;
    auto itConf = cj.load_File(jsonPath_);
    jsonState_ = std::make_unique<mkfit::ConfigJsonPatcher>(true);
    jsonState_->load(*itConf);
    os << " jsonPath=" << jsonPath_ << " jsonBytes=" << text.size() << " jsonHash=" << hex16(fnv1a(text));
  } else {
    os << " jsonPath=-";
  }
  if (interactive) {
    namespace fs = std::filesystem;
    fs::create_directories(fs::path(controlDir_) / "commands");
    controlDir_ = fs::absolute(controlDir_).lexically_normal().string();
    // A command left from an earlier job is not taken: it is moved aside.
    std::error_code ec;
    fs::path const stale = fs::path(controlDir_) / "command";
    if (fs::exists(stale, ec)) {
      fs::path const to = fs::path(controlDir_) / "commands" / ("stale-" + std::to_string(::getpid()) + ".py");
      fs::rename(stale, to, ec);
      os << " staleCommand=" << to.string();
    }
    os << " controlDir=" << controlDir_ << " pollMs=" << pollMs_ << " idleTimeout=" << idleTimeout_;
  }
  mark("construct", os.str());
  for (unsigned int i = 0; i < edits_.size(); ++i) {
    auto const& e = edits_[i];
    std::cout << "FWFFLOOP   edit " << i << " -> loop " << i + 1 << ": name='" << e.name << "' module='"
              << e.moduleLabel << "' overrides=" << e.overrides.getParameterNames().size() << " jsonPatch='"
              << e.jsonPatch << "'\n";
  }
  std::cout << std::flush;
  if (interactive)
    writeStatus("starting", -1);
}

void FWFFMkFitLooper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription overrides;
  overrides.setAllowAnything();

  edm::ParameterSetDescription edit;
  edit.add<std::string>("name", "")->setComment("a name for the loop this edit starts, printed on its FWFFLOOP lines");
  edit.add<std::string>("moduleLabel", "")
      ->setComment("ED module to rebuild through ModuleChanger; empty: no module change");
  edit.add<edm::ParameterSetDescription>("overrides", overrides)
      ->setComment(
          "parameters copied over the module's current PSet; each must exist there with the same type and "
          "trackiness; nested PSets are merged, VPSets replaced; empty: rebuild the module unchanged");
  edit.add<std::string>("jsonPatch", "")
      ->setComment(
          "patch in mkFit's ConfigJsonPatcher syntax (a JSON object; keys may be JSON pointers, '[i..j]' selects "
          "array ranges), applied to the iteration JSON at jsonPath; non-empty (also '{}') resets "
          "TrackerRecoGeometryRecord before the next loop");

  edm::ParameterSetDescription desc;
  desc.addVPSet("edits", edit, std::vector<edm::ParameterSet>())
      ->setComment("scripted mode: one entry per loop after the first; entry k is applied at the end of loop k");
  desc.add<std::string>("jsonPath", "")
      ->setComment(
          "the iteration JSON that MkFitIterationConfigESProducer reads: absolute, or relative as its FileInPath "
          "'config'; must be a scratch copy, not in $CMSSW_BASE/src, the release, the data area or cvmfs");
  desc.add<std::string>("controlDir", "")
      ->setComment(
          "interactive mode when non-empty (then 'edits' must be empty): after each loop the looper writes "
          "status.json here and waits for a file 'command'");
  desc.add<unsigned int>("pollMs", 100)->setComment("interactive mode: interval between looks for the command file");
  desc.add<double>("idleTimeout", 0.)
      ->setComment("interactive mode: stop when no command came for this many seconds; 0: wait forever");
  descriptions.add("fwffMkFitLooper", desc);
}

std::set<edm::eventsetup::EventSetupRecordKey> FWFFMkFitLooper::modifyingRecords() const {
  // EventProcessor::prepareForNextLoop asks for this right after endOfLoop, before the next loop starts, and resets
  // each record returned with its dependent records.
  if (!resetNext_)
    return {};
  mark("modifyingRecords", std::string("reset=") + recordKey_.name());
  return {recordKey_};
}

void FWFFMkFitLooper::mark(char const* what, std::string const& extra) const {
  auto const now = Clock::now();
  std::ostringstream os;
  os << "FWFFLOOP loop=" << loop_ << " " << what << " t_ms=" << std::fixed << std::setprecision(3)
     << msBetween(t0_, now) << " dt_ms=" << msBetween(tLast_, now);
  if (!extra.empty())
    os << " " << extra;
  os << "\n";
  std::cout << os.str() << std::flush;
  tLast_ = now;
}

void FWFFMkFitLooper::startingNewLoop(unsigned int iLoop) {
  loop_ = static_cast<int>(iLoop);
  fwff::currentLoop() = loop_;
  fwff::takeDigestLines();
  nEvents_ = 0;
  tLoopStart_ = Clock::now();
  tBeginLumi_ = Clock::time_point();
  std::string extra = "edit=";
  if (iLoop == 0)
    extra += "-";
  else if (controlDir_.empty())
    extra += "'" + edits_[iLoop - 1].name + "'";
  else
    extra += "'" + loopEditName_ + "'";
  mark("startingNewLoop", extra);
}

edm::EDLooper::Status FWFFMkFitLooper::duringLoop(edm::Event const& iEvent, edm::EventSetup const&) {
  ++nEvents_;
  tEventDone_ = Clock::now();
  std::ostringstream os;
  os << "event=" << iEvent.id().run() << ":" << iEvent.id().luminosityBlock() << ":" << iEvent.id().event();
  mark("eventDone", os.str());
  return kContinue;
}

void FWFFMkFitLooper::beginRun(edm::Run const& iRun, edm::EventSetup const&) {
  mark("beginRun", "run=" + std::to_string(iRun.run()));
}

void FWFFMkFitLooper::endRun(edm::Run const& iRun, edm::EventSetup const&) {
  mark("endRun", "run=" + std::to_string(iRun.run()));
}

void FWFFMkFitLooper::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
  if (tBeginLumi_ == Clock::time_point())
    tBeginLumi_ = Clock::now();
  mark("beginLumi", "lumi=" + std::to_string(iLumi.luminosityBlock()));
}

void FWFFMkFitLooper::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
  mark("endLumi", "lumi=" + std::to_string(iLumi.luminosityBlock()));
}

edm::EDLooper::Status FWFFMkFitLooper::endOfLoop(edm::EventSetup const&, unsigned int iLoop) {
  resetNext_ = false;
  tEndOfLoopBegin_ = Clock::now();
  mark("endOfLoop.begin", "events=" + std::to_string(nEvents_));
  if (!controlDir_.empty())
    return interactiveEndOfLoop(iLoop);

  if (iLoop >= edits_.size()) {
    mark("endOfLoop.end", "status=stop");
    return kStop;
  }
  Edit const& edit = edits_[iLoop];
  applyEdit(edit);
  mark("endOfLoop.end",
       "status=continue next='" + edit.name + "' resetRecord=" + std::string(resetNext_ ? recordKey_.name() : "-"));
  return kContinue;
}

edm::ParameterSet FWFFMkFitLooper::prepareModuleEdit(Edit const& edit, std::vector<std::string>& changes) {
  edm::ParameterSet const* current = nullptr;
  auto it = currentPSets_.find(edit.moduleLabel);
  if (it != currentPSets_.end()) {
    current = &it->second;
  } else {
    current = scheduleInfo()->parametersForModule(edit.moduleLabel);
    if (current == nullptr) {
      throw cms::Exception("FWFFMkFitLooper") << "no ED module with label '" << edit.moduleLabel << "' in the job";
    }
  }
  // ModuleHolderFactory::makeReplacementModule looks the maker up by the configured type as it is. A type resolved
  // at run time ("X@alpaka" becomes e.g. "alpaka_serial_sync::X") is not found there, the factory returns no
  // module, and ModuleRegistry::replaceModule dereferences it (measured: a segmentation fault when rebuilding
  // hltLST, LSTProducer@alpaka). Refuse such modules instead.
  std::string const type = current->getParameter<std::string>("@module_type");
  if (type.find('@') != std::string::npos) {
    throw cms::Exception("FWFFMkFitLooper")
        << "module '" << edit.moduleLabel << "' has the type '" << type
        << "', which the framework resolves at run time; ModuleChanger cannot rebuild such a module (the "
           "replacement factory looks up the unresolved type and returns nothing)";
  }

  std::ostringstream log;
  edm::ParameterSet next;
  next.copyForModify(*current);
  applyOverrides(next, edit.overrides, edit.moduleLabel, log);
  next.registerIt();
  std::istringstream lines(log.str());
  for (std::string l; std::getline(lines, l);) {
    std::string const prefix = "FWFFLOOP   ";
    changes.push_back(l.compare(0, prefix.size(), prefix) == 0 ? l.substr(prefix.size()) : l);
  }
  return next;
}

void FWFFMkFitLooper::commitModuleEdit(Edit const& edit, edm::ParameterSet const& next) {
  mark("changeModule.begin", "module=" + edit.moduleLabel);
  bool ok = false;
  try {
    ok = moduleChanger()->changeModule(edit.moduleLabel, next);
  } catch (cms::Exception& e) {
    e.addContext("FWFFMkFitLooper: changeModule('" + edit.moduleLabel + "') for the edit '" + edit.name + "'");
    throw;
  }
  if (!ok) {
    throw cms::Exception("FWFFMkFitLooper")
        << "changeModule('" << edit.moduleLabel << "') found no such ED module in the schedule";
  }
  currentPSets_[edit.moduleLabel] = next;
  std::ostringstream os;
  os << "module=" << edit.moduleLabel << " pset=" << next.id();
  mark("changeModule.end", os.str());
}

std::vector<std::string> FWFFMkFitLooper::applyEdit(Edit const& edit) {
  std::vector<std::string> changes;
  bool const doModule = !edit.moduleLabel.empty();
  bool const doJson = !edit.jsonPatch.empty();
  if (doJson && !jsonState_) {
    throw cms::Exception("FWFFMkFitLooper")
        << "the edit '" << edit.name << "' has a jsonPatch, but the looper has no jsonPath";
  }

  // 1. Prepare both parts; nothing in the job changes yet.
  edm::ParameterSet next;
  if (doModule)
    next = prepareModuleEdit(edit, changes);

  nlohmann::json backup;
  std::string text;
  int nReplaced = 0;
  if (doJson) {
    mark("jsonPatch.begin", "file=" + jsonPath_);
    jsonState_->cd_top();
    backup = jsonState_->get("");
    try {
      nlohmann::json const patch = nlohmann::json::parse(edit.jsonPatch);
      jsonState_->cd_top();
      // The patcher is verbose: it prints each changed value as "  <path>: old -> new".
      nReplaced = jsonState_->replace(patch);
      text = jsonState_->dump(2);
      text += "\n";
      // Check that the result loads as an IterationConfig, as MkFitIterationConfigESProducer will load it.
      std::string const check = jsonPath_ + ".check";
      writeFile(check, text);
      mkfit::ConfigJson cj;
      cj.load_File(check);
      std::remove(check.c_str());
    } catch (std::exception const& e) {
      jsonState_->cd_top();
      jsonState_->get("") = backup;
      throw cms::Exception("FWFFMkFitLooper")
          << "jsonPatch of the edit '" << edit.name << "' failed: " << e.what() << "\npatch: " << edit.jsonPatch;
    }
    changes.push_back("json " + edit.jsonPatch + ": replaced " + std::to_string(nReplaced));
  }

  // 2. Commit: the module first (its rebuild can still fail), then the JSON file.
  if (doModule) {
    for (auto const& c : changes) {
      if (c.compare(0, 4, "set ") == 0)
        std::cout << "FWFFLOOP   " << c << "\n";
    }
    std::cout << std::flush;
    try {
      commitModuleEdit(edit, next);
    } catch (...) {
      if (doJson) {
        jsonState_->cd_top();
        jsonState_->get("") = backup;
      }
      throw;
    }
  }
  if (doJson) {
    // Written in place, so that a symbolic link at jsonPath_ is followed rather than replaced.
    writeFile(jsonPath_, text);
    resetNext_ = true;
    mark("jsonPatch.end",
         "replaced=" + std::to_string(nReplaced) + " bytes=" + std::to_string(text.size()) +
             " jsonHash=" + hex16(fnv1a(text)));
  }
  return changes;
}

// ---------------------------------------------------------------------------------------------------------------------
// Interactive mode
// ---------------------------------------------------------------------------------------------------------------------

void FWFFMkFitLooper::fillLastLoop(unsigned int iLoop) {
  nlohmann::json times;
  if (iLoop == 0) {
    times["loop"] = msBetween(t0_, tEndOfLoopBegin_);
  } else if (haveCommandTimes_) {
    times["wait"] = lastWaitMs_;
    times["apply"] = msBetween(tCommandSeen_, tApplied_);
    times["next"] = msBetween(tApplied_, tLoopStart_);
    times["loop"] = msBetween(tApplied_, tEndOfLoopBegin_);
    times["result"] = msBetween(tCommandSeen_, tEndOfLoopBegin_);
  }
  if (tBeginLumi_ != Clock::time_point()) {
    times["runLumi"] = msBetween(tLoopStart_, tBeginLumi_);
    if (nEvents_ > 0)
      times["event"] = msBetween(tBeginLumi_, tEventDone_);
  }
  if (nEvents_ > 0)
    times["end"] = msBetween(tEventDone_, tEndOfLoopBegin_);
  lastLoop_ = nlohmann::json::object();
  lastLoop_["loop"] = iLoop;
  lastLoop_["edit"] = iLoop == 0 ? std::string("-") : loopEditName_;
  lastLoop_["events"] = nEvents_;
  lastLoop_["times_ms"] = times;
  lastLoop_["digest"] = fwff::takeDigestLines();
}

void FWFFMkFitLooper::writeStatus(std::string const& state, int loop, std::string const& error) {
  namespace fs = std::filesystem;
  auto const now = Clock::now();
  nlohmann::json s;
  s["seq"] = ++statusSeq_;
  s["state"] = state;
  s["waiting"] = (state == "waiting" || state == "error");
  s["loop"] = loop;
  s["pid"] = static_cast<long>(::getpid());
  s["t_ms"] = msBetween(t0_, now);
  s["epoch"] = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
  s["commands"] = nCommands_;
  s["jsonPath"] = jsonPath_;
  if (!lastCommand_.is_null())
    s["command"] = lastCommand_;
  if (!error.empty())
    s["error"] = error;
  s["applied"] = applied_;
  if (!lastLoop_.is_null())
    s["lastLoop"] = lastLoop_;
  std::string const path = (fs::path(controlDir_) / "status.json").string();
  std::string const tmp = path + ".tmp";
  try {
    writeFile(tmp, s.dump(1) + "\n");
    fs::rename(tmp, path);
  } catch (std::exception const& e) {
    std::cout << "FWFFLOOP status write failed: " << e.what() << "\n" << std::flush;
  }
  mark("status", "state=" + state + " seq=" + std::to_string(statusSeq_));
}

std::string FWFFMkFitLooper::waitForCommand(std::string& reason) {
  namespace fs = std::filesystem;
  fs::path const cmd = fs::path(controlDir_) / "command";
  auto const start = Clock::now();
  while (true) {
    std::error_code ec;
    if (fs::exists(cmd, ec)) {
      char name[32];
      std::snprintf(name, sizeof(name), "cmd-%04u.py", nCommands_ + 1);
      fs::path const claimed = fs::path(controlDir_) / "commands" / name;
      fs::rename(cmd, claimed, ec);
      if (!ec)
        return claimed.string();
    }
    if (edm::shutdown_flag.load()) {
      reason = "shutdown signal";
      return std::string();
    }
    if (idleTimeout_ > 0 && msBetween(start, Clock::now()) > 1000. * idleTimeout_) {
      reason = "idle timeout";
      return std::string();
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(pollMs_));
  }
}

FWFFMkFitLooper::Command FWFFMkFitLooper::readCommand(std::string const& file,
                                                      std::string const& text,
                                                      unsigned int id) const {
  Command c;
  std::string const t = trim(text);
  if (t == "rerun") {
    c.kind = CommandKind::kRerun;
    return c;
  }
  if (t == "stop") {
    c.kind = CommandKind::kStop;
    return c;
  }
  if (t.empty())
    throw cms::Exception("FWFFMkFitLooper") << "empty command";

  // readPSetsFrom starts a Python interpreter of its own and ends it before it returns; cmsRun's interpreter ended
  // after reading the configuration. A name ending in ".py" makes it read the file.
  std::unique_ptr<edm::ParameterSet> top = edm::cmspybind11::readPSetsFrom(file);
  std::vector<std::string> found;
  for (char const* n : {"edit", "rerun", "stop"}) {
    if (top->existsAs<edm::ParameterSet>(n, true) || top->existsAs<edm::ParameterSet>(n, false))
      found.emplace_back(n);
  }
  if (found.size() != 1) {
    std::ostringstream os;
    for (auto const& n : top->getParameterNames())
      os << " " << n;
    throw cms::Exception("FWFFMkFitLooper")
        << "the command must define exactly one of the PSets 'edit', 'rerun', 'stop'; its PSets:" << os.str();
  }
  if (found[0] == "rerun") {
    c.kind = CommandKind::kRerun;
  } else if (found[0] == "stop") {
    c.kind = CommandKind::kStop;
  } else {
    c.kind = CommandKind::kEdit;
    edm::ParameterSet const p = top->existsAs<edm::ParameterSet>("edit", true)
                                    ? top->getParameter<edm::ParameterSet>("edit")
                                    : top->getUntrackedParameter<edm::ParameterSet>("edit");
    c.edit = editFromCommandPSet(p, "cmd" + std::to_string(id));
  }
  return c;
}

edm::EDLooper::Status FWFFMkFitLooper::interactiveEndOfLoop(unsigned int iLoop) {
  fillLastLoop(iLoop);
  writeStatus("waiting", static_cast<int>(iLoop));
  tStatusWaiting_ = Clock::now();

  while (true) {
    std::string reason;
    std::string const file = waitForCommand(reason);
    auto const seen = Clock::now();
    if (file.empty()) {
      writeStatus("stopped", static_cast<int>(iLoop), "no command: " + reason);
      mark("endOfLoop.end", "status=stop reason='" + reason + "'");
      return kStop;
    }
    unsigned int const id = ++nCommands_;
    std::string text;
    lastCommand_ = nlohmann::json::object();
    lastCommand_["id"] = id;
    lastCommand_["file"] = file;
    mark("command.begin", "id=" + std::to_string(id) + " file=" + file);
    std::exception_ptr failure;
    try {
      text = readFile(file);
      lastCommand_["text"] = text.size() > 4000 ? text.substr(0, 4000) + "\n..." : text;
      Command const cmd = readCommand(file, text, id);
      if (cmd.kind == CommandKind::kStop) {
        lastCommand_["kind"] = "stop";
        writeStatus("stopped", static_cast<int>(iLoop));
        mark("endOfLoop.end", "status=stop reason='stop command'");
        return kStop;
      }
      std::vector<std::string> changes;
      std::string name = "rerun";
      if (cmd.kind == CommandKind::kEdit) {
        lastCommand_["kind"] = "edit";
        lastCommand_["name"] = cmd.edit.name;
        name = cmd.edit.name;
        changes = applyEdit(cmd.edit);
      } else {
        lastCommand_["kind"] = "rerun";
      }
      lastCommand_["name"] = name;
      tApplied_ = Clock::now();
      tCommandSeen_ = seen;
      lastWaitMs_ = msBetween(tStatusWaiting_, seen);
      haveCommandTimes_ = true;
      loopEditName_ = name;
      nlohmann::json a;
      a["loop"] = iLoop + 1;
      a["command"] = id;
      a["kind"] = lastCommand_["kind"];
      a["name"] = name;
      if (cmd.kind == CommandKind::kEdit) {
        a["moduleLabel"] = cmd.edit.moduleLabel;
        a["jsonPatch"] = cmd.edit.jsonPatch;
        a["changes"] = changes;
      }
      applied_.push_back(a);
      writeStatus("running", static_cast<int>(iLoop) + 1);
      mark("endOfLoop.end",
           "status=continue next='" + name + "' resetRecord=" + std::string(resetNext_ ? recordKey_.name() : "-"));
      return kContinue;
    } catch (...) {
      failure = std::current_exception();
    }
    std::string const error = errorText(failure);
    mark("command.error", "id=" + std::to_string(id));
    std::cout << "FWFFLOOP   " << error << "\n" << std::flush;
    if (!lastCommand_.contains("kind"))
      lastCommand_["kind"] = "invalid";
    // The looper waits again; the time of the next command counts from here.
    writeStatus("error", static_cast<int>(iLoop), error);
    tStatusWaiting_ = Clock::now();
  }
}

void FWFFMkFitLooper::endOfJob() {
  mark("endOfJob");
  if (!controlDir_.empty())
    writeStatus("ended", loop_);
}

DEFINE_FWK_LOOPER(FWFFMkFitLooper);
