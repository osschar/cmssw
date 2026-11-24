// -*- C++ -*-
//
// Package:    MkFitFWFF/Looper
// Class:      FWFFMkFitLooper
//
// An EDLooper that runs the whole input once per loop and changes the configuration between loops, for re-running
// one event with edited mkFit configuration.
//
// Loop 0 runs the job as configured. At the end of loop k, entry k of 'edits' is applied and loop k+1 starts; when
// the entries are used up, the looper stops. An entry can
//  - rebuild an ED module (usually MkFitProducer) from its current PSet with 'overrides' copied over it, through
//    ModuleChanger (moduleLabel non-empty; empty overrides rebuild the module unchanged);
//  - patch the mkFit iteration JSON at 'jsonPath' with 'jsonPatch' (mkFit's ConfigJsonPatcher syntax) and have the
//    framework reset TrackerRecoGeometryRecord before the next loop, so that MkFitIterationConfigESProducer reads the
//    file again;
//  - do neither (re-run unchanged).
// Module and JSON edits are cumulative: each applies to the state left by the previous ones.
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
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/FileInPath.h"

#include "RecoTracker/MkFitCore/interface/IterationConfig.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"

#include "FWFFLoopCounter.h"

#include "nlohmann/json.hpp"

#include <chrono>
#include <cstdint>
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
#include <vector>

namespace {
  using Clock = std::chrono::steady_clock;

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
}  // namespace

class FWFFMkFitLooper : public edm::EDLooper {
public:
  explicit FWFFMkFitLooper(edm::ParameterSet const& iConfig);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  std::set<edm::eventsetup::EventSetupRecordKey> modifyingRecords() const override;

private:
  void startingNewLoop(unsigned int iLoop) override;
  Status duringLoop(edm::Event const& iEvent, edm::EventSetup const&) override;
  Status endOfLoop(edm::EventSetup const&, unsigned int iLoop) override;
  void beginRun(edm::Run const& iRun, edm::EventSetup const&) override;
  void endRun(edm::Run const& iRun, edm::EventSetup const&) override;
  void beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) override;
  void endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) override;
  void endOfJob() override;

  void mark(char const* what, std::string const& extra = std::string()) const;
  void applyModuleEdit(Edit const& edit);
  void applyJsonEdit(Edit const& edit);

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
};

FWFFMkFitLooper::FWFFMkFitLooper(edm::ParameterSet const& iConfig)
    : recordKey_{edm::eventsetup::EventSetupRecordKey::makeKey<TrackerRecoGeometryRecord>()},
      t0_{Clock::now()},
      tLast_{t0_} {
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

  std::string const jsonPath = iConfig.getParameter<std::string>("jsonPath");
  if (anyJson && jsonPath.empty()) {
    throw cms::Exception("Configuration") << "FWFFMkFitLooper: an edit has a jsonPatch but jsonPath is empty";
  }

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
  mark("construct", os.str());
  for (unsigned int i = 0; i < edits_.size(); ++i) {
    auto const& e = edits_[i];
    std::cout << "FWFFLOOP   edit " << i << " -> loop " << i + 1 << ": name='" << e.name << "' module='"
              << e.moduleLabel << "' overrides=" << e.overrides.getParameterNames().size() << " jsonPatch='"
              << e.jsonPatch << "'\n";
  }
  std::cout << std::flush;
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
      ->setComment("one entry per loop after the first; entry k is applied at the end of loop k");
  desc.add<std::string>("jsonPath", "")
      ->setComment(
          "the iteration JSON that MkFitIterationConfigESProducer reads: absolute, or relative as its FileInPath "
          "'config'; must be a scratch copy, not in $CMSSW_BASE/src, the release, the data area or cvmfs");
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
  auto const ms = [](Clock::duration d) { return std::chrono::duration<double, std::milli>(d).count(); };
  std::ostringstream os;
  os << "FWFFLOOP loop=" << loop_ << " " << what << " t_ms=" << std::fixed << std::setprecision(3) << ms(now - t0_)
     << " dt_ms=" << ms(now - tLast_);
  if (!extra.empty())
    os << " " << extra;
  os << "\n";
  std::cout << os.str() << std::flush;
  tLast_ = now;
}

void FWFFMkFitLooper::startingNewLoop(unsigned int iLoop) {
  loop_ = static_cast<int>(iLoop);
  fwff::currentLoop() = loop_;
  nEvents_ = 0;
  std::string extra = "edit=";
  if (iLoop == 0)
    extra += "-";
  else
    extra += "'" + edits_[iLoop - 1].name + "'";
  mark("startingNewLoop", extra);
}

edm::EDLooper::Status FWFFMkFitLooper::duringLoop(edm::Event const& iEvent, edm::EventSetup const&) {
  ++nEvents_;
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
  mark("beginLumi", "lumi=" + std::to_string(iLumi.luminosityBlock()));
}

void FWFFMkFitLooper::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
  mark("endLumi", "lumi=" + std::to_string(iLumi.luminosityBlock()));
}

edm::EDLooper::Status FWFFMkFitLooper::endOfLoop(edm::EventSetup const&, unsigned int iLoop) {
  resetNext_ = false;
  mark("endOfLoop.begin", "events=" + std::to_string(nEvents_));
  if (iLoop >= edits_.size()) {
    mark("endOfLoop.end", "status=stop");
    return kStop;
  }
  Edit const& edit = edits_[iLoop];
  if (!edit.moduleLabel.empty())
    applyModuleEdit(edit);
  if (!edit.jsonPatch.empty()) {
    applyJsonEdit(edit);
    resetNext_ = true;
  }
  mark("endOfLoop.end",
       "status=continue next='" + edit.name + "' resetRecord=" + std::string(resetNext_ ? recordKey_.name() : "-"));
  return kContinue;
}

void FWFFMkFitLooper::applyModuleEdit(Edit const& edit) {
  auto it = currentPSets_.find(edit.moduleLabel);
  if (it == currentPSets_.end()) {
    edm::ParameterSet const* original = scheduleInfo()->parametersForModule(edit.moduleLabel);
    if (original == nullptr) {
      throw cms::Exception("FWFFMkFitLooper") << "no ED module with label '" << edit.moduleLabel << "' in the job";
    }
    it = currentPSets_.emplace(edit.moduleLabel, *original).first;
  }

  std::ostringstream log;
  edm::ParameterSet next;
  next.copyForModify(it->second);
  applyOverrides(next, edit.overrides, edit.moduleLabel, log);
  next.registerIt();
  std::cout << log.str() << std::flush;

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
  it->second = next;
  std::ostringstream os;
  os << "module=" << edit.moduleLabel << " pset=" << next.id();
  mark("changeModule.end", os.str());
}

void FWFFMkFitLooper::applyJsonEdit(Edit const& edit) {
  mark("jsonPatch.begin", "file=" + jsonPath_);
  int nReplaced = 0;
  std::string text;
  try {
    nlohmann::json const patch = nlohmann::json::parse(edit.jsonPatch);
    jsonState_->cd_top();
    // The patcher is verbose: it prints each changed value as "  <path>: old -> new".
    nReplaced = jsonState_->replace(patch);
    text = jsonState_->dump(2);
    text += "\n";
  } catch (std::exception const& e) {
    throw cms::Exception("FWFFMkFitLooper")
        << "jsonPatch of the edit '" << edit.name << "' failed: " << e.what() << "\npatch: " << edit.jsonPatch;
  }
  {
    // Written in place, so that a symbolic link at jsonPath_ is followed rather than replaced.
    std::ofstream out(jsonPath_, std::ios::binary | std::ios::trunc);
    out << text;
    out.close();
    if (!out)
      throw cms::Exception("FWFFMkFitLooper") << "cannot write " << jsonPath_;
  }
  try {
    // Read back as MkFitIterationConfigESProducer will.
    mkfit::ConfigJson cj;
    cj.load_File(jsonPath_);
  } catch (std::exception const& e) {
    throw cms::Exception("FWFFMkFitLooper")
        << "the patched JSON " << jsonPath_ << " does not load as an IterationConfig: " << e.what();
  }
  mark("jsonPatch.end",
       "replaced=" + std::to_string(nReplaced) + " bytes=" + std::to_string(text.size()) +
           " jsonHash=" + hex16(fnv1a(text)));
}

void FWFFMkFitLooper::endOfJob() { mark("endOfJob"); }

DEFINE_FWK_LOOPER(FWFFMkFitLooper);
