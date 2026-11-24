// -*- C++ -*-
//
// Package:     HelpfulWatchers
// Class  :     G4Snitch
//
/**\class SimTracer SimTracer.h SimG4Core/HelpfulWatchers/interface/G4Snitch.h

Description: Outputs a ROOT tree of kinematics.

*/
//
// Original Author: Matevz Tadel
//         Created: Mon Jan 22, 2024
//

// system include files
#include <cmath>
#include <cstdlib>
#include <cxxabi.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <typeinfo>
#include <unistd.h>

// user include files
#include "SimG4Core/HelpfulWatchers/interface/G4Snitch.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/EndOfTrack.h"

#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4Version.hh"
#include "G4Threading.hh"
#include "G4Event.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VTouchable.hh"
#include "G4Material.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "SimG4Core/SensitiveDetector/interface/SensitiveDetector.h"
#include "Randomize.hh"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <algorithm>
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"


#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

namespace {

  // The energy cuts that decide what is tracked and the last state are parameters: track_ekin_cut and
  // last_state_ekin_cut (G4Snitch::m_track_ekin_cut, m_last_state_ekin_cut).

  // Conversion to TTree units, all printouts are in G4 units
  constexpr double Pos_fac = 0.1;   // mm => cm
  constexpr double Mom_fac = 0.001; // MeV => GeV
  // time stays in ns

  void g4_pos_to_cms(const G4Track *t, G4S_Particle::Vec4D &x) {
    auto &pos = t->GetPosition();
    x.Set(Pos_fac * pos.x(), Pos_fac * pos.y(), Pos_fac * pos.z());
    x.fT = t->GetGlobalTime();
  }
  void g4_mom_to_cms(const G4Track *t, G4S_Particle::Vec4D &p) {
    auto  mom = t->GetMomentum();
    p.Set(Mom_fac * mom.x(), Mom_fac * mom.y(), Mom_fac * mom.z());
    p.fT = Mom_fac * t->GetTotalEnergy();
  }
  void g4_pos_mom_to_cms(const G4Track *t, G4S_Particle::Vec4D &x, G4S_Particle::Vec4D &p) {
    g4_pos_to_cms(t, x);
    g4_mom_to_cms(t, p);
  }
  void g4_pos_mom_to_cms_birth(const G4Track *t, G4S_Particle &part) {
    g4_pos_mom_to_cms(t, part.m_x_birth, part.m_p_birth);
  }
  void g4_pos_mom_to_cms_last(const G4Track *t, G4S_Particle &part) {
    g4_pos_mom_to_cms(t, part.m_x_last, part.m_p_last);
  }
  void g4_pos_mom_to_cms_death(const G4Track *t, G4S_Particle &part) {
    g4_pos_mom_to_cms(t, part.m_x_death, part.m_p_death);
  }

  void g4_pos_to_cms_begin(const G4StepPoint *sp, G4S_Step &step) {
    auto &pos = sp->GetPosition();
    step.m_x_beg.Set(Pos_fac * pos.x(), Pos_fac * pos.y(), Pos_fac * pos.z());
    step.m_x_beg.fT = sp->GetGlobalTime();
  }
  void g4_pos_to_cms_end(const G4StepPoint *sp, G4S_Step &step) {
    auto &pos = sp->GetPosition();
    step.m_x_end.Set(Pos_fac * pos.x(), Pos_fac * pos.y(), Pos_fac * pos.z());
    step.m_x_end.fT = sp->GetGlobalTime();
  }

  const char* b2yn(bool b) { return b ? "yes" : "no"; }
}

//------------------------------------------------------------------------------

G4Snitch::G4Snitch(const edm::ParameterSet &pSet)
      : m_verbose(pSet.getUntrackedParameter<bool>("verbose", false)),
        m_verbose_stack_level(pSet.getUntrackedParameter<bool>("verbose_stack_level", false)),
        m_verbose_transport(pSet.getUntrackedParameter<bool>("verbose_transport", false)),
        m_verbose_skip(pSet.getUntrackedParameter<bool>("verbose_skip", false)),
        m_verbose_skip_with_ids(pSet.getUntrackedParameter<bool>("verbose_skip_with_ids", false)),

        m_output_sensitive_steps(pSet.getUntrackedParameter<bool>("output_sensitive_steps", true)),
        m_output_inert_steps(pSet.getUntrackedParameter<bool>("output_inert_steps", true)),
        m_sensitive_step_ecut(pSet.getUntrackedParameter<double>("sensitive_step_ecut", 0.0)),
        m_output_sensitive_zero_steps(pSet.getUntrackedParameter<bool>("output_sensitive_zero_steps", false)),
        m_inert_step_ecut(pSet.getUntrackedParameter<double>("inert_step_ecut", 1e-5)),
        m_inert_aggregate_descendants(pSet.getUntrackedParameter<bool>("inert_aggregate_descendants", true)),
        m_sensitive_det_id(pSet.getUntrackedParameter<bool>("sensitive_det_id", true)),
        m_no_det_id_volumes(pSet.getUntrackedParameter<std::vector<std::string>>(
            "no_det_id_volumes", std::vector<std::string>{"HVQF", "OQUA"})),

        m_output_path(pSet.getUntrackedParameter<bool>("output_path", true)),
        m_path_max_angle(pSet.getUntrackedParameter<double>("path_max_angle", 0.05)),
        m_path_max_length(pSet.getUntrackedParameter<double>("path_max_length", 5.0)),
        m_output_path_bfield(pSet.getUntrackedParameter<bool>("output_path_bfield", false))
{
  m_path_cos_max_angle = std::cos(m_path_max_angle);
  m_select_track_id = pSet.getUntrackedParameter<int>("select_track_id", 0);
  m_select_cone = pSet.getUntrackedParameter<std::vector<double>>("select_cone", std::vector<double>());
  m_track_ekin_cut = pSet.getUntrackedParameter<double>("track_ekin_cut", 10.0);
  m_last_state_ekin_cut = pSet.getUntrackedParameter<double>("last_state_ekin_cut", 1.0);
  m_track_unstable_always = pSet.getUntrackedParameter<bool>("track_unstable_always", true);
  m_positron_by_total_energy = pSet.getUntrackedParameter<bool>("positron_by_total_energy", true);
  m_output_file = pSet.getUntrackedParameter<std::string>("output_file", "G4Snitch.root");
  m_output_per_thread = pSet.getUntrackedParameter<bool>("output_per_thread", false);
  // CMS runs a watcher on more than one thread only if it says it can (RunManagerMTWorker::createWatchers); the
  // snitch can when each thread writes its own file.
  setMT(m_output_per_thread);
  m_pset_dump = pSet.dump();
  // 2025-04-18 Linking through BuildFile.xml no longer works (-Wl,--as-needed)
  // so we just pull the dicts in here.
  gSystem->Load("libRootG4Snitch.so");
}

void G4Snitch::add_path_point(int id, const double *pos_mm, double t_ns, const double *mom_MeV, double e_MeV) {
  G4S_PathPoint pp;
  pp.m_x.Set(Pos_fac * pos_mm[0], Pos_fac * pos_mm[1], Pos_fac * pos_mm[2]);
  pp.m_x.fT = t_ns;
  pp.m_p.Set(Mom_fac * mom_MeV[0], Mom_fac * mom_MeV[1], Mom_fac * mom_MeV[2]);
  pp.m_p.fT = Mom_fac * e_MeV;
  pp.m_x0 = m_mat_x0;
  pp.m_lambda = m_mat_lambda;
  particle_steps(id).m_path.emplace_back(pp);
  if (m_output_path_bfield) {
    // The CMSSW field that Geant4's field (sim::Field) wraps, asked directly: sim::Field caches its last point, and
    // asking it here would change what the stepper gets back next, and so the simulation.
    float bx = 0, by = 0, bz = 0;
    if (m_bfield) {
      const GlobalVector v = m_bfield->inTesla(GlobalPoint(Pos_fac * pos_mm[0], Pos_fac * pos_mm[1], Pos_fac * pos_mm[2]));
      bx = v.x(); by = v.y(); bz = v.z();
    }
    (*m_pbfield_vec)[id].m_b.emplace_back(bx, by, bz);
  }
}

int G4Snitch::process_index(const G4VProcess *proc) {
  if (!proc)
    return -1;
  auto it = m_proc_index.find(proc);
  if (it != m_proc_index.end())
    return it->second;
  const int idx = (int)m_proc_names.size();
  m_proc_index.emplace(proc, idx);
  m_proc_names.push_back(proc->GetProcessName());
  m_proc_types.push_back((int)proc->GetProcessType());
  m_proc_subtypes.push_back(proc->GetProcessSubType());
  return idx;
}

int G4Snitch::sd_index(const G4VSensitiveDetector *sd) {
  if (!sd)
    return -1;
  auto it = m_sd_index.find(sd);
  if (it != m_sd_index.end())
    return it->second;
  const int idx = (int)m_sd_names.size();
  m_sd_index.emplace(sd, idx);
  m_sd_names.push_back(sd->GetName());
  int status = 0;
  char *dm = abi::__cxa_demangle(typeid(*sd).name(), nullptr, nullptr, &status);
  m_sd_classes.push_back(status == 0 && dm ? dm : typeid(*sd).name());
  std::free(dm);
  return idx;
}

int G4Snitch::volume_index(const G4StepPoint *sp) {
  const G4VTouchable *th = sp->GetTouchable();
  const int depth = th->GetHistoryDepth();
  uint64_t h = 1469598103934665603ULL;
  auto mix = [&h](uint64_t v) {
    h ^= v;
    h *= 1099511628211ULL;
  };
  for (int d = 0; d <= depth; ++d) {
    mix((uint64_t)(uintptr_t)th->GetVolume(d));
    mix((uint64_t)(uint32_t)th->GetCopyNumber(d));
  }
  auto it = m_vol_index.find(h);
  if (it != m_vol_index.end())
    return it->second;
  const int idx = (int)m_vol_paths.size();
  m_vol_index.emplace(h, idx);
  std::string path;
  bool no_det_id = false;
  for (int d = depth; d >= 0; --d) {
    if (!path.empty())
      path += '/';
    for (const std::string &n : m_no_det_id_volumes)
      if (th->GetVolume(d)->GetName() == n)
        no_det_id = true;
    path += th->GetVolume(d)->GetName();
    path += ':';
    path += std::to_string(th->GetCopyNumber(d));
  }
  m_vol_paths.push_back(path);
  m_vol_no_det_id.push_back(no_det_id);
  const G4LogicalVolume *lv = th->GetVolume(0)->GetLogicalVolume();
  m_vol_sd.push_back(sd_index(lv->GetSensitiveDetector()));
  const G4Material *mat = lv->GetMaterial();
  m_vol_materials.push_back(mat->GetName());
  m_vol_x0.push_back(Pos_fac * mat->GetRadlen());
  m_vol_lambda.push_back(Pos_fac * mat->GetNuclearInterLength());
  m_vol_density.push_back(mat->GetDensity() / (CLHEP::g / CLHEP::cm3));
  const G4ThreeVector &o = th->GetTranslation(0);
  m_vol_x.push_back(Pos_fac * o.x());
  m_vol_y.push_back(Pos_fac * o.y());
  m_vol_z.push_back(Pos_fac * o.z());
  return idx;
}

void G4Snitch::add_inert_step(int slot, const G4Step *step, double tot, double niel) {
  const G4StepPoint *sp1 = step->GetPreStepPoint(), *sp2 = step->GetPostStepPoint();
  const int vol = volume_index(sp1);
  const uint64_t key = ((uint64_t)(uint32_t)slot << 32) | (uint64_t)(uint32_t)vol;
  const G4ThreeVector &p1 = sp1->GetPosition(), &p2 = sp2->GetPosition();
  // A charged track's deposit is spread along its step: the midpoint, and the segment's own spread D D^T / 12. A
  // neutral track has no continuous loss: its deposit is local, at the step's end.
  const bool charged = step->GetTrack()->GetDefinition()->GetPDGCharge() != 0;
  double x[4], dd[3] = {0, 0, 0};
  if (charged) {
    x[0] = 0.5 * Pos_fac * (p1.x() + p2.x());
    x[1] = 0.5 * Pos_fac * (p1.y() + p2.y());
    x[2] = 0.5 * Pos_fac * (p1.z() + p2.z());
    x[3] = 0.5 * (sp1->GetGlobalTime() + sp2->GetGlobalTime());
    dd[0] = Pos_fac * (p2.x() - p1.x());
    dd[1] = Pos_fac * (p2.y() - p1.y());
    dd[2] = Pos_fac * (p2.z() - p1.z());
  } else {
    x[0] = Pos_fac * p2.x();
    x[1] = Pos_fac * p2.y();
    x[2] = Pos_fac * p2.z();
    x[3] = sp2->GetGlobalTime();
  }
  int ai;
  auto it = m_inert_acc_index.find(key);
  if (it == m_inert_acc_index.end()) {
    ai = (int)m_inert_acc.size();
    m_inert_acc_index.emplace(key, ai);
    InertAcc a{slot, vol, 0, 0, 0, {x[0], x[1], x[2], x[3]}, {0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, 0};
    m_inert_acc.push_back(a);
  } else {
    ai = it->second;
  }
  InertAcc &a = m_inert_acc[ai];
  ++a.n;
  a.e += tot;
  a.niel += niel;
  double d[4];
  for (int k = 0; k < 4; ++k) {
    d[k] = x[k] - a.ref[k];
    a.s1[k] += tot * d[k];
  }
  static const int ci[6] = {0, 1, 2, 0, 0, 1}, cj[6] = {0, 1, 2, 1, 2, 2};
  for (int c = 0; c < 6; ++c)
    a.s2[c] += tot * (d[ci[c]] * d[cj[c]] + dd[ci[c]] * dd[cj[c]] / 12.0);
  a.s2t += tot * d[3] * d[3];
}

void G4Snitch::flush_inert_sums() {
  // Every sum has e > 0: only steps with a deposit above zero are added.
  for (const InertAcc &a : m_inert_acc) {
    G4S_InertSum is;
    is.m_volume = a.vol;
    is.m_n_steps = a.n;
    is.m_edep.set_edep(a.e, a.niel);
    double mu[4];
    for (int k = 0; k < 4; ++k)
      mu[k] = a.s1[k] / a.e;
    is.m_x_mean.Set(a.ref[0] + mu[0], a.ref[1] + mu[1], a.ref[2] + mu[2]);
    is.m_x_mean.fT = a.ref[3] + mu[3];
    static const int ci[6] = {0, 1, 2, 0, 0, 1}, cj[6] = {0, 1, 2, 1, 2, 2};
    for (int c = 0; c < 6; ++c) {
      double v = a.s2[c] / a.e - mu[ci[c]] * mu[cj[c]];
      if (c < 3 && v < 0)
        v = 0;
      is.m_cov[c] = v;
    }
    is.m_t_rms = std::sqrt(std::max(0.0, a.s2t / a.e - mu[3] * mu[3]));
    particle_steps(a.slot).m_inert_sums.emplace_back(is);
  }
  m_inert_acc.clear();
  m_inert_acc_index.clear();
}

void G4Snitch::registerConsumes(edm::ConsumesCollector cc) {
  if (m_output_path_bfield)
    m_bfield_token = cc.esConsumes<MagneticField, IdealMagneticFieldRecord, edm::Transition::BeginRun>();
}

void G4Snitch::beginRun(edm::EventSetup const &es) {
  if (m_output_path_bfield)
    m_bfield = &es.getData(m_bfield_token);
}

void G4Snitch::steps_emplace_back() {
  m_psteps_vec->emplace_back(G4S_ParticleSteps());
  if (m_output_path_bfield)
    m_pbfield_vec->emplace_back(G4S_PathBField());
}
void G4Snitch::steps_resize(int n) {
  m_psteps_vec->resize(n);
  if (m_output_path_bfield)
    m_pbfield_vec->resize(n);
}
void G4Snitch::steps_reserve(int n) {
  m_psteps_vec->reserve(n);
  if (m_output_path_bfield)
    m_pbfield_vec->reserve(n);
}
void G4Snitch::steps_clear() {
  m_psteps_vec->clear();
  if (m_output_path_bfield)
    m_pbfield_vec->clear();
}

// this never gets called, sigh
G4Snitch::~G4Snitch() {
  std::cout <<"G4Snitch destructed\n";
}

void G4Snitch::open_file_tree() {
  m_thread_id = G4Threading::G4GetThreadId();
  m_output_name = m_output_file;
  if (m_output_per_thread) {
    const std::string tag = "_t" + std::to_string(m_thread_id);
    const auto dot = m_output_name.rfind(".root");
    m_output_name = dot == std::string::npos ? m_output_name + tag : m_output_name.insert(dot, tag);
  }
  m_file = TFile::Open(m_output_name.c_str(), "RECREATE");
  m_tree = new TTree("T", "G4 Dump Tree of some kind");
  // No ROOT implicit multi-threading for the snitch's trees. CMS runs Geant4 on hand-off threads outside TBB while the
  // TBB workers wait for them in OscarMTProducer::produce, so a basket compression handed to the TBB arena from
  // TTree::Fill never runs and every thread waits: with 4 threads the job hung at the first EndOfEvent (2026-11-23).
  m_tree->SetImplicitMT(false);
  m_tree->Branch("ev", &m_ev_number, "ev/I");

  m_info.reset(new G4S_Info);
  m_tree->Branch("i", m_info.get());

  m_vec_capacity = 16384;

  m_part_vec.reset(new std::vector<G4S_Particle>);
  m_part_vec->reserve(m_vec_capacity);
  m_tree->Branch("p", m_part_vec.get());

  m_psteps_vec.reset(new std::vector<G4S_ParticleSteps>);
  m_tree->Branch("s", m_psteps_vec.get());
  if (m_output_path_bfield) {
    m_pbfield_vec.reset(new std::vector<G4S_PathBField>);
    m_tree->Branch("b", m_pbfield_vec.get());
  }
  steps_reserve(m_vec_capacity);
}

void G4Snitch::write_tree_close_file() {
  m_file->cd();
  m_tree->Write();
  delete m_tree;
  // The process table that m_birth_proc and m_death_proc index.
  {
    TTree ptree("P", "G4Snitch process table: idx, name, type (G4ProcessType), subtype");
    ptree.SetImplicitMT(false);
    int idx = 0, type = 0, subtype = 0;
    std::string name;
    ptree.Branch("idx", &idx);
    ptree.Branch("name", &name);
    ptree.Branch("type", &type);
    ptree.Branch("subtype", &subtype);
    for (idx = 0; idx < (int)m_proc_names.size(); ++idx) {
      name = m_proc_names[idx];
      type = m_proc_types[idx];
      subtype = m_proc_subtypes[idx];
      ptree.Fill();
    }
    ptree.Write();
  }
  // The volume table that G4S_InertSum::m_volume indexes.
  {
    TTree vtree("V", "G4Snitch volume table: idx, path (name:copy from the world down), material, x0 (cm), density "
                     "(g/cm3), x, y, z (cm, the volume's origin in global coordinates), sd (index into tree D, -1 for "
                     "none)");
    vtree.SetImplicitMT(false);
    int idx = 0, sd = -1;
    std::string path, material;
    double x0 = 0, lambda = 0, density = 0, x = 0, y = 0, z = 0;
    vtree.Branch("idx", &idx);
    vtree.Branch("path", &path);
    vtree.Branch("material", &material);
    vtree.Branch("x0", &x0);
    vtree.Branch("lambda", &lambda);
    vtree.Branch("density", &density);
    vtree.Branch("x", &x);
    vtree.Branch("y", &y);
    vtree.Branch("z", &z);
    vtree.Branch("sd", &sd);
    for (idx = 0; idx < (int)m_vol_paths.size(); ++idx) {
      path = m_vol_paths[idx];
      material = m_vol_materials[idx];
      x0 = m_vol_x0[idx];
      lambda = m_vol_lambda[idx];
      density = m_vol_density[idx];
      x = m_vol_x[idx];
      y = m_vol_y[idx];
      z = m_vol_z[idx];
      sd = m_vol_sd[idx];
      vtree.Fill();
    }
    vtree.Write();
  }
  // The sensitive-detector table that the volume table's column "sd" indexes.
  {
    TTree dtree("D", "G4Snitch sensitive-detector table: idx, name, class");
    dtree.SetImplicitMT(false);
    int idx = 0;
    std::string name, cls;
    dtree.Branch("idx", &idx);
    dtree.Branch("name", &name);
    dtree.Branch("class", &cls);
    for (idx = 0; idx < (int)m_sd_names.size(); ++idx) {
      name = m_sd_names[idx];
      cls = m_sd_classes[idx];
      dtree.Fill();
    }
    dtree.Write();
  }
  write_config_tree();
  m_file->Close();
  delete m_file;
}

// The configuration tree "C" (key, value): what decides what this file holds. The parameters as used (defaults
// included), the rules of filter() in words, the units, the parameter set as given, the code (the stamps the
// Ceph area's build.sh and sync20.sh leave in $CMSSW_BASE: COMMIT.built at the last successful build, COMMIT at the
// last sync, each the commit and the dirty files), the G4SNITCH_* environment of the job, and the versions.
void G4Snitch::write_config_tree() {
  TTree ctree("C", "G4Snitch configuration: key, value");
  ctree.SetImplicitMT(false);
  std::string key, value;
  ctree.Branch("key", &key);
  ctree.Branch("value", &value);
  auto put = [&](const std::string &k, const std::string &v) {
    key = k;
    value = v;
    ctree.Fill();
  };
  auto num = [](double x) {
    std::ostringstream o;
    o.precision(10);
    o << x;
    return o.str();
  };
  auto vec = [&](const auto &v) {
    std::string r;
    for (const auto &x : v) {
      if (!r.empty())
        r += ",";
      if constexpr (std::is_same_v<std::decay_t<decltype(x)>, std::string>)
        r += x;
      else
        r += num(x);
    }
    return r;
  };
  auto file = [](const std::string &path) {
    std::ifstream f(path);
    if (!f)
      return std::string("(missing: ") + path + ")";
    std::ostringstream o;
    o << f.rdbuf();
    return o.str();
  };
  put("format", "trees T (events: ev, i, p, s, b), P (processes), V (volumes), D (sensitive detectors), C (this)");
  put("verbose", b2yn(m_verbose));
  put("verbose_stack_level", b2yn(m_verbose_stack_level));
  put("verbose_transport", b2yn(m_verbose_transport));
  put("verbose_skip", b2yn(m_verbose_skip));
  put("verbose_skip_with_ids", b2yn(m_verbose_skip_with_ids));
  put("output_sensitive_steps", b2yn(m_output_sensitive_steps));
  put("output_inert_steps", b2yn(m_output_inert_steps));
  put("sensitive_step_ecut", num(m_sensitive_step_ecut) + " MeV");
  put("output_sensitive_zero_steps", b2yn(m_output_sensitive_zero_steps));
  put("inert_step_ecut", num(m_inert_step_ecut) + " MeV");
  put("inert_aggregate_descendants", b2yn(m_inert_aggregate_descendants));
  put("sensitive_det_id", b2yn(m_sensitive_det_id));
  put("no_det_id_volumes", vec(m_no_det_id_volumes));
  put("output_path", b2yn(m_output_path));
  put("path_max_angle", num(m_path_max_angle) + " rad");
  put("path_max_length", num(m_path_max_length) + " cm");
  put("output_path_bfield", b2yn(m_output_path_bfield));
  put("select_track_id", std::to_string(m_select_track_id));
  put("select_cone", vec(m_select_cone));
  put("track_ekin_cut", num(m_track_ekin_cut) + " MeV");
  put("last_state_ekin_cut", num(m_last_state_ekin_cut) + " MeV");
  put("output_file", m_output_file);
  put("output_per_thread", b2yn(m_output_per_thread));
  put("output_name", m_output_name);
  put("thread_id", std::to_string(m_thread_id));
  put("track_unstable_always", b2yn(m_track_unstable_always));
  put("positron_by_total_energy", b2yn(m_positron_by_total_energy));
  put("rule.tracked",
      "a particle gets its own G4S_Particle if it is unstable and not a neutron (any energy; only with "
      "track_unstable_always), or a positron with total energy >= track_ekin_cut (only with "
      "positron_by_total_energy), or has kinetic energy >= track_ekin_cut at creation; otherwise its steps are filed "
      "under its nearest tracked ancestor (filter())");
  put("units", "positions cm, time ns, momenta and energies GeV, deposits MeV, field T, x0 and lambda cm");
  put("pset", m_pset_dump);
  const char *base = std::getenv("CMSSW_BASE");
  put("code.built", base ? file(std::string(base) + "/COMMIT.built") : "(no CMSSW_BASE)");
  put("code.synced", base ? file(std::string(base) + "/COMMIT") : "(no CMSSW_BASE)");
  const char *ver = std::getenv("CMSSW_VERSION");
  put("cmssw_version", ver ? ver : "");
  put("geant4_version", G4Version);
  char host[256] = {0};
  gethostname(host, sizeof(host) - 1);
  put("host", host);
  for (char **e = environ; e && *e; ++e) {
    const std::string kv(*e);
    if (kv.rfind("G4SNITCH_", 0) == 0) {
      const auto eq = kv.find('=');
      put("env." + kv.substr(0, eq), eq == std::string::npos ? "" : kv.substr(eq + 1));
    }
  }
  ctree.Write();
}

void G4Snitch::reset_output_structs() {
  m_part_vec->clear();
  steps_clear();
  m_vec_size = 0;
  m_gtp2vid.clear();
  m_stack.clear();
  m_stack_vid.clear();
  m_inert_acc.clear();
  m_inert_acc_index.clear();
  m_info->reset();
}

//------------------------------------------------------------------------------

bool G4Snitch::filter(const G4Track *t) const {
  // Returns true if the particle is to be filtered out.
  // ? filter neutrinos
  // ? filter low-E neutrons (considered NOT stable)
  const G4ParticleDefinition* pd = t->GetParticleDefinition();
  if (m_track_unstable_always && !pd->GetPDGStable() && pd->GetPDGEncoding() != 2112)
    return false;
  if (m_positron_by_total_energy && pd->GetPDGEncoding() == -11 && t->GetTotalEnergy() >= m_track_ekin_cut)
    return false;
  return t->GetKineticEnergy() < m_track_ekin_cut;
}

//------------------------------------------------------------------------------

void G4Snitch::update(const BeginOfJob* boj) {
  std::cout <<"++ signal BeginOfJob\n";
}

void G4Snitch::update(const BeginOfRun* bor) {
  std::cout << "++ signal BeginOfRun --- creating root file / tree\n";
  open_file_tree();
  reset_output_structs();
  m_total_tracks_accepted = m_total_tracks_skipped = 0;
}

void G4Snitch::update(const BeginOfEvent* boe) {
  std::cout << "++ signal BeginOfEvent\n";
  ++m_event;
  m_ev_number = (*boe)() ? (*boe)()->GetEventID() : -1;
  m_tracking = true;

  m_stack.clear();
  m_stack_vid.clear();
  m_starting_new_event = true; // delay initializatio until expected number of primaries are known
}

void G4Snitch::update(const BeginOfTrack* bot)
{
  const G4Track *iTrk = (*bot)();
  const G4ParticleDefinition* pd = iTrk->GetParticleDefinition();
  const int gid = iTrk->GetTrackID();
  const int pid = iTrk->GetParentID();

  // Every track, tracked or not, sums its material from its birth.
  m_mat_x0 = m_mat_lambda = 0;

  if (m_starting_new_event) {
    // G4Id of incoming track is the final primary id ... setup primary slots.
    // The first primary (processed last) will have id 1 (not 0).
    // They will be reversed into id order, ie, primary i --> slot i.
    // Slot 0 is filled with an "administrative mother" - this is a g4snitch thing.
    m_g4_num_primaries = gid;
    m_id_current_primary = 0;

    m_part_vec->resize(m_g4_num_primaries + 1);
    steps_resize(m_g4_num_primaries + 1);
    m_vec_size = m_g4_num_primaries + 1;
    m_info->m_n_primaries = m_g4_num_primaries; // XXX might not be true !!! fix at end !!!
    m_event_tracks_accepted = m_event_tracks_skipped = 0;

    m_primary_n_accepted_tracks = 0;

    m_last_g4_id = -1;

    /// insert 0, setup mamma particle
    m_stack.push_back(0);
    m_stack_vid.push_back(0);
    auto &p = particle(0);
    p.m_daughters_begin = 1;
    p.m_daughters_end = 1;
    p.m_g4_id = 0;
    p.m_g4_level = stack_level();

    m_starting_new_event = false;
  }

  m_last_g4_id = std::max(m_last_g4_id, gid); // XXX not needed? store to info ?

  // Figure out stack stuff
  int prev_level = stack_level();
  while (pid != stack_top()) {
    if (m_tracking == false && stack_top() == m_non_tracking_gid) {
      // Finish skipping !!!
      if (m_verbose_skip)
        printf("TRACKING RESUMED, skipped %d tracks, %d steps, daughters=%d\n",
               m_tracks_skipped, m_steps_skipped, m_daughters_skipped);
      m_event_tracks_skipped += m_tracks_skipped;
      m_tracking = true;
    }
    m_stack.pop_back();
    m_stack_vid.pop_back();
    if (m_verbose_stack_level && (m_tracking || m_verbose_skip_with_ids))
      printf("--- UP 1 level to %d\n", stack_level());
  }
  m_stack.push_back(gid);
  m_stack_vid.push_back(-1);
  if (m_verbose_stack_level && (m_tracking || m_verbose_skip_with_ids)) {
    printf("+++ DOWN 1 level to %d\n", stack_level());
  } else if (m_verbose && (m_tracking || m_verbose_skip_with_ids)) {
    int delta = stack_level() - prev_level;
    if (delta != 0)
      printf("%s %d level%s to %d\n", (delta > 0 ? "+++ DOWN" : "--- UP"),
             std::abs(delta), (std::abs(delta) > 1 ? "s" : ""), stack_level());
  }

  if (m_tracking == false) {
    if (m_verbose_skip && m_verbose_skip_with_ids) {
      printf("  skipping id=%d pid=%d, E-k=%.1f ...\n",
             gid, iTrk->GetParentID(), iTrk->GetKineticEnergy());
    }
    ++m_tracks_skipped;
    return;
  }

  if (gid > m_g4_num_primaries && filter(iTrk)) {
    if (m_verbose_skip)
      printf("TRACKING PAUSED at lvl=%d id=%d pid=%d pdg=%d stable=%s E_t=%.1f E_k=%.1f\n",
             stack_level(), gid, iTrk->GetParentID(), pd->GetPDGEncoding(), b2yn(pd->GetPDGStable()),
             iTrk->GetTotalEnergy(), iTrk->GetKineticEnergy());
    m_non_tracking_gid = gid;
    m_tracks_skipped = 1;
    m_steps_skipped = 0;
    m_daughters_skipped = 0;
    m_tracking = false;
    // The skipped subtree's steps go to the parent, the entry below the top: tracked, as tracking was on.
    m_id = m_stack_vid[m_stack_vid.size() - 2];
    assert(m_id >= 0);
    return;
  }

  if (m_verbose) {
    auto &pos = iTrk->GetPosition();
    auto  mom = iTrk->GetMomentum();
    printf("BEG lvl=%d id=%d pid=%d  pdg=%d stable=%s chg=%.1f m=%.1f pos=%.1f,%.1f,%.1f  mom=%.1f,%.1f,%.1f  t=%.3f E_kin=%.1f E_tot=%.1f v=%.4f\n",
      stack_level(), gid, iTrk->GetParentID(),
      pd->GetPDGEncoding(), b2yn(pd->GetPDGStable()), pd->GetPDGCharge(), pd->GetPDGMass(),
      pos.x(), pos.y(), pos.z(),
      mom.x(), mom.y(), mom.z(),
      iTrk->GetGlobalTime(), iTrk->GetKineticEnergy(), iTrk->GetTotalEnergy(), iTrk->GetVelocity()
    );
  }

  if (gid <= m_g4_num_primaries) {
    m_g4_id_current_primary = gid;
    ++m_id_current_primary;

    // no need to schlep finished primaries and their daughters along.
    m_gtp2vid.clear();

    m_id = m_id_current_primary;
    // Primaries require extra initialization otherwise done in step processing
    auto &p = particle(m_id);
    p.m_parent = 0;
    p.m_g4_level = stack_level();
    ++particle(0).m_daughters_end;

    m_primary_n_accepted_tracks = 1;
  } else {
    // For secondaries, there should be an assigned location in the output vector.
    auto tpmi = m_gtp2vid.find(iTrk);
    assert(tpmi != m_gtp2vid.end());
    m_id = tpmi->second;
    m_gtp2vid.erase(tpmi);

    // However, there might be mechanisms that produce "filter-passing-particles" for whoose mothers
    // we do not account for properly in the filter() function, i.e.:
    // - unstable low E_k particles that will decay into high-enough E_k ones;
    // - positrons
    // - deep inelastic scattering (no secondary, particle is just replaced with a new one).
  }
  m_stack_vid.back() = m_id;
  auto &p = particle(m_id);
  p.m_daughters_begin = p.m_daughters_end = m_vec_size;
  p.m_g4_id = gid;
  p.m_was_tracked = true;

  m_num_accepted_daugters_for_track = 0;
  m_num_total_daugters_for_track = 0;

  m_info->update_min_time(iTrk->GetGlobalTime());
  m_step_n = 0;

  g4_pos_mom_to_cms_birth(iTrk, p);
  p.m_birth_proc = process_index(iTrk->GetCreatorProcess());
  // Also set last ... in case there is no step info. In End momentum = 0!
  g4_pos_mom_to_cms_last(iTrk, p);
  p.m_mass = Mom_fac * pd->GetPDGMass();
  p.m_pdg = pd->GetPDGEncoding();
  p.m_charge = (int) pd->GetPDGCharge();

  if (m_output_path) {
    const G4ThreeVector &pos = iTrk->GetPosition();
    const G4ThreeVector mom = iTrk->GetMomentum();
    const G4ThreeVector &dir = iTrk->GetMomentumDirection();
    const double x[3] = {pos.x(), pos.y(), pos.z()}, q[3] = {mom.x(), mom.y(), mom.z()};
    add_path_point(m_id, x, iTrk->GetGlobalTime(), q, iTrk->GetTotalEnergy());
    m_path_dir[0] = dir.x(); m_path_dir[1] = dir.y(); m_path_dir[2] = dir.z();
    m_path_len_since = 0;
  }
}

void G4Snitch::update(const G4Step *iStep)
{
  G4Track *iTrk = iStep->GetTrack();

  G4StepPoint *sp1 = iStep->GetPreStepPoint();
  G4VPhysicalVolume *pv1 = sp1->GetPhysicalVolume();
  G4LogicalVolume *lv1 = pv1->GetLogicalVolume();
  // If the following is non-null, the step is in sensitive material.
  G4VSensitiveDetector *sd1 = lv1->GetSensitiveDetector();

  // The post-step point has no volume when the step leaves the world.
  G4StepPoint *sp2 = iStep->GetPostStepPoint();
  G4VPhysicalVolume *pv2 = sp2->GetPhysicalVolume();
  G4LogicalVolume *lv2 = pv2 ? pv2->GetLogicalVolume() : nullptr;
  G4VSensitiveDetector *sd2 = lv2 ? lv2->GetSensitiveDetector() : nullptr;

  // The material of this step (its pre-step point's): stored steps get the sums at their start, path points after it.
  double step_x0 = 0, step_lambda = 0;
  {
    const G4Material *mat = sp1->GetMaterial();
    const double len = iStep->GetStepLength();
    step_x0 = len / mat->GetRadlen();
    step_lambda = len / mat->GetNuclearInterLength();
  }

  // We process energy deposits in any case.
  {
    double tot  = iStep->GetTotalEnergyDeposit();
    double niel = iStep->GetNonIonizingEnergyDeposit();
    bool is_sensitive = (sd1 != nullptr);

    auto &psteps = particle_steps(m_id);
    psteps.add_edep(tot, niel, is_sensitive);

    // Store all, eventually merge inert ones up to some
    // thresholds on distance / summed up energy / volume-id.

    if (!is_sensitive && !m_tracking && m_inert_aggregate_descendants) {
      // An untracked descendant's inert step: summed per volume under the tracked ancestor.
      if (tot > 0 && tot >= m_inert_step_ecut)
        add_inert_step(m_id, iStep, tot, niel);
    } else if (( is_sensitive && m_output_sensitive_steps) ||
               (!is_sensitive && m_output_inert_steps)) {
      G4S_Step step;
      g4_pos_to_cms_begin(sp1, step);
      g4_pos_to_cms_end(sp2, step);
      {
        const G4ThreeVector mom = sp2->GetMomentum();
        step.m_p_end.Set(Mom_fac * mom.x(), Mom_fac * mom.y(), Mom_fac * mom.z());
        step.m_p_end.fT = Mom_fac * sp2->GetTotalEnergy();
      }
      step.m_edep.set_edep(tot, niel);
      step.m_g4_pid = iTrk->GetTrackID();
      step.m_pdg = iTrk->GetDefinition()->GetPDGEncoding();
      step.m_proc = process_index(sp2->GetProcessDefinedStep());
      step.m_x0_beg = m_mat_x0;
      step.m_lambda_beg = m_mat_lambda;
      if (is_sensitive) {
        // A step with no deposit is stored only on request; one with a deposit when it reaches the cut.
        if (m_output_sensitive_steps && (tot > 0 ? tot >= m_sensitive_step_ecut : m_output_sensitive_zero_steps)) {
          step.m_volume = volume_index(sp1);
          // Not under no_det_id_volumes: CMS makes HF's hits from a shower library and fibre light, not from steps,
          // so no step there has a hit's id. (HCalSD's id for an HF step would also draw a random number, below.)
          // The sensitive detector has processed this step already. For CaloSD the members its setDetUnitId() sets
          // (HGCalSD's fraction_ and calibCell_) are set again by the next ProcessHits() before they are read.
          // HCalSD's setDetUnitId() draws a random number (HF: depth 1 or 3, 2 or 4) and changed the simulation from
          // there on; the engine's state is kept around the call, so the snitch takes nothing from it.
          if (m_sensitive_det_id && !m_vol_no_det_id[step.m_volume]) {
            if (auto *csd = dynamic_cast<SensitiveDetector *>(sd1)) {
              CLHEP::HepRandomEngine *engine = G4Random::getTheEngine();
              const std::vector<unsigned long> engine_state = engine->put();
              step.m_det_id = csd->setDetUnitId(iStep);
              engine->get(engine_state);
            }
          }
          psteps.m_steps_sensitive.emplace_back(step);
        }
      } else {
        if (m_output_inert_steps && tot >= m_inert_step_ecut) {
          step.m_volume = volume_index(sp1);
          psteps.m_steps_inert.emplace_back(step);
        }
      }
    }
  }

  m_mat_x0 += step_x0;
  m_mat_lambda += step_lambda;

  // The secondaries created in this step. GetNumberOfSecondariesInCurrentStep() can count a secondary that is not
  // appended to the track's secondary vector (seen with eBrem in Geant4 11.4, CMSSW_20_1_0_pre5), so indexing the
  // cumulative vector GetSecondary() with a running count goes wrong; this vector holds exactly the step's own.
  const std::vector<const G4Track *> *step_secondaries = iStep->GetSecondaryInCurrentStep();
  const int n_step_secondaries = step_secondaries ? (int)step_secondaries->size() : 0;

  if (m_tracking == false) {
    ++m_steps_skipped;
    int nd_skipped = n_step_secondaries;
    if (nd_skipped > 0) {
      m_daughters_skipped += nd_skipped;
      particle(m_id).m_n_sub_daughters_skipped += nd_skipped;
    }
    return;
  }

  ++m_step_n;
  // const G4VProcess *p1 = sp1->GetProcessDefinedStep();
  // std::string pn1(p1 ? p1->GetProcessName() : "null");
  const G4VProcess *p2 = sp2->GetProcessDefinedStep();
  G4ProcessType pt2(p2 ? p2->GetProcessType() : fNotDefined);
  std::string pn2(p2 ? p2->GetProcessName() : "null");

  if (iTrk->GetKineticEnergy() > m_last_state_ekin_cut) {
    g4_pos_mom_to_cms_last(iTrk, particle(m_id));
  }

  // Thinned path: the step's end point when the direction turned by more than path_max_angle or the path grew by
  // more than path_max_length since the last stored point. Done before the secondaries below grow the vectors.
  if (m_output_path) {
    const G4ThreeVector &dir = sp2->GetMomentumDirection();
    m_path_len_since += iStep->GetStepLength();
    const double cosang = dir.x() * m_path_dir[0] + dir.y() * m_path_dir[1] + dir.z() * m_path_dir[2];
    if (cosang < m_path_cos_max_angle || Pos_fac * m_path_len_since > m_path_max_length) {
      const G4ThreeVector &pos = sp2->GetPosition();
      const G4ThreeVector mom = sp2->GetMomentum();
      const double x[3] = {pos.x(), pos.y(), pos.z()}, q[3] = {mom.x(), mom.y(), mom.z()};
      add_path_point(m_id, x, sp2->GetGlobalTime(), q, sp2->GetTotalEnergy());
      m_path_dir[0] = dir.x(); m_path_dir[1] = dir.y(); m_path_dir[2] = dir.z();
      m_path_len_since = 0;
    }
  }

  if (m_verbose && (m_verbose_transport || pt2 != fTransportation)) {
    printf("  %3d. t=%.3f E_kin=%.1f E_tot=%.1f v=%.4f  --  vol-%s '%s' sensitive=%d/%d Edep=%f Eniel=%f N_sec_curr=%d, N_sec=%d\n",
        m_step_n,
        iTrk->GetGlobalTime(), iTrk->GetKineticEnergy(), iTrk->GetTotalEnergy(), iTrk->GetVelocity(),
        (pv1 == pv2) ? "same" : "chng", pn2.c_str(),
        (sd1) ? 1 : 0, (sd2) ? 1 : 0,
        iStep->GetTotalEnergyDeposit(), iStep->GetNonIonizingEnergyDeposit(),
        n_step_secondaries,
        (int) iStep->GetSecondary()->size()
    );
  }

  // Process any new secondaries, prepare map entries for those that pass. Resize is done in EndOfTrack.
  int nss = n_step_secondaries;
  if (nss > 0) {
    for (int i = 0; i < nss; ++i) {
      const G4Track *gt = (*step_secondaries)[i];
      auto *gp = gt->GetParticleDefinition();
      int new_gid = m_last_g4_id + 1 + m_num_total_daugters_for_track;
      if (filter(gt)) {
        ++particle(m_id).m_n_daughters_skipped;
        if (m_verbose_skip) {
          printf("------> NOT Scheduling insertion of secondary gid %d, stable=%s pdg=%d E_k=%.1f\n",
                 new_gid, b2yn(gp->GetPDGStable()), gp->GetPDGEncoding(), gt->GetKineticEnergy());
        }
      }
      else
      {
        if (m_verbose)
          printf("------> Scheduling insertion of secondary gid %d to %d, stable=%s pdg=%d E_k = %.1f\n",
                new_gid, m_vec_size,
                b2yn(gp->GetPDGStable()), gp->GetPDGEncoding(), gt->GetKineticEnergy());
        m_gtp2vid.insert(std::make_pair(gt, m_vec_size));

        if (m_vec_size >= m_vec_capacity) {
          if (2 * m_id_current_primary < m_g4_num_primaries || m_id_current_primary <= 1) {
            m_vec_capacity *= 2;
          } else {
            // The extrapolation from the finished primaries goes negative when the current primary alone has more
            // accepted tracks than the capacity (a TeV hadron into HF among ten primaries), and vector::reserve
            // threw; double instead whenever it does not grow the vectors.
            const int est = (m_vec_capacity - m_primary_n_accepted_tracks) / (m_id_current_primary - 1) * m_g4_num_primaries;
            m_vec_capacity = (est > m_vec_size) ? est : 2 * m_vec_capacity;
          }
          if (m_verbose)
            printf("MMMMMMM Growing vec memory from %d to %d\n", m_vec_size, m_vec_capacity);
          m_part_vec->reserve(m_vec_capacity);
          steps_reserve(m_vec_capacity);
        }
        G4S_Particle p;
        g4_pos_mom_to_cms_birth(gt, p);
        p.m_birth_proc = process_index(gt->GetCreatorProcess());
        p.m_mass = Mom_fac * gp->GetPDGMass();
        p.m_pdg = gp->GetPDGEncoding();
        p.m_charge = (int) gp->GetPDGCharge();
        p.m_parent = m_id;
        p.m_g4_level = stack_level() + 1;
        m_part_vec->emplace_back(p);
        steps_emplace_back();
        ++m_vec_size;
        ++m_num_accepted_daugters_for_track;
        ++m_primary_n_accepted_tracks;
        ++m_event_tracks_accepted;
      }
      ++m_num_total_daugters_for_track;
    }
  }
}

void G4Snitch::update(const EndOfTrack* eot)
{
  if (m_tracking == false) {
    return;
  }

  auto &p = particle(m_id);
  if (m_num_accepted_daugters_for_track > 0) {
    p.m_daughters_end = m_vec_size;
  }

  const G4Track *iTrk = (*eot)();
  auto &pos = iTrk->GetPosition();
  auto  mom = iTrk->GetMomentum();

  // Death: the track's final state, and the process of its last step.
  g4_pos_mom_to_cms_death(iTrk, p);
  {
    const G4Step *last_step = iTrk->GetStep();
    const G4StepPoint *post = last_step ? last_step->GetPostStepPoint() : nullptr;
    p.m_death_proc = process_index(post ? post->GetProcessDefinedStep() : nullptr);
    p.m_death_volume = last_step ? volume_index(last_step->GetPreStepPoint()) : -1;
  }

  // The path's last point: the track's end, unless the last step already stored it.
  if (m_output_path && m_path_len_since > 0) {
    const double x[3] = {pos.x(), pos.y(), pos.z()}, q[3] = {mom.x(), mom.y(), mom.z()};
    add_path_point(m_id, x, iTrk->GetGlobalTime(), q, iTrk->GetTotalEnergy());
    m_path_len_since = 0;
  }
  if (m_verbose)
    printf("  END         pos=%.1f,%.1f,%.1f  mom=%.1f,%.1f,%.1f  t=%.3f E_kin=%.1f E_tot=%.1f v=%.4f; daughters: %d -> %d\n",
      pos.x(), pos.y(), pos.z(),
      mom.x(), mom.y(), mom.z(),
      iTrk->GetGlobalTime(), iTrk->GetKineticEnergy(), iTrk->GetTotalEnergy(), iTrk->GetVelocity(),
      p.m_daughters_begin, p.m_daughters_end
    );

  m_info->update_max_time(iTrk->GetGlobalTime());
}

void G4Snitch::update(const EndOfEvent *eoe)
{
  std::cout << "++ signal EndOfEvent -- fill the tree\n";
  printf(" tracks_accepted=%d, m_event_tracks_skipped=%d\n",
         m_event_tracks_accepted, m_event_tracks_skipped);
  m_total_tracks_accepted += m_event_tracks_accepted;
  m_total_tracks_skipped += m_event_tracks_skipped;

  flush_inert_sums();
  select_output();
  m_tree->Fill();
  reset_output_structs();
}

void G4Snitch::select_output() {
  if (m_select_track_id > 0) {
    // The subtree of the selected particle, in the same layout: slot 0 the administrative mother with the selected
    // particle as its only daughter, then each particle's daughters as one contiguous block.
    std::vector<G4S_Particle> &pv = *m_part_vec;
    int root = -1;
    for (int i = 1; i < (int)pv.size(); ++i)
      if (pv[i].m_was_tracked && pv[i].m_g4_id == m_select_track_id) {
        root = i;
        break;
      }
    std::vector<G4S_Particle> np;
    std::vector<G4S_ParticleSteps> ns;
    std::vector<G4S_PathBField> nb;
    std::vector<int> order;  // old index of each new slot
    np.push_back(pv[0]);
    ns.emplace_back();
    if (m_output_path_bfield)
      nb.emplace_back();
    order.push_back(0);
    np[0].m_daughters_begin = 1;
    np[0].m_daughters_end = 1;
    if (root >= 0) {
      np[0].m_daughters_end = 2;
      order.push_back(root);
      np.push_back(pv[root]);
      np[1].m_parent = 0;
      for (size_t k = 1; k < order.size(); ++k) {
        const int o = order[k];
        ns.push_back(std::move((*m_psteps_vec)[o]));
        if (m_output_path_bfield)
          nb.push_back(std::move((*m_pbfield_vec)[o]));
        const int db = pv[o].m_daughters_begin, de = pv[o].m_daughters_end;
        if (db < 0 || de <= db) {
          if (db >= 0)
            np[k].m_daughters_begin = np[k].m_daughters_end = (int)np.size();
          continue;
        }
        np[k].m_daughters_begin = (int)np.size();
        for (int j = db; j < de; ++j) {
          order.push_back(j);
          np.push_back(pv[j]);
          np.back().m_parent = (int)k;
        }
        np[k].m_daughters_end = (int)np.size();
      }
    } else {
      edm::LogWarning("G4Snitch") << "select_track_id " << m_select_track_id << ": no tracked particle with that id";
    }
    m_part_vec->swap(np);
    m_psteps_vec->swap(ns);
    if (m_output_path_bfield)
      m_pbfield_vec->swap(nb);
  }
  if (m_select_cone.size() == 3) {
    const double eta0 = m_select_cone[0], phi0 = m_select_cone[1], dr2 = m_select_cone[2] * m_select_cone[2];
    auto inside = [&](double x, double y, double z) {
      const double r = std::hypot(x, y);
      if (r <= 0)
        return false;
      const double deta = std::asinh(z / r) - eta0;
      const double dphi = std::remainder(std::atan2(y, x) - phi0, 2 * M_PI);
      return deta * deta + dphi * dphi <= dr2;
    };
    for (G4S_ParticleSteps &ps : *m_psteps_vec) {
      for (std::vector<G4S_Step> *v : {&ps.m_steps_sensitive, &ps.m_steps_inert})
        v->erase(std::remove_if(v->begin(), v->end(),
                                [&](const G4S_Step &st) {
                                  return !inside(0.5 * (st.m_x_beg.fX + st.m_x_end.fX),
                                                 0.5 * (st.m_x_beg.fY + st.m_x_end.fY),
                                                 0.5 * (st.m_x_beg.fZ + st.m_x_end.fZ));
                                }),
                 v->end());
      ps.m_inert_sums.erase(std::remove_if(ps.m_inert_sums.begin(), ps.m_inert_sums.end(),
                                           [&](const G4S_InertSum &is) {
                                             return !inside(is.m_x_mean.fX, is.m_x_mean.fY, is.m_x_mean.fZ);
                                           }),
                            ps.m_inert_sums.end());
    }
  }
}

void G4Snitch::update(const EndOfRun *eor)
{
  std::cout << "++ signal EndOfRun  -- closing the root file\n";
  printf(" m_total_tracks_accepted=%d, m_total_tracks_skipped=%d\n",
         m_total_tracks_accepted, m_total_tracks_skipped);
  write_tree_close_file();
}

// this never gets called
void G4Snitch::update(const EndOfJob* eoj)
{
  std::cout <<"++ signal EndOfJob";
}
