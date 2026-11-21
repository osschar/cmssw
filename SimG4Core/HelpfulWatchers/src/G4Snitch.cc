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
#include <iostream>
#include <typeinfo>

// user include files
#include "SimG4Core/HelpfulWatchers/interface/G4Snitch.h"
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/EndOfTrack.h"

#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VTouchable.hh"
#include "G4Material.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "SimG4Core/SensitiveDetector/interface/SensitiveDetector.h"
#include "Randomize.hh"
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

  // Cuts in G4 units
  // Minimum kinetic energy at particle creation for that particle to be tracked.
  // Unstable particles are treated correctly, see filter().
  constexpr double E_cut = 10.0;
  // Kinetic energy below which the end position/momentum of a particle
  // are no longer updated in the G4S_Particle structure.
  constexpr double E_final_cut = 1.0;

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
        m_inert_step_ecut(pSet.getUntrackedParameter<double>("inert_step_ecut", 1e-6)),
        m_inert_aggregate_descendants(pSet.getUntrackedParameter<bool>("inert_aggregate_descendants", true)),
        m_sensitive_det_id(pSet.getUntrackedParameter<bool>("sensitive_det_id", true)),

        m_output_path(pSet.getUntrackedParameter<bool>("output_path", true)),
        m_path_max_angle(pSet.getUntrackedParameter<double>("path_max_angle", 0.05)),
        m_path_max_length(pSet.getUntrackedParameter<double>("path_max_length", 5.0)),
        m_output_path_bfield(pSet.getUntrackedParameter<bool>("output_path_bfield", false))
{
  m_path_cos_max_angle = std::cos(m_path_max_angle);
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
  for (int d = depth; d >= 0; --d) {
    if (!path.empty())
      path += '/';
    path += th->GetVolume(d)->GetName();
    path += ':';
    path += std::to_string(th->GetCopyNumber(d));
  }
  m_vol_paths.push_back(path);
  const G4LogicalVolume *lv = th->GetVolume(0)->GetLogicalVolume();
  m_vol_sd.push_back(sd_index(lv->GetSensitiveDetector()));
  const G4Material *mat = lv->GetMaterial();
  m_vol_materials.push_back(mat->GetName());
  m_vol_x0.push_back(Pos_fac * mat->GetRadlen());
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
  const double x[4] = {0.5 * Pos_fac * (p1.x() + p2.x()), 0.5 * Pos_fac * (p1.y() + p2.y()),
                       0.5 * Pos_fac * (p1.z() + p2.z()), 0.5 * (sp1->GetGlobalTime() + sp2->GetGlobalTime())};
  int ai;
  auto it = m_inert_acc_index.find(key);
  if (it == m_inert_acc_index.end()) {
    ai = (int)m_inert_acc.size();
    m_inert_acc_index.emplace(key, ai);
    InertAcc a{slot, vol, 0, 0, 0, {x[0], x[1], x[2], x[3]}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    m_inert_acc.push_back(a);
  } else {
    ai = it->second;
  }
  InertAcc &a = m_inert_acc[ai];
  ++a.n;
  a.e += tot;
  a.niel += niel;
  for (int k = 0; k < 4; ++k) {
    const double d = x[k] - a.ref[k];
    a.sw[k] += tot * d;
    a.sw2[k] += tot * d * d;
  }
}

void G4Snitch::flush_inert_sums() {
  // Every sum has e > 0: only steps with a deposit above zero are added.
  for (const InertAcc &a : m_inert_acc) {
    G4S_InertSum is;
    is.m_volume = a.vol;
    is.m_n_steps = a.n;
    is.m_edep.set_edep(a.e, a.niel);
    double m[4], r[4];
    for (int k = 0; k < 4; ++k) {
      const double mu = a.sw[k] / a.e;
      m[k] = a.ref[k] + mu;
      r[k] = std::sqrt(std::max(0.0, a.sw2[k] / a.e - mu * mu));
    }
    is.m_x_mean.Set(m[0], m[1], m[2]);
    is.m_x_mean.fT = m[3];
    is.m_x_rms.Set(r[0], r[1], r[2]);
    is.m_x_rms.fT = r[3];
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
  m_file = TFile::Open("G4Snitch.root", "RECREATE");
  m_tree = new TTree("T", "G4 Dump Tree of some kind");

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
    int idx = 0, sd = -1;
    std::string path, material;
    double x0 = 0, density = 0, x = 0, y = 0, z = 0;
    vtree.Branch("idx", &idx);
    vtree.Branch("path", &path);
    vtree.Branch("material", &material);
    vtree.Branch("x0", &x0);
    vtree.Branch("density", &density);
    vtree.Branch("x", &x);
    vtree.Branch("y", &y);
    vtree.Branch("z", &z);
    vtree.Branch("sd", &sd);
    for (idx = 0; idx < (int)m_vol_paths.size(); ++idx) {
      path = m_vol_paths[idx];
      material = m_vol_materials[idx];
      x0 = m_vol_x0[idx];
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
  m_file->Close();
  delete m_file;
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
  if ( ! pd->GetPDGStable() && pd->GetPDGEncoding() != 2112)
    return false;
  if (pd->GetPDGEncoding() == -11 && t->GetTotalEnergy() >= E_cut)
    return false;
  return t->GetKineticEnergy() < E_cut;
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
      if (is_sensitive) {
        // A step with no deposit is stored only on request; one with a deposit when it reaches the cut.
        if (m_output_sensitive_steps && (tot > 0 ? tot >= m_sensitive_step_ecut : m_output_sensitive_zero_steps)) {
          step.m_volume = volume_index(sp1);
          // The sensitive detector has processed this step already. For CaloSD the members its setDetUnitId() sets
          // (HGCalSD's fraction_ and calibCell_) are set again by the next ProcessHits() before they are read.
          // HCalSD's setDetUnitId() draws a random number (HF: depth 1 or 3, 2 or 4) and changed the simulation from
          // there on; the engine's state is kept around the call, so the snitch takes nothing from it.
          if (m_sensitive_det_id) {
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

  if (iTrk->GetKineticEnergy() > E_final_cut) {
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
          if (2 * m_id_current_primary < m_g4_num_primaries || m_id_current_primary <= 1)
            m_vec_capacity *= 2;
          else
            m_vec_capacity = (m_vec_capacity - m_primary_n_accepted_tracks) / (m_id_current_primary - 1) * m_g4_num_primaries;
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
  m_tree->Fill();
  reset_output_structs();
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
