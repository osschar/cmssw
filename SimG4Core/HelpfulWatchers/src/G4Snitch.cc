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
#include <iostream>

// user include files
#include "SimG4Core/HelpfulWatchers/interface/G4Snitch.h"
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/EndOfTrack.h"

#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"

#include <TFile.h>
#include <TTree.h>

namespace {

  // Cuts in G4 units
  constexpr double E_cut = 10.0;
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
  void g4_pos_mom_to_cms_begin(const G4Track *t, G4S_Particle &part) {
    g4_pos_mom_to_cms(t, part.m_x_beg, part.m_p_beg);
  }
  void g4_pos_mom_to_cms_end(const G4Track *t, G4S_Particle &part) {
    g4_pos_mom_to_cms(t, part.m_x_end, part.m_p_end);
  }

  const char* b2yn(bool b) { return b ? "yes" : "no"; }
}

//------------------------------------------------------------------------------

G4Snitch::G4Snitch(const edm::ParameterSet &pSet)
      : m_verbose(pSet.getUntrackedParameter<bool>("verbose", false)),
        m_verbose_stack_level(pSet.getUntrackedParameter<bool>("verbose_stack_level", false)),
        m_verbose_transport(pSet.getUntrackedParameter<bool>("verbose_transport", false)),
        m_verbose_skip(pSet.getUntrackedParameter<bool>("verbose_skip", false)),
        m_verbose_skip_with_ids(pSet.getUntrackedParameter<bool>("verbose_skip_with_ids", false)) {
}

// this never gets called, sigh
G4Snitch::~G4Snitch() {
  std::cout <<"G4Snitch destructed\n";
}

void G4Snitch::open_file_tree() {
  m_file = TFile::Open("G4Snitch.root", "RECREATE");
  m_tree = new TTree("T", "G4 Dump Tree of some kind");

  m_part_vec.reset(new std::vector<G4S_Particle>);
  m_part_vec->reserve(16384);
  m_vec_capacity = 16384;
  m_tree->Branch("p", m_part_vec.get());
  m_info.reset(new G4S_Info);
  m_tree->Branch("i", m_info.get());
}

void G4Snitch::write_tree_close_file() {
  m_tree->Write();
  delete m_tree;
  m_file->Close();
  delete m_file;
}

void G4Snitch::reset_output_structs() {
  m_gid2vid.clear();
  m_part_vec->clear();
  m_vec_size = 0;
  m_gtp2vid.clear();
  m_stack.clear();
  m_info->reset();
}

//------------------------------------------------------------------------------

bool G4Snitch::filter(const G4Track *t) const {
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
    // They will be reversed into id order, ie, primary i --> slot i - 1;
    m_g4_num_primaries = gid;
    m_id_current_primary = 0;

    m_part_vec->resize(m_g4_num_primaries + 1);
    m_vec_size = m_g4_num_primaries + 1;
    m_info->m_n_primaries = m_g4_num_primaries; // XXX might not be true !!! fix at end !!!
    m_event_tracks_accepted = m_event_tracks_skipped = 0;

    m_primary_n_accepted_tracks = 0;

    m_last_g4_id = -1;

    /// insert 0, setup mamma particle
    m_stack.push_back(0);
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
    if (m_verbose_stack_level && (m_tracking || m_verbose_skip_with_ids))
      printf("--- UP 1 level to %d\n", stack_level());
  }
  m_stack.push_back(gid);
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

    m_gid2vid.clear(); // no need to schlep finished primaries' ids along.
    m_gid2vid.insert(std::make_pair(gid, m_id_current_primary));
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
  auto &p = particle(m_id);
  p.m_daughters_begin = p.m_daughters_end = m_vec_size;
  p.m_g4_id = gid;
  p.m_was_tracked = true;

  m_num_accepted_daugters_for_track = 0;
  m_num_total_daugters_for_track = 0;

  m_info->update_min_time(iTrk->GetGlobalTime());
  m_step_n = 0;

  g4_pos_mom_to_cms_begin(iTrk, p);
  // Also set end ... in case there is no step info. In End momentum = 0!
  g4_pos_mom_to_cms_end(iTrk, p);
  p.m_mass = Mom_fac * pd->GetPDGMass();
  p.m_pdg = pd->GetPDGEncoding();
  p.m_charge = (int) pd->GetPDGCharge();
}

void G4Snitch::update(const G4Step *iStep)
{
  if (m_tracking == false) {
    ++m_steps_skipped;
    int nd_skipped = iStep->GetNumberOfSecondariesInCurrentStep();
    if (nd_skipped > 0) {
      m_daughters_skipped += nd_skipped;
      auto &p = particle(m_id);
      if (stack_top() == m_id)
        p.m_n_daughters_skipped += nd_skipped;
      else
        p.m_n_sub_daughters_skipped += nd_skipped;
    }
    return;
  }

  ++m_step_n;
  G4Track *iTrk = iStep->GetTrack();
  G4StepPoint *sp1 = iStep->GetPreStepPoint();
  G4VPhysicalVolume *pv1 = sp1->GetPhysicalVolume();
  // const G4VProcess *p1 = sp1->GetProcessDefinedStep();
  // std::string pn1(p1 ? p1->GetProcessName() : "null");

  G4StepPoint *sp2 = iStep->GetPostStepPoint();
  G4VPhysicalVolume *pv2 = sp2->GetPhysicalVolume();
  const G4VProcess *p2 = sp2->GetProcessDefinedStep();
  G4ProcessType pt2(p2 ? p2->GetProcessType() : fNotDefined);
  std::string pn2(p2 ? p2->GetProcessName() : "null");

  if (iTrk->GetKineticEnergy() > E_final_cut) {
    g4_pos_mom_to_cms_end(iTrk, particle(m_id));
  }

  if (m_verbose && (m_verbose_transport || pt2 != fTransportation))
    printf("  %3d. t=%.3f E_kin=%.1f E_tot=%.1f v=%.4f  --  vol-%s '%s' Edep=%f N_sec_curr=%d, N_sec=%d\n",
        m_step_n,
        iTrk->GetGlobalTime(), iTrk->GetKineticEnergy(), iTrk->GetTotalEnergy(), iTrk->GetVelocity(),
        (pv1 == pv2) ? "same" : "chng", pn2.c_str(),
        iStep->GetTotalEnergyDeposit(),
        (int) iStep->GetNumberOfSecondariesInCurrentStep(),
        (int) iStep->GetSecondary()->size()
    );

  // Process any new secondaries, prepare map entries for those that pass. Resize is done in EndOfTrack.
  int nss = (int) iStep->GetNumberOfSecondariesInCurrentStep();
  if (nss > 0) {
    auto &gsv = * iStep->GetSecondary();
    for (int i = 0; i < nss; ++i) {
      auto *gt = gsv[m_num_total_daugters_for_track];
      auto *gp = gt->GetParticleDefinition();
      int new_gid = m_last_g4_id + 1 + m_num_total_daugters_for_track;
      if (filter(gt)) {
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
        }
        G4S_Particle p;
        g4_pos_mom_to_cms_begin(gt, p);
        p.m_mass = Mom_fac * gp->GetPDGMass();
        p.m_pdg = gp->GetPDGEncoding();
        p.m_charge = (int) gp->GetPDGCharge();
        p.m_parent = m_id;
        p.m_g4_level = stack_level() + 1;
        m_part_vec->emplace_back(p);
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
