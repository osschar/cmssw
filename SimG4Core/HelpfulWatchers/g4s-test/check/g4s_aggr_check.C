// g4s_aggr_check.C(fileA, fileB): A ran with inert_aggregate_descendants off, B with it on, same events.
// Per tracked particle: B's inert sums must equal A's descendant inert steps (m_g4_pid != the particle's m_g4_id) in
// energy and step count; B's own inert steps and both files' sensitive steps must equal A's.
#include <cmath>
#include <cstdio>
#include <vector>
#include "SimG4Core/HelpfulWatchers/interface/G4SnitchReader.h"
void g4s_aggr_check(const char *fa, const char *fb) {
  G4SnitchFile a(fa), b(fb);  // TTree or RNTuple files
  std::vector<G4S_Particle> *&pa = a.p, *&pb = b.p;
  std::vector<G4S_ParticleSteps> *&sa = a.s, *&sb = b.s;
  long n_part = 0, n_bad_e = 0, n_bad_n = 0, n_bad_own = 0, n_bad_sens = 0, n_sums = 0, n_desc_steps = 0;
  double e_desc = 0, e_sums = 0, max_rel = 0;
  for (long ev = 0; ev < a.n_events(); ++ev) {
    a.load(ev); b.load(ev);
    if (a.ev != b.ev) {
      printf("entry %ld: event %d in A, %d in B; the files must hold the same events in the same order\n", ev, a.ev, b.ev);
      return;
    }
    if (pa->size() != pb->size()) { printf("event %ld: particle counts differ %zu %zu\n", ev, pa->size(), pb->size()); return; }
    for (size_t i = 1; i < pa->size(); ++i) {
      if (!(*pa)[i].m_was_tracked) continue;
      ++n_part;
      const int gid = (*pa)[i].m_g4_id;
      double ea = 0; long na = 0, own_a = 0;
      for (const G4S_Step &st : (*sa)[i].m_steps_inert) {
        if (st.m_g4_pid == gid) ++own_a;
        else { ea += st.m_edep.m_total; ++na; }
      }
      double eb = 0; long nb = 0;
      for (const G4S_InertSum &is : (*sb)[i].m_inert_sums) { eb += is.m_edep.m_total; nb += is.m_n_steps; ++n_sums; }
      long own_b = (*sb)[i].m_steps_inert.size();
      e_desc += ea; e_sums += eb; n_desc_steps += na;
      double rel = ea > 0 ? std::abs(eb - ea) / ea : std::abs(eb);
      if (rel > max_rel) max_rel = rel;
      if (rel > 1e-9) ++n_bad_e;
      if (na != nb) ++n_bad_n;
      if (own_a != own_b) ++n_bad_own;
      if ((*sa)[i].m_steps_sensitive.size() != (*sb)[i].m_steps_sensitive.size()) ++n_bad_sens;
    }
  }
  printf("%lld events, %ld tracked particles\n", (long long)a.n_events(), n_part);
  printf("descendant inert steps in A: %ld, %.1f MeV; sums in B: %ld entries, %.1f MeV\n", n_desc_steps, e_desc, n_sums, e_sums);
  printf("particles differing: energy %ld (max rel %.2e), step count %ld, own inert steps %ld, sensitive steps %ld\n",
         n_bad_e, max_rel, n_bad_n, n_bad_own, n_bad_sens);
  printf("volume table of B: %lld volumes\n", (long long)b.column<int>("V", "idx").size());
}
