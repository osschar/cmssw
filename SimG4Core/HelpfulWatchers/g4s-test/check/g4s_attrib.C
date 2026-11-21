// g4s_attrib.C(file): are the steps filed under a tracked particle P from P itself or its descendants?
// Geant4 numbers a track's daughters when the track ends, so every descendant of P has a larger track id than P.
// A step under P whose m_g4_pid is smaller than P's m_g4_id cannot come from a descendant: misattributed.
#include <cstdio>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "SimG4Core/HelpfulWatchers/interface/G4SnitchDataFormat.h"
void g4s_attrib(const char *fn) {
  TFile f(fn);
  TTree *t = (TTree *)f.Get("T");
  std::vector<G4S_Particle> *p = nullptr;
  std::vector<G4S_ParticleSteps> *s = nullptr;
  t->SetBranchAddress("p", &p);
  t->SetBranchAddress("s", &s);
  const char *kind[2] = {"sensitive", "inert"};
  long n[2][3] = {}; double e[2][3] = {};  // [kind][own, larger id, smaller id]
  long n_part_bad = 0, n_part = 0;
  for (long ev = 0; ev < t->GetEntries(); ++ev) {
    t->GetEntry(ev);
    for (size_t i = 1; i < p->size(); ++i) {
      const G4S_Particle &pp = (*p)[i];
      if (!pp.m_was_tracked) continue;
      ++n_part;
      bool bad = false;
      for (int k = 0; k < 2; ++k) {
        const std::vector<G4S_Step> &v = k == 0 ? (*s)[i].m_steps_sensitive : (*s)[i].m_steps_inert;
        for (const G4S_Step &st : v) {
          int c = st.m_g4_pid == pp.m_g4_id ? 0 : (st.m_g4_pid > pp.m_g4_id ? 1 : 2);
          ++n[k][c]; e[k][c] += st.m_edep.m_total;
          if (c == 2) bad = true;
        }
      }
      if (bad) ++n_part_bad;
    }
  }
  printf("%ld events, %ld tracked particles, %ld with steps of a smaller track id\n", (long)t->GetEntries(), n_part, n_part_bad);
  printf("%-10s %12s %12s %12s   %10s %10s %10s (MeV)\n", "steps", "own", "larger id", "smaller id", "own", "larger", "smaller");
  for (int k = 0; k < 2; ++k)
    printf("%-10s %12ld %12ld %12ld   %10.1f %10.1f %10.1f\n", kind[k], n[k][0], n[k][1], n[k][2], e[k][0], e[k][1], e[k][2]);
}
