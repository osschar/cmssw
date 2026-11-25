#ifndef HelpfulWatchers_G4SnitchReader_h
#define HelpfulWatchers_G4SnitchReader_h
// G4SnitchReader.h: a G4Snitch output file read the same way whatever its format: RNTuples (since 2026-11-24) or the
// TTrees of older files. Header only, for macros (ACLiC, after gSystem->Load("libRootG4Snitch.so")) and programs.
//
//   G4SnitchFile sf("G4Snitch.root");
//   for (long k = 0; k < sf.n_events(); ++k) {
//     sf.load(k);
//     sf.ev, *sf.i, *sf.p, *sf.s, sf.b (nullptr when the file has no field)
//   }
//   auto names = sf.column<std::string>("P", "name");  // a column of a table P, V, D, C, by row index
//
// T's entries are in the order the threads finished their events: match events by sf.ev.

#include "SimG4Core/HelpfulWatchers/interface/G4SnitchDataFormat.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

class G4SnitchFile {
public:
  int ev = -1;
  G4S_Info *i = nullptr;
  std::vector<G4S_Particle> *p = nullptr;
  std::vector<G4S_ParticleSteps> *s = nullptr;
  std::vector<G4S_PathBField> *b = nullptr;

  explicit G4SnitchFile(const std::string &fn) : m_fn(fn), m_file(TFile::Open(fn.c_str())) {
    if (!m_file || m_file->IsZombie())
      throw std::runtime_error("G4SnitchFile: cannot open " + fn);
    m_tree = m_file->Get<TTree>("T");
    if (m_tree) {
      m_tree->SetBranchAddress("i", &i);
      m_tree->SetBranchAddress("p", &p);
      m_tree->SetBranchAddress("s", &s);
      if (m_tree->GetBranch("ev"))
        m_tree->SetBranchAddress("ev", &ev);
      if (m_tree->GetBranch("b"))
        m_tree->SetBranchAddress("b", &b);
      m_n = m_tree->GetEntries();
    } else {
      m_reader = ROOT::RNTupleReader::Open("T", fn);
      m_entry = m_reader->CreateEntry();
      m_ev = m_entry->GetPtr<int>("ev").get();
      i = m_entry->GetPtr<G4S_Info>("i").get();
      p = m_entry->GetPtr<std::vector<G4S_Particle>>("p").get();
      s = m_entry->GetPtr<std::vector<G4S_ParticleSteps>>("s").get();
      if (m_reader->GetDescriptor().FindFieldId("b") != ROOT::kInvalidDescriptorId)
        b = m_entry->GetPtr<std::vector<G4S_PathBField>>("b").get();
      m_n = m_reader->GetNEntries();
    }
  }

  bool is_ttree() const { return m_tree != nullptr; }
  long n_events() const { return m_n; }
  TFile &file() { return *m_file; }

  // Load event k into ev, i, p, s, b.
  void load(long k) {
    if (m_tree) {
      m_tree->GetEntry(k);
    } else {
      m_reader->LoadEntry(k, *m_entry);
      ev = *m_ev;
    }
  }

  // Column col of table P, V, D or C, indexed by row (the tables' idx is the row number).
  template <typename T>
  std::vector<T> column(const char *table, const char *col) {
    std::vector<T> out;
    if (auto *t = m_file->Get<TTree>(table)) {
      TTreeReader r(t);
      TTreeReaderValue<T> v(r, col);
      while (r.Next())
        out.push_back(*v);
    } else {
      auto r = ROOT::RNTupleReader::Open(table, m_fn);
      auto view = r->GetView<T>(col);
      for (auto k : r->GetEntryRange())
        out.push_back(view(k));
    }
    return out;
  }

private:
  std::string m_fn;
  std::unique_ptr<TFile> m_file;
  TTree *m_tree = nullptr;
  std::unique_ptr<ROOT::RNTupleReader> m_reader;
  std::unique_ptr<ROOT::REntry> m_entry;
  int *m_ev = nullptr;
  long m_n = 0;
};

#endif
