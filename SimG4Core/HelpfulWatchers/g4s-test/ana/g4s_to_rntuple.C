// g4s_to_rntuple.C(in, out): a feasibility test of RNTuple for the snitch's output. Converts a snitch TTree file
// (tree T: branches ev (optional), i, p, s, b (optional)) into an RNTuple "T" with the same classes, one entry and one
// cluster per event, the TTree's branches bound to the RNTuple model's objects (the pattern of mkFit-external's
// RntDumper: one set of objects for both formats). The tables P, V, D, C are copied as TTrees. Then reports the
// compression settings and sizes, reads every event back from both formats and compares per-event summaries (counts
// of particles, sensitive and inert steps, inert sums, path points, field points, and the summed deposits), and times
// reading one event (the middle one) from each, warm (file already opened once).
// Run in the CMSSW area's environment after gSystem->Load("libRootG4Snitch.so"), with ACLiC and
// -I$CMSSW_BASE/src.
#include "SimG4Core/HelpfulWatchers/interface/G4SnitchDataFormat.h"

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <ROOT/RNTupleWriter.hxx>

#include "TFile.h"
#include "TStopwatch.h"
#include "TTree.h"

#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace {
  struct Summary {
    long np = 0, nsens = 0, ninert = 0, nsums = 0, npath = 0, nb = 0;
    double esens = 0, einert = 0, esums = 0;
    bool operator==(const Summary &o) const {
      return np == o.np && nsens == o.nsens && ninert == o.ninert && nsums == o.nsums && npath == o.npath &&
             nb == o.nb && esens == o.esens && einert == o.einert && esums == o.esums;
    }
  };
  Summary summarise(const std::vector<G4S_Particle> &p, const std::vector<G4S_ParticleSteps> &s,
                    const std::vector<G4S_PathBField> *b) {
    Summary r;
    r.np = p.size();
    for (const auto &ps : s) {
      r.nsens += ps.m_steps_sensitive.size();
      r.ninert += ps.m_steps_inert.size();
      r.nsums += ps.m_inert_sums.size();
      r.npath += ps.m_path.size();
      for (const auto &x : ps.m_steps_sensitive)
        r.esens += x.m_edep.m_total;
      for (const auto &x : ps.m_steps_inert)
        r.einert += x.m_edep.m_total;
      for (const auto &x : ps.m_inert_sums)
        r.esums += x.m_edep.m_total;
    }
    if (b)
      for (const auto &x : *b)
        r.nb += x.m_b.size();
    return r;
  }
}  // namespace

void g4s_to_rntuple(const char *in, const char *out) {
  std::unique_ptr<TFile> fin(TFile::Open(in));
  TTree *t = (TTree *)fin->Get("T");
  const bool has_ev = t->GetBranch("ev"), has_b = t->GetBranch("b");

  auto model = ROOT::RNTupleModel::Create();
  auto ev = model->MakeField<int>("ev");
  auto info = model->MakeField<G4S_Info>("i");
  auto p = model->MakeField<std::vector<G4S_Particle>>("p");
  auto s = model->MakeField<std::vector<G4S_ParticleSteps>>("s");
  std::shared_ptr<std::vector<G4S_PathBField>> b;
  if (has_b)
    b = model->MakeField<std::vector<G4S_PathBField>>("b");

  G4S_Info *pi = info.get();
  std::vector<G4S_Particle> *pp = p.get();
  std::vector<G4S_ParticleSteps> *ps = s.get();
  std::vector<G4S_PathBField> *pb = b.get();
  if (has_ev)
    t->SetBranchAddress("ev", ev.get());
  t->SetBranchAddress("i", &pi);
  t->SetBranchAddress("p", &pp);
  t->SetBranchAddress("s", &ps);
  if (has_b)
    t->SetBranchAddress("b", &pb);

  std::unique_ptr<TFile> fout(TFile::Open(out, "RECREATE"));
  ROOT::RNTupleWriteOptions opts;
  opts.SetCompression(fin->GetCompressionSettings());
  const Long64_t n = t->GetEntries();
  TStopwatch sw;
  {
    auto w = ROOT::RNTupleWriter::Append(std::move(model), "T", *fout, opts);
    for (Long64_t i = 0; i < n; ++i) {
      t->GetEntry(i);
      if (!has_ev)
        *ev = (int)i;
      w->Fill();
      w->CommitCluster();
    }
  }
  const double t_write = sw.RealTime();
  fout->cd();
  for (const char *name : {"P", "V", "D", "C"}) {
    if (auto *tab = (TTree *)fin->Get(name)) {
      TTree *c = tab->CloneTree(-1, "fast");
      c->Write();
    }
  }
  fout->Close();

  // sizes and settings
  std::unique_ptr<TFile> fo(TFile::Open(out));
  long nclusters = 0;
  {
    auto rd = ROOT::RNTupleReader::Open("T", out);
    nclusters = rd->GetDescriptor().GetNClusters();
  }
  printf("input  %s: %lld events, TTree T zipped %.1f MB (unzipped %.1f MB), file %.1f MB, compression %d\n", in, n,
         t->GetZipBytes() / 1e6, t->GetTotBytes() / 1e6, fin->GetSize() / 1e6, fin->GetCompressionSettings());
  printf("output %s: RNTuple T in %ld clusters, file %.1f MB (tables P V D C copied as TTrees), compression %d, "
         "written in %.1f s\n",
         out, nclusters, fo->GetSize() / 1e6, opts.GetCompression(), t_write);

  // read back and compare
  auto rd = ROOT::RNTupleReader::Open("T", out);
  auto e = rd->GetModel().CreateEntry();
  auto rp = e->GetPtr<std::vector<G4S_Particle>>("p");
  auto rs = e->GetPtr<std::vector<G4S_ParticleSteps>>("s");
  std::shared_ptr<std::vector<G4S_PathBField>> rb;
  if (has_b)
    rb = e->GetPtr<std::vector<G4S_PathBField>>("b");
  int ndiff = 0;
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    rd->LoadEntry(i, *e);
    const Summary a = summarise(*pp, *ps, has_b ? pb : nullptr), c = summarise(*rp, *rs, has_b ? rb.get() : nullptr);
    if (!(a == c)) {
      ++ndiff;
      printf("event %lld differs: particles %ld/%ld steps %ld/%ld inert %ld/%ld\n", i, a.np, c.np, a.nsens, c.nsens,
             a.ninert, c.ninert);
    }
  }
  printf("read back: %lld events, %d differ in the per-event summary\n", n, ndiff);

  // time reading one event, warm
  const Long64_t mid = n / 2;
  sw.Start();
  t->GetEntry(mid);
  const double t_tree = sw.RealTime();
  sw.Start();
  rd->LoadEntry(mid, *e);
  const double t_rnt = sw.RealTime();
  printf("reading event %lld: TTree %.3f s, RNTuple %.3f s (%ld particles)\n", mid, t_tree, t_rnt, (long)rp->size());
}
