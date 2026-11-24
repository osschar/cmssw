"""g4s_read.py: a snitch output file, its events T a TTree or an RNTuple, read the same way.

    from g4s_read import SnitchFile
    sf = SnitchFile("G4Snitch.root")
    P = [str(x.name) for x in sf.f.Get("P")]   # the tables P, V, D, C are TTrees in both formats
    for ev in sf:                              # or sf.entry(k)
        ev.ev, ev.i, ev.p, ev.s, ev.b          # as the TTree's branches: one object, reloaded for each event

An RNTuple's entries are in the order the threads finished their events: match events by ev. The objects of an
entry are valid until the next entry is loaded, as with a TTree. Needs libRootG4Snitch.so (loaded here)."""
import ROOT

ROOT.gSystem.Load("libRootG4Snitch.so")


class _Entry:
    """An RNTuple entry with the fields as attributes, like a TTree's branches."""
    __slots__ = ("_e",)

    def __init__(self, e):
        self._e = e

    def __getattr__(self, name):
        return self._e[name]


class SnitchFile:
    def __init__(self, fn):
        self.fn = fn
        self.f = ROOT.TFile.Open(fn)
        if not self.f or self.f.IsZombie():
            raise OSError(f"cannot open {fn}")
        t = self.f.Get("T")
        self.is_ttree = isinstance(t, ROOT.TTree)
        if self.is_ttree:
            self.t = t
            self.n = t.GetEntries()
            self.fields = {str(b.GetName()) for b in t.GetListOfBranches()}
        else:
            self.r = ROOT.RNTupleReader.Open("T", fn)
            self.e = self.r.CreateEntry()
            self.n = self.r.GetNEntries()
            self.fields = {str(fld.GetFieldName()) for fld in self.r.GetDescriptor().GetTopLevelFields()}
            self._entry = _Entry(self.e)

    def __len__(self):
        return self.n

    def entry(self, k):
        """Load entry k and return it (the same object for every entry, as a TTree)."""
        if self.is_ttree:
            self.t.GetEntry(k)
            return self.t
        self.r.LoadEntry(k, self.e)
        return self._entry

    def __iter__(self):
        for k in range(self.n):
            yield self.entry(k)
