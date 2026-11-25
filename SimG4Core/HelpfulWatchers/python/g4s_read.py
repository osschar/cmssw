"""g4s_read.py: a G4Snitch output file, read the same way whatever its format: RNTuples (since 2026-11-24) or the
TTrees of older files.

    from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
    sf = SnitchFile("G4Snitch.root")
    for ev in sf:                              # or sf.entry(k); len(sf) events
        ev.ev, ev.i, ev.p, ev.s, ev.b          # one object, reloaded for each event, as a TTree's branches
    P = [str(x.name) for x in sf.rows("P")]   # the tables P, V, D, C, row by row, the same way
    C = sf.config()                            # tree C as a dict

T's entries are in the order the threads finished their events: match events by ev. The objects of an entry or a
row are valid until the next one is loaded. Needs libRootG4Snitch.so (loaded here)."""
import ROOT

ROOT.gSystem.Load("libRootG4Snitch.so")


class _Entry:
    """An RNTuple entry with the fields as attributes, like a TTree's branches."""
    __slots__ = ("_e",)

    def __init__(self, e):
        self._e = e

    def __getattr__(self, name):
        return self._e[name]


def _rntuple_rows(name, fn):
    r = ROOT.RNTupleReader.Open(name, fn)
    e = r.CreateEntry()
    w = _Entry(e)
    for k in range(r.GetNEntries()):
        r.LoadEntry(k, e)
        yield w


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
        """Load event k and return it (the same object for every event, as a TTree)."""
        if self.is_ttree:
            self.t.GetEntry(k)
            return self.t
        self.r.LoadEntry(k, self.e)
        return self._entry

    def __iter__(self):
        for k in range(self.n):
            yield self.entry(k)

    def rows(self, name):
        """The rows of table name (P, V, D, C), each with the columns as attributes (getattr(x, "class") for D's)."""
        t = self.f.Get(name)
        if not t:
            return iter(())
        if isinstance(t, ROOT.TTree):
            return iter(t)
        return _rntuple_rows(name, self.fn)

    def config(self):
        """Table C as a dict, key -> value."""
        return {str(x.key): str(x.value) for x in self.rows("C")}
