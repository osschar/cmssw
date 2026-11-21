#!/usr/bin/env python3
"""g4s_hits_check.py STEP1_A STEP1_B [G4SNITCH_B]: in the cmsRun area's environment.
1. Every PCaloHit and PSimHit collection of the two event files must be equal, hit by hit (id or detUnitId, energy,
   time, track, and for PSimHits the entry and exit points).
2. With a snitch file (written in the same job as STEP1_B): how many of its sensitive steps with a deposit and a
   DetId have that DetId among the hits of their event, by DetId::Detector."""
import sys, ROOT
ROOT.gSystem.Load("libFWCoreFWLite"); ROOT.FWLiteEnabler.enable()
for lib in ("libDataFormatsCommon", "libSimDataFormatsCaloHit", "libSimDataFormatsTrackingHit"):
    ROOT.gSystem.Load(lib)
from DataFormats.FWLite import Events, Handle

def collections(fn):
    f = ROOT.TFile.Open(fn); t = f.Get("Events"); out = []
    for b in t.GetListOfBranches():
        n = b.GetName().rstrip(".")
        if n.startswith("PCaloHits_") or n.startswith("PSimHits_"):
            typ, mod, inst, proc = n.split("_")
            cls = b.GetClassName()  # edm::Wrapper<vector<io_v1::PCaloHit> >: the handle takes the wrapped type
            out.append((cls[len("edm::Wrapper<"):cls.rindex(">")].strip(), mod, inst, proc))
    f.Close(); return sorted(out)

def key(h, calo):
    if calo:
        return (h.id(), h.energy(), h.time(), h.geantTrackId(), h.depth())
    e, x = h.entryPoint(), h.exitPoint()
    return (h.detUnitId(), h.energyLoss(), h.tof(), h.trackId(), e.x(), e.y(), e.z(), x.x(), x.y(), x.z())

ca, cb = collections(sys.argv[1]), collections(sys.argv[2])
assert ca == cb, "different collections"
ea, eb = Events(sys.argv[1]), Events(sys.argv[2])
handles = [(c, Handle(c[0])) for c in ca]
n_hits, n_diff, ids = 0, 0, []
for eva, evb in zip(ea, eb):
    ev_ids = set()
    for c, h in handles:
        calo = "PCaloHit" in c[0]
        eva.getByLabel(c[1], c[2], c[3], h); la = [key(x, calo) for x in h.product()] if h.isValid() else []
        evb.getByLabel(c[1], c[2], c[3], h); lb = [key(x, calo) for x in h.product()] if h.isValid() else []
        n_hits += len(la)
        if la != lb:
            n_diff += 1; print("DIFFERENT:", c[1:], len(la), len(lb))
        ev_ids.update(k[0] for k in lb)
    ids.append(ev_ids)
print(f"{len(ids)} events, {len(handles)} hit collections, {n_hits} hits in A; collections differing: {n_diff}")

if len(sys.argv) > 3:
    ROOT.gSystem.Load("libRootG4Snitch.so")
    f = ROOT.TFile.Open(sys.argv[3]); t = f.Get("T")
    vol_sd = [int(v.sd) for v in f.Get("V")]
    dt = f.Get("D"); sd_name = {}
    for d in dt:
        sd_name[int(d.idx)] = f"{str(getattr(d, 'class'))} ({str(d.name)})"
    tot, found, zero = {}, {}, {}
    for i, ev in enumerate(t):
        for ps in ev.s:
            for st in ps.m_steps_sensitive:
                if st.m_edep.m_total <= 0: continue
                k = sd_name.get(vol_sd[st.m_volume], "none")
                if k.startswith("HCalSD"):  # HCAL test numbering: the top four bits are the HCAL part
                    k += " " + ({1: "HB", 2: "HE", 3: "HO", 4: "HF"}.get(st.m_det_id >> 28, "?") if st.m_det_id else "(no id)")
                if st.m_det_id == 0: zero[k] = zero.get(k, 0) + 1; continue
                tot[k] = tot.get(k, 0) + 1
                if st.m_det_id in ids[i]: found[k] = found.get(k, 0) + 1
    print("snitch sensitive steps with a deposit, by sensitive detector (class (name)):")
    for k in sorted(set(tot) | set(zero)):
        n = tot.get(k, 0)
        print(f"  {k:58s} id 0: {zero.get(k, 0):6d}  with id: {n:7d}  among the event's hits: {found.get(k, 0):7d}" + (f" ({100.0 * found.get(k, 0) / n:.2f} %)" if n else ""))
