#!/usr/bin/env python3
"""g4s_hits_check.py A.root B.root [B_G4Snitch.root]: in the cmsRun area's environment.
1. Every PCaloHit, PSimHit, SimTrack and SimVertex collection of g4SimHits in the event files A and B must be equal,
   item by item (ids, energies, times, tracks, points, momenta, vertices), for the events present in both (matched by
   run, lumi and event; the process name is ignored, so a re-run of a few events compares against the original).
2. With a snitch file written in the same job as B: how many of its sensitive steps with a deposit and a DetId have
   that DetId among the hits of their event, by sensitive detector (HCAL by part: its ids are test-numbering indices)."""
import sys, ROOT
ROOT.gSystem.Load("libFWCoreFWLite"); ROOT.FWLiteEnabler.enable()
for lib in ("libDataFormatsCommon", "libSimDataFormatsCaloHit", "libSimDataFormatsTrackingHit"):
    ROOT.gSystem.Load(lib)
from DataFormats.FWLite import Events, Handle

KINDS = ("PCaloHits", "PSimHits", "SimTracks", "SimVertexs")

def collections(fn):
    """{(kind, module, instance): (wrapped type, process)} of g4SimHits."""
    f = ROOT.TFile.Open(fn); t = f.Get("Events"); out = {}
    for b in t.GetListOfBranches():
        n = b.GetName().rstrip(".")
        parts = n.split("_")
        if len(parts) != 4 or parts[0] not in KINDS or parts[1] != "g4SimHits":
            continue
        cls = b.GetClassName()  # edm::Wrapper<vector<io_v1::PCaloHit> >: the handle takes the wrapped type
        out[(parts[0], parts[1], parts[2])] = (cls[len("edm::Wrapper<"):cls.rindex(">")].strip(), parts[3])
    f.Close()
    return out

def key(kind, h):
    if kind == "PCaloHits":
        return (h.id(), h.energy(), h.time(), h.geantTrackId(), h.depth())
    if kind == "PSimHits":
        e, x = h.entryPoint(), h.exitPoint()
        return (h.detUnitId(), h.energyLoss(), h.tof(), h.trackId(), e.x(), e.y(), e.z(), x.x(), x.y(), x.z())
    if kind == "SimTracks":
        m = h.momentum()
        return (h.trackId(), h.type(), m.x(), m.y(), m.z(), m.t(), h.vertIndex(), h.genpartIndex())
    p = h.position()
    return (p.x(), p.y(), p.z(), p.t(), h.parentIndex(), h.vertexId())

def event_ids(evs):
    return [(e.eventAuxiliary().id().run(), e.eventAuxiliary().id().luminosityBlock(), e.eventAuxiliary().id().event())
            for e in evs]

ca, cb = collections(sys.argv[1]), collections(sys.argv[2])
common = sorted(set(ca) & set(cb))
only = sorted(set(ca) ^ set(cb))
ea, eb = Events(sys.argv[1]), Events(sys.argv[2])
ida, idb = event_ids(ea), event_ids(eb)
ia = {x: i for i, x in enumerate(ida)}
handles = {c: (Handle(ca[c][0]), Handle(cb[c][0])) for c in common}
n_items, n_diff, n_ev, ids_b = 0, 0, 0, []
for jb, x in enumerate(idb):
    ev_ids = set()
    eb.to(jb)
    ja = ia.get(x)
    if ja is not None:
        ea.to(ja); n_ev += 1
    for c in common:
        hA, hB = handles[c]
        eb.getByLabel(c[1], c[2], cb[c][1], hB); lb = [key(c[0], y) for y in hB.product()] if hB.isValid() else []
        if c[0] in ("PCaloHits", "PSimHits"):
            ev_ids.update(k[0] for k in lb)
        if ja is None:
            continue
        ea.getByLabel(c[1], c[2], ca[c][1], hA); la = [key(c[0], y) for y in hA.product()] if hA.isValid() else []
        n_items += len(la)
        if la != lb:
            n_diff += 1; print("DIFFERENT:", x, c, len(la), len(lb))
    ids_b.append(ev_ids)
print(f"{n_ev} events in both, {len(common)} collections (hits, SimTracks, SimVertices), {n_items} items in A; "
      f"collections differing: {n_diff}" + (f"; in only one file: {only}" if only else ""))

if len(sys.argv) > 3:
    from g4s_read import SnitchFile
    sf = SnitchFile(sys.argv[3]); f = sf.f
    jb_of_event = {x[2]: jb for jb, x in enumerate(idb)}  # the snitch's ev is the CMSSW event number
    vol_sd = [int(v.sd) for v in f.Get("V")]
    sd_name = {int(d.idx): f"{str(getattr(d, 'class'))} ({str(d.name)})" for d in f.Get("D")}
    tot, found, zero = {}, {}, {}
    for k, ev in enumerate(sf):
        i = jb_of_event[int(ev.ev)] if "ev" in sf.fields else k
        for ps in ev.s:
            for st in ps.m_steps_sensitive:
                if st.m_edep.m_total <= 0: continue
                k = sd_name.get(vol_sd[st.m_volume], "none")
                if k.startswith("HCalSD"):  # HCAL test numbering: the top four bits are the HCAL part
                    k += " " + ({1: "HB", 2: "HE", 3: "HO", 4: "HF"}.get(st.m_det_id >> 28, "?") if st.m_det_id else "(no id)")
                if st.m_det_id == 0: zero[k] = zero.get(k, 0) + 1; continue
                tot[k] = tot.get(k, 0) + 1
                if st.m_det_id in ids_b[i]: found[k] = found.get(k, 0) + 1
    print("snitch sensitive steps with a deposit, by sensitive detector (class (name)):")
    for k in sorted(set(tot) | set(zero)):
        n = tot.get(k, 0)
        print(f"  {k:58s} id 0: {zero.get(k, 0):6d}  with id: {n:7d}  among the event's hits: {found.get(k, 0):7d}" + (f" ({100.0 * found.get(k, 0) / n:.2f} %)" if n else ""))
