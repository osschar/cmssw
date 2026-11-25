#!/usr/bin/env python3
"""g4s_disc_split.py GEO ENVELOPES FILE...: where a forward pixel disc's material lies relative to its sensors, from the
primary's true path. For each step between hits on adjacent discs of one family (TFPX, mkFit layers 16-23 and 38-45;
TEPX, 24-27 and 46-49; a hit's disc from mkFit's module dump GEO), and between hits on adjacent barrel layers (IT1-IT4,
OT1-OT6, the two sensors of an OT stack being one layer), the true X0 between the two hits' mid points is
split by position into: inside the source disc's envelope (after its hit), inside the destination's (before its
hit), and elsewhere. ENVELOPES as for g4s_tgeo_cmp.py (the TGeo material session's discs at zpad 0.5, barrel layers at bare extent). Each Geant4 step
inside the tracker volume (r < 120 cm, |z| < 290 cm at either end) is cut into pieces of at most 0.1 cm, timed by
linear interpolation. A line "SUMS family n after before elsewhere" per family lets runs on several files be added. The share after the hit, src / (src + dst), is
what "half before, half after" assumes to be 0.5."""
import sys, math, collections, ROOT
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
ROOT.gSystem.Load("libRootG4Snitch.so")

geo = ROOT.TFile.Open(sys.argv[1])
mods = {int(m.detid): int(m.layer) for m in geo.Get("M")}
ENV = {}
for line in open(sys.argv[2]):
    a = line.split()
    if a[0] == "D" and (a[1].startswith("TFPX") or a[1].startswith("TEPX")):
        z, hz = float(a[4]), float(a[5]) + 0.5
        ENV[a[1]] = (float(a[2]), float(a[3]), z - hz, z + hz)
    elif a[0] == "L" and (a[1].startswith("IT") or a[1].startswith("OT")):
        ENV[a[1]] = (float(a[2]), float(a[3]), 0.0, float(a[4]))


BARREL = ["IT1", "IT2", "IT3", "IT4", "OT1", "OT2", "OT3", "OT4", "OT5", "OT6"]


def disc(lay):
    """(family, index) of an mkFit pixel-disc layer, both sides folded; ("BRL", i) for the barrel, i indexing BARREL"""
    if 0 <= lay <= 3:
        return ("BRL", lay)
    if 4 <= lay <= 15:
        return ("BRL", 4 + (lay - 4) // 2)
    l = lay if lay < 38 else lay - 22
    if 16 <= l <= 23:
        return ("TFPX", l - 16 + 1)
    if 24 <= l <= 27:
        return ("TEPX", l - 24 + 9)
    return None


def env_of(r, z):
    for nm, (r0, r1, z0, z1) in ENV.items():
        if r0 <= r <= r1 and z0 <= abs(z) <= z1:
            return nm
    return None


acc = collections.defaultdict(collections.Counter)
nstep = collections.Counter()
for fn in sys.argv[3:]:
    sf = SnitchFile(fn)
    sdcls = [str(getattr(d, "class")) for d in sf.rows("D")]
    trk = {i for i, c in enumerate(sdcls) if "TkAccumulating" in c}
    V = [(float(v.x0), int(v.sd) in trk) for v in sf.rows("V")]
    for ev in sf:
        p, s = ev.p, ev.s
        k = p[0].m_daughters_begin
        q = p[k]
        g = q.m_g4_id
        own = sorted([st for st in list(s[k].m_steps_sensitive) + list(s[k].m_steps_inert) if st.m_g4_pid == g],
                     key=lambda st: st.m_x_beg.fT)
        pieces, hits = [], []
        for st in own:
            a = (st.m_x_beg.fX, st.m_x_beg.fY, st.m_x_beg.fZ)
            b = (st.m_x_end.fX, st.m_x_end.fY, st.m_x_end.fZ)
            L = math.dist(a, b)
            x0 = V[st.m_volume][0] if 0 <= st.m_volume < len(V) else 0
            in_trk = (math.hypot(a[0], a[1]) < 120 and abs(a[2]) < 290) or (math.hypot(b[0], b[1]) < 120 and abs(b[2]) < 290)
            if L > 0 and x0 > 0 and in_trk:
                nn = max(1, int(math.ceil(L / 0.1)))
                for j in range(nn):
                    u = (j + 0.5) / nn
                    pt_ = [a[c] + u * (b[c] - a[c]) for c in range(3)]
                    pieces.append((st.m_x_beg.fT + u * (st.m_x_end.fT - st.m_x_beg.fT), math.hypot(pt_[0], pt_[1]),
                                   pt_[2], L / nn / x0))
            if st.m_det_id != 0 and 0 <= st.m_volume < len(V) and V[st.m_volume][1]:
                if hits and hits[-1][2] == st.m_det_id:
                    continue
                d = disc(mods.get(int(st.m_det_id), -1))
                hits.append((0.5 * (st.m_x_beg.fT + st.m_x_end.fT), d, st.m_det_id))
        for (ta, da, _), (tb, db, _) in zip(hits, hits[1:]):
            if da is None or db is None or da[0] != db[0] or db[1] != da[1] + 1:
                continue
            fam = da[0]
            if fam == "BRL":
                src, dst = BARREL[da[1]], BARREL[db[1]]
                fam = f"BRL {src}->{dst}"
            else:
                src, dst = f"{fam}{da[1]}", f"{fam}{db[1]}"
            nstep[fam] += 1
            for (tt, r, z, x) in pieces:
                if ta <= tt <= tb:
                    e = env_of(r, z)
                    acc[fam]["after the source hit" if e == src else "before the destination hit" if e == dst
                             else "elsewhere"] += x
for fam in sorted(acc):
    a = acc[fam]
    n = nstep[fam]
    s_, d_ = a["after the source hit"], a["before the destination hit"]
    print(f"SUMS {fam} {n} {s_:.6f} {d_:.6f} {a['elsewhere']:.6f}")
    print(f"{fam}: {n} steps between hits on adjacent discs; X0 per step: after the source hit {s_ / n:.4f}, before the "
          f"destination hit {d_ / n:.4f}, elsewhere {a['elsewhere'] / n:.4f}; share of a disc's material after its "
          f"hit {s_ / (s_ + d_):.2f}")
