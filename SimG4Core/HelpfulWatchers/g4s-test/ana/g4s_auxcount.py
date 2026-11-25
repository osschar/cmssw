#!/usr/bin/env python3
"""g4s_auxcount.py FILE...: what auxiliary points on curved steps would add to the snitch's paths, counted from
existing files (needs the field branch "b"); and how many of the large turns between path points are field turns.

For each segment between consecutive path points of a charged particle: the path length L = beta c dt from the
stored times (beta = |p| / E, the mean of the two ends), the field turn phi_f = k |q| |B x t| L / p (k = 0.0029979
GeV / (T cm), B the mean of the stored field at the two ends, t and p from the mean of the two momenta), and the
observed turn phi_o, the angle between the two momenta. Auxiliary points every TOL degrees of field turn would add
ceil(phi_f / TOL) - 1 points to the segment (counted for TOL 10, 20, 30). A segment that turns much in total turned
in its last Geant4 step (the thinning stores a point after 0.05 rad), so phi_f approximates that step's field turn.
A segment with phi_o > 10 degrees is called field-like if |phi_o - phi_f| < 0.2 phi_f + 2 degrees.
Classes by |PDG|: e (11), mu (13), hadron (211, 321, 2212), other charged; neutral particles have path points but
no turns. For e+- at 20 degrees also: the extra points by their number per segment, by region, by kinetic energy,
and the totals with a cap of 17 (one turn) or 35 extra points per segment."""
import sys, math, collections, ROOT
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
ROOT.gSystem.Load("libRootG4Snitch.so")

K = 0.299792458e-2
C = 29.9792458  # cm / ns
TOLS = (10, 20, 30)


def cls(pdg, q):
    a = abs(pdg)
    if q == 0:
        return "neutral"
    return {11: "e", 13: "mu", 211: "hadron", 321: "hadron", 2212: "hadron"}.get(a, "other")


nev = 0
pts = collections.Counter()
segs = collections.Counter()
extra = {t: collections.Counter() for t in TOLS}
big = collections.Counter()      # segments with phi_o > 10 deg
fieldlike = collections.Counter()
lowfield = collections.Counter()  # phi_o > 10 deg with phi_f < 2 deg: not the field
maxextra = collections.Counter()
# e+- only, at 20 degrees: extra points by points per segment, by region, by kinetic energy; and with a cap of
# CAP extra points per segment (one turn = 17)
SIZE = [(1, 2), (3, 17), (18, 100), (101, 10 ** 9)]
bysize = collections.Counter()
byreg = collections.Counter()
byek = collections.Counter()
CAPS = (17, 35)
capped = collections.Counter()
EK = [0.001, 0.01, 0.1, 1.0, 1e9]  # GeV
for fn in sys.argv[1:]:
    sf = SnitchFile(fn); f = sf.f
    for ev in sf:
        nev += 1
        p, s, bb = ev.p, ev.s, ev.b
        for k in range(p.size()):
            q = p[k].m_charge
            c = cls(p[k].m_pdg, q)
            path = s[k].m_path
            n = path.size()
            pts[c] += n
            if q == 0 or n < 2:
                continue
            bv = bb[k].m_b
            if bv.size() != n:
                continue
            for i in range(n - 1):
                a, b = path[i], path[i + 1]
                pa = (a.m_p.fX, a.m_p.fY, a.m_p.fZ)
                pb = (b.m_p.fX, b.m_p.fY, b.m_p.fZ)
                na = math.sqrt(sum(x * x for x in pa))
                nb = math.sqrt(sum(x * x for x in pb))
                if na <= 0 or nb <= 0:
                    continue
                segs[c] += 1
                cosang = sum(x * y for x, y in zip(pa, pb)) / (na * nb)
                phio = math.degrees(math.acos(max(-1.0, min(1.0, cosang))))
                pm = [0.5 * (x + y) for x, y in zip(pa, pb)]
                pmn = math.sqrt(sum(x * x for x in pm))
                if pmn <= 0:
                    continue
                t = [x / pmn for x in pm]
                beta = 0.5 * (na / a.m_p.fT + nb / b.m_p.fT) if a.m_p.fT > 0 and b.m_p.fT > 0 else 1.0
                L = beta * C * abs(b.m_x.fT - a.m_x.fT)
                B = [0.5 * (bv[i].fX + bv[i + 1].fX), 0.5 * (bv[i].fY + bv[i + 1].fY), 0.5 * (bv[i].fZ + bv[i + 1].fZ)]
                bxt = (B[1] * t[2] - B[2] * t[1], B[2] * t[0] - B[0] * t[2], B[0] * t[1] - B[1] * t[0])
                bperp = math.sqrt(sum(x * x for x in bxt))
                phif = math.degrees(K * abs(q) * bperp * L / (0.5 * (na + nb)))
                for tol in TOLS:
                    if phif > tol:
                        e = int(math.ceil(phif / tol)) - 1
                        extra[tol][c] += e
                        if tol == 20:
                            maxextra[c] = max(maxextra[c], e)
                            for cap in CAPS:
                                capped[(c, cap)] += min(e, cap)
                            if c == "e":
                                bysize[next(i for i, (lo, hi) in enumerate(SIZE) if lo <= e <= hi)] += e
                                xm = [0.5 * (a.m_x.fX + b.m_x.fX), 0.5 * (a.m_x.fY + b.m_x.fY), 0.5 * (a.m_x.fZ + b.m_x.fZ)]
                                r = math.hypot(xm[0], xm[1])
                                reg = ("beam pipe r<2.2" if r < 2.2 and abs(xm[2]) < 300 else
                                       "tracker r<118 |z|<282" if r < 118 and abs(xm[2]) < 282 else
                                       "r<300 |z|<600 (calo, HGCal)" if r < 300 and abs(xm[2]) < 600 else "outside")
                                byreg[reg] += e
                                ek = 0.5 * (na + nb)
                                ek = math.sqrt(ek * ek + 0.000511 ** 2) - 0.000511
                                byek[next(i for i, x in enumerate(EK) if ek < x)] += e
                if phio > 10:
                    big[c] += 1
                    if abs(phio - phif) < 0.2 * phif + 2:
                        fieldlike[c] += 1
                    elif phif < 2:
                        lowfield[c] += 1
    f.Close()

classes = ["e", "mu", "hadron", "other", "neutral"]
tot = sum(pts.values())
print(f"{nev} events; path points per event today: {tot / nev:.0f}")
print(f"{'class':>8} {'points/ev':>10} {'segments/ev':>11} " + " ".join(f"{'+aux ' + str(t) + 'deg':>11}" for t in TOLS)
      + f" {'(+% of all pts at 20)':>21} {'max/seg 20':>10} {'turn>10deg/ev':>13} {'field-like':>10} {'not field':>9}")
for c in classes:
    if not pts[c]:
        continue
    row = f"{c:>8} {pts[c] / nev:10.1f} {segs[c] / nev:11.1f} " + " ".join(f"{extra[t][c] / nev:11.1f}" for t in TOLS)
    row += f" {100 * extra[20][c] / tot:20.2f}% {maxextra[c]:10d} {big[c] / nev:13.1f}"
    row += f" {100 * fieldlike[c] / big[c]:9.0f}% {100 * lowfield[c] / big[c]:8.0f}%" if big[c] else ""
    print(row)
print(f"{'all':>8} {tot / nev:10.1f} {sum(segs.values()) / nev:11.1f} " + " ".join(f"{sum(extra[t].values()) / nev:11.1f}" for t in TOLS)
      + f" {100 * sum(extra[20].values()) / tot:20.2f}%")

ne = sum(bysize.values())
if ne:
    print(f"e+- extra points at 20 deg ({ne / nev:.0f} per event): by extra points per segment " +
          ", ".join(f"{lo}-{hi if hi < 10 ** 8 else 'inf'}: {100 * bysize[i] / ne:.0f}%" for i, (lo, hi) in enumerate(SIZE)))
    print("  by region: " + ", ".join(f"{k}: {100 * v / ne:.0f}%" for k, v in byreg.most_common()))
    lo = [0.0] + EK[:-1]
    print("  by kinetic energy: " + ", ".join(f"{1000 * lo[i]:g}-{1000 * EK[i]:g} MeV: {100 * byek[i] / ne:.0f}%"
                                              for i in range(len(EK))))
for cap in CAPS:
    print(f"with at most {cap} extra points per segment: " + ", ".join(
        f"{c} {capped[(c, cap)] / nev:.0f}" for c in classes if capped[(c, cap)]) +
          f"; all {sum(v for (c, k), v in capped.items() if k == cap) / nev:.0f} per event "
          f"(+{100 * sum(v for (c, k), v in capped.items() if k == cap) / tot:.1f}% of the points)")
