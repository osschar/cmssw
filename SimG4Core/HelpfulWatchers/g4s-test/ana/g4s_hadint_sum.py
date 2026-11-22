#!/usr/bin/env python3
"""g4s_hadint_sum.py OUT...: tables of g4s_hadint.py's output, per pT and |eta| region.
1. Fates inside the tracker: inelastic before the first hit, inelastic after at least one hit, decay; for the
   inelastic ones with hits, the median number of hits before and its share of the median of the survivors in the same
   bin, and the share with at least 4 hits (seeded, then stopped).
2. The leading charged daughter of an inelastic interaction after at least one hit: share above 0.3 GeV; share with
   x = p_lead / p_before above 0.5 and 0.8; share with the angle below 10 and 50 mrad; "kink-like" (x > 0.5 and angle
   < 50 mrad); share leaving at least 4 tracker hits of its own; charged daughters above 0.3 GeV.
3. Decays inside (pi -> mu nu): median x, median angle, median muon hits.
4. Kinks from hadElastic and CoulombScat along the primary between its first and last hit: share of tracks with at
   least one above 20 and 100 mrad.
5. Where the inelastic interactions after the first hit happen: radius (|eta| < 1.2) and |z| (|eta| > 1.7)."""
import sys, math, collections, ROOT

H = ROOT.TChain("H")
for a in sys.argv[1:]:
    H.Add(a)
REG = [(0.0, 0.9), (0.9, 1.5), (1.5, 2.5), (2.5, 3.0)]
rows = collections.defaultdict(list)
for e in H:
    r = next((i for i, (lo, hi) in enumerate(REG) if lo <= abs(e.eta) < hi), None)
    if r is None:
        continue
    rows[(round(e.pt), r)].append({c: getattr(e, c) for c in ("eta", "nhits", "fate", "death_r", "death_z", "p_before",
                                                              "lead_pdg", "lead_x", "lead_ang", "lead_hits", "n_ch03",
                                                              "nk20", "nk100")})


def med(v):
    v = sorted(v)
    return v[len(v) // 2] if v else float("nan")


def pct(n, d):
    return f"{100 * n / d:5.1f}" if d else "    -"


pts = sorted({k[0] for k in rows})
keys = [(pt, r) for pt in pts for r in range(len(REG))] + [(None, r) for r in range(len(REG))]


def sel(key):
    pt, r = key
    return [x for k, v in rows.items() if k[1] == r and (pt is None or k[0] == pt) for x in v]


def lab(key):
    pt, r = key
    return f"{'all' if pt is None else pt:>4} {REG[r][0]:.1f}-{REG[r][1]:.1f}"


print("1. fates inside the tracker (r < 118 cm, |z| < 282 cm, outgoing leg), %")
print(f"{'pT':>4} {'|eta|':>7} {'n':>5} {'inel, 0 hits':>12} {'inel, >=1 hit':>13} {'decay':>6} {'hits before (med)':>17} "
      f"{'of survivors':>12} {'>=4 hits':>8}")
for key in keys:
    v = sel(key)
    if not v:
        continue
    n = len(v)
    i0 = [x for x in v if x["fate"] == 1 and x["nhits"] == 0]
    i1 = [x for x in v if x["fate"] == 1 and x["nhits"] > 0]
    de = [x for x in v if x["fate"] == 2]
    surv = med([x["nhits"] for x in v if x["fate"] == 0])
    mh = med([x["nhits"] for x in i1])
    print(f"{lab(key)} {n:5d} {pct(len(i0), n):>12} {pct(len(i1), n):>13} {pct(len(de), n):>6} {mh:17.0f} "
          f"{(mh / surv if surv else float('nan')):12.2f} {pct(sum(1 for x in i1 if x['nhits'] >= 4), len(i1)):>8}")

print("\n2. leading charged daughter of an inelastic interaction after >= 1 hit, % of those interactions")
print(f"{'pT':>4} {'|eta|':>7} {'n':>5} {'p>0.3':>6} {'x>0.5':>6} {'x>0.8':>6} {'<10mrad':>7} {'<50mrad':>7} {'kink-like':>9} "
      f"{'>=4 hits':>8} {'nch>0.3 (med)':>13}")
for key in keys:
    v = [x for x in sel(key) if x["fate"] == 1 and x["nhits"] > 0]
    if not v:
        continue
    n = len(v)
    lead = [x for x in v if x["lead_x"] > 0 and x["lead_x"] * x["p_before"] > 0.3]
    print(f"{lab(key)} {n:5d} {pct(len(lead), n):>6} {pct(sum(1 for x in v if x['lead_x'] > 0.5), n):>6} "
          f"{pct(sum(1 for x in v if x['lead_x'] > 0.8), n):>6} "
          f"{pct(sum(1 for x in v if 0 <= x['lead_ang'] < 0.010), n):>7} "
          f"{pct(sum(1 for x in v if 0 <= x['lead_ang'] < 0.050), n):>7} "
          f"{pct(sum(1 for x in v if x['lead_x'] > 0.5 and 0 <= x['lead_ang'] < 0.050), n):>9} "
          f"{pct(sum(1 for x in v if x['lead_hits'] >= 4), n):>8} {med([x['n_ch03'] for x in v]):13.0f}")

print("\n3. decays inside the tracker: median x, angle (mrad), muon hits")
for key in keys:
    v = [x for x in sel(key) if x["fate"] == 2 and x["lead_x"] > 0]
    if v:
        print(f"{lab(key)} {len(v):5d}  x {med([x['lead_x'] for x in v]):.2f}  angle {1000 * med([x['lead_ang'] for x in v]):6.1f}"
              f"  hits {med([x['lead_hits'] for x in v]):.0f}")

print("\n4. kinks (hadElastic, CoulombScat) between the first and the last hit, % of tracks with >= 2 hits")
print(f"{'pT':>4} {'|eta|':>7} {'n':>5} {'>20 mrad':>8} {'>100 mrad':>9}")
for key in keys:
    v = [x for x in sel(key) if x["nhits"] >= 2]
    if v:
        print(f"{lab(key)} {len(v):5d} {pct(sum(1 for x in v if x['nk20'] > 0), len(v)):>8} "
              f"{pct(sum(1 for x in v if x['nk100'] > 0), len(v)):>9}")

print("\n5. where inelastic interactions after >= 1 hit happen (all pT), % of them")
v = [x for k, vv in rows.items() for x in vv if x["fate"] == 1 and x["nhits"] > 0]
rb = [0, 2.5, 17.0, 21.5, 33.0, 60.0, 118.0]
vb = [x for x in v if abs(x["eta"]) < 1.2]
print(f"  |eta| < 1.2, {len(vb)}: " + ", ".join(
    f"r {rb[i]:g}-{rb[i + 1]:g}: {pct(sum(1 for x in vb if rb[i] <= x['death_r'] < rb[i + 1]), len(vb)).strip()}"
    for i in range(len(rb) - 1)))
zb = [0, 25, 60, 120, 170, 282]
ve = [x for x in v if abs(x["eta"]) > 1.7]
print(f"  |eta| > 1.7, {len(ve)}: " + ", ".join(
    f"|z| {zb[i]:g}-{zb[i + 1]:g}: {pct(sum(1 for x in ve if zb[i] <= abs(x['death_z']) < zb[i + 1]), len(ve)).strip()}"
    for i in range(len(zb) - 1)))
