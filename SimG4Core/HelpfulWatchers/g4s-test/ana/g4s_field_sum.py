#!/usr/bin/env python3
"""g4s_field_sum.py FILE... [-- PTMIN PTMAX]: tables of g4s_field.py's output per |eta| bin: for each field model the
mean and the rms over tracks of dphi (relative error of the azimuthal bending integral, = the relative pT bias of a
fit with that model, in %), the mean dth (polar bending error over the azimuthal bending, %), and the median res
(micron at pT = 1 GeV, the azimuthal residual left after the curvature fit; it scales as 1 / pT), and the mean and
rms of dpt, the pT bias of that position fit (%)."""
import sys, math, ROOT

args = sys.argv[1:]
ptr = (0.0, 1e9)
if "--" in args:
    i = args.index("--")
    ptr = (float(args[i + 1]), float(args[i + 2]))
    args = args[:i]
ch = ROOT.TChain("F")
for a in args:
    ch.Add(a)
M = ["const", "build", "fit", "build_step", "fit_mid", "fit_s1", "fit_s2", "fit_s4"]
edges = [0.0, 0.5, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0]
rows = {i: [] for i in range(len(edges) - 1)}
for e in ch:
    if not ptr[0] <= e.pt < ptr[1]:
        continue
    a = abs(e.eta)
    for i in range(len(edges) - 1):
        if edges[i] <= a < edges[i + 1]:
            rows[i].append((e.nhits, e.ltr, e.bzmean_map, e.gth_map / e.gphi_map,
                            [getattr(e, "dphi_" + m) for m in M], [getattr(e, "dth_" + m) for m in M],
                            [getattr(e, "res_" + m) for m in M], [getattr(e, "dpt_" + m) for m in M]))


def ms(v):
    n = len(v)
    m = sum(v) / n
    return m, math.sqrt(max(0.0, sum((x - m) ** 2 for x in v) / n))


def med(v):
    v = sorted(v)
    return v[len(v) // 2]


print(f"{ch.GetEntries()} tracks, pT {ptr[0]}..{ptr[1]}")
print("dphi: mean (rms) over tracks of G_phi(model)/G_phi(map) - 1, in %")
print(f"{'|eta|':>9} {'n':>5} {'nh':>4} {'L_T':>5} {'<Bz>':>6} {'th/phi':>6} | " + " | ".join(f"{m:>15}" for m in M))
for i, r in rows.items():
    if not r:
        continue
    cells = []
    for j in range(len(M)):
        m, s = ms([x[4][j] * 100 for x in r])
        cells.append(f"{m:+7.3f} ({s:5.3f})")
    print(f"{edges[i]:.1f}-{edges[i + 1]:.1f} {len(r):5d} {med([x[0] for x in r]):4.0f} {med([x[1] for x in r]):5.0f} "
          f"{ms([x[2] for x in r])[0]:6.3f} {ms([x[3] * 100 for x in r])[0]:+6.2f} | " + " | ".join(cells))
print("dpt: mean (rms) over tracks of the pT bias of a position fit (1/c - 1), in %")
print(f"{'|eta|':>9} | " + " | ".join(f"{m:>15}" for m in M))
for i, r in rows.items():
    if not r:
        continue
    cells = []
    for j in range(len(M)):
        m, s = ms([x[7][j] * 100 for x in r])
        cells.append(f"{m:+7.3f} ({s:5.3f})")
    print(f"{edges[i]:.1f}-{edges[i + 1]:.1f} | " + " | ".join(cells))
print("dth: mean of (G_th(model) - G_th(map)) / G_phi(map), in %;  res: median micron at pT = 1 GeV")
print(f"{'|eta|':>9} | " + " | ".join(f"{m:>17}" for m in M))
for i, r in rows.items():
    if not r:
        continue
    cells = []
    for j in range(len(M)):
        cells.append(f"{ms([x[5][j] * 100 for x in r])[0]:+7.3f} {med([x[6][j] for x in r]):7.1f}")
    print(f"{edges[i]:.1f}-{edges[i + 1]:.1f} | " + " | ".join(cells))
