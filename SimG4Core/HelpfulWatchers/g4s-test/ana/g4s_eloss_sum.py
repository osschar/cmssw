#!/usr/bin/env python3
"""g4s_eloss_sum.py OUT...: tables of g4s_eloss.py's output.
1. Per pT and |eta| bin, per track between the first and the last hit: the mean true energy loss, mkFit's (module
   bbxi), their ratio, the ratio with the PS/2S stack's value applied once (no model loss on the step between the
   two sensors of one stack), and the true loss over p.
2. Per arrival at a physical layer (pT up to --ptmax, all by default): mean and median of the true loss against mkFit's mean, and the ratio of
   the means; the median shows how far the most probable loss lies below the mean (Landau tail)."""
import sys, math, collections, ROOT


def phys(lay):
    if lay < 0:
        return None
    if lay <= 3 or 16 <= lay <= 27 or 38 <= lay <= 49:
        return lay
    if lay <= 15:
        return 100 + (lay - 4) // 2
    return 200 + (lay - 28) // 2 if lay <= 37 else 300 + (lay - 50) // 2


def name(ph):
    if ph <= 3:
        return f"IT L{ph + 1}"
    if ph < 100:
        l = ph if ph < 38 else ph - 22
        return f"TFPX D{l - 15}" if l <= 23 else f"TEPX D{l - 23}"
    if ph < 200:
        return f"OT L{ph - 99}"
    return f"TEDD D{(ph % 100) + 1}"


args = sys.argv[1:]
ptmax_arr = 1e9
if "--ptmax" in args:
    i = args.index("--ptmax")
    ptmax_arr = float(args[i + 1])
    args = args[:i] + args[i + 2:]
S = ROOT.TChain("S")
for a in args:
    S.Add(a)
edges = [0.0, 0.9, 1.5, 2.1, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5]
trk = collections.defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])   # (pt, eta, ibin) -> de, mod, once, p
arr = collections.defaultdict(lambda: [[], 0.0])                # name -> ([true per arrival], model sum)
cur, acc = None, None
for e in S:
    i = next((i for i in range(len(edges) - 1) if edges[i] <= abs(e.eta) < edges[i + 1]), None)
    if i is None or e.lay0 < 0 or e.lay1 < 0:
        continue
    key = (round(e.pt), e.eta, i)
    stack = int(e.lay0) != int(e.lay1) and phys(int(e.lay0)) == phys(int(e.lay1))
    t = trk[key]
    t[0] += e.de
    t[1] += e.de_mod
    t[2] += 0.0 if stack else e.de_mod
    t[3] = e.p1
    nm = name(phys(int(e.lay1)))
    if e.pt > ptmax_arr:
        continue
    if acc is not None and (acc[0] != key or acc[1] != nm):
        arr[acc[1]][0].append(acc[2])
        arr[acc[1]][1] += acc[3]
        acc = None
    if acc is None:
        acc = [key, nm, 0.0, 0.0]
    acc[2] += e.de
    acc[3] += e.de_mod
if acc is not None:
    arr[acc[1]][0].append(acc[2])
    arr[acc[1]][1] += acc[3]

pts = sorted({k[0] for k in trk})
print("per track between the first and the last hit: true energy loss (MeV) with its statistical error, mkFit's / "
      "true, the same without the 5 % of tracks with the largest true loss, mkFit with the stack value once / true, "
      "true loss / p (%)")
print(f"{'pT':>4} {'|eta|':>9} {'n':>5} {'true MeV':>8} {'+-':>5} {'mkFit/tr':>8} {'trunc95':>7} {'once/tr':>8} "
      f"{'dE/p %':>7}")
for pt in pts:
    for i in range(len(edges) - 1):
        v = [x for k, x in trk.items() if k[0] == pt and k[2] == i]
        if not v:
            continue
        de, mo, on = sum(x[0] for x in v), sum(x[1] for x in v), sum(x[2] for x in v)
        rel = sum(x[0] / x[3] for x in v if x[3] > 0) / len(v)
        mean = de / len(v)
        err = math.sqrt(sum((x[0] - mean) ** 2 for x in v) / len(v)) / math.sqrt(len(v))
        vt = sorted(v, key=lambda x: x[0])[:max(1, int(0.95 * len(v)))]
        tr = sum(x[1] for x in vt) / sum(x[0] for x in vt)
        print(f"{pt:4d} {edges[i]:.1f}-{edges[i + 1]:.1f} {len(v):5d} {1000 * mean:8.2f} {1000 * err:5.2f} {mo / de:8.2f} "
              f"{tr:7.2f} {on / de:8.2f} {100 * rel:7.3f}")


def order(n):
    a = n.split()
    return (["IT", "OT", "TFPX", "TEPX", "TEDD"].index(a[0]), int(a[1][1:]))


print(f"\nper arrival at a physical layer (pT <= {ptmax_arr:g}): true loss mean and median (keV), mkFit's mean (keV), "
      "mkFit / true")
print(f"{'layer':>9} {'n':>6} {'true mean':>9} {'median':>7} {'mkFit':>7} {'ratio':>6}")
for nm in sorted(arr, key=order):
    v, msum = arr[nm]
    if len(v) < 20:
        continue
    vs = sorted(v)
    m = sum(v) / len(v)
    print(f"{nm:>9} {len(v):6d} {1e6 * m:9.1f} {1e6 * vs[len(vs) // 2]:7.1f} {1e6 * msum / len(v):7.1f} "
          f"{msum / len(v) / m:6.2f}")
