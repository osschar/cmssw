#!/usr/bin/env python3
"""g4s_tkmat_sum.py FILE... [-- PDG]: tables of g4s_tkmat.py's output.
1. Per |eta| bin, per track: truth X0 and lambda_I between the first and the last hit, mkFit's applied material (per
   module, grid), their ratio, the truth X0 from the origin to the first hit, and the fraction of hadrons that
   interact inelastically inside the tracker.
2. Per arrival at a physical layer (the steps whose destination hit is in that layer, summed per track: the stereo
   pair's second sensor and overlaps included, so the sum is the material since the last hit in an earlier layer):
   truth X0 against the model, mean over arrivals, and the spread (rms / mean) of truth and of the ratio per arrival.
   TEDD split at r = 60 cm (inner rings PS, outer 2S, approximately)."""
import sys, math, collections, ROOT

args = sys.argv[1:]
pdg = None
if "--" in args:
    i = args.index("--")
    pdg = int(args[i + 1])
    args = args[:i]


def phys(lay, r):
    if lay < 0:
        return "?"
    if lay <= 3:
        return f"IT L{lay + 1}"
    if lay <= 15:
        return f"OT L{(lay - 4) // 2 + 1}"
    l = lay if lay < 38 else lay - 22
    if l <= 23:
        return f"TFPX D{l - 15}"
    if l <= 27:
        return f"TEPX D{l - 23}"
    return f"TEDD D{(l - 28) // 2 + 1} " + ("r<60" if r < 60 else "r>60")


def order(name):
    keys = ["IT", "OT", "TFPX", "TEPX", "TEDD"]
    a = name.split()
    return (keys.index(a[0]) if a[0] in keys else 9, int(a[1][1:]) if len(a) > 1 else 0, a[2] if len(a) > 2 else "")


K = ROOT.TChain("K")
S = ROOT.TChain("S")
for a in args:
    K.Add(a)
    S.Add(a)
edges = [0.0, 0.5, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5]
rows = collections.defaultdict(list)
for e in K:
    if pdg is not None and int(e.pdg) != pdg:
        continue
    for i in range(len(edges) - 1):
        if edges[i] <= abs(e.eta) < edges[i + 1]:
            rows[i].append((e.nhits, e.x0, e.lam, e.rl_mod, e.rl_grid, e.x0_first, e.inel))


def mean(v):
    return sum(v) / len(v) if v else float("nan")


def rms_rel(v):
    m = mean(v)
    return math.sqrt(mean([(x - m) ** 2 for x in v])) / m if v and m else float("nan")


print(f"{K.GetEntries()} tracks" + (f", pdg {pdg}" if pdg is not None else ""))
print("per track, between the first and the last hit (X0, lambda_I); model = sum of destination radl * p/|p.n|")
print(f"{'|eta|':>9} {'n':>5} {'nh':>5} {'X0 true':>8} {'mod':>7} {'grid':>7} {'mod/tr':>7} {'grid/tr':>7} "
      f"{'lam true':>8} {'X0 0->1':>8} {'inel':>6}")
for i in sorted(rows):
    r = rows[i]
    xt, xm, xg = mean([x[1] for x in r]), mean([x[3] for x in r]), mean([x[4] for x in r])
    print(f"{edges[i]:.1f}-{edges[i + 1]:.1f} {len(r):5d} {mean([x[0] for x in r]):5.1f} {xt:8.3f} {xm:7.3f} {xg:7.3f} "
          f"{xm / xt:7.2f} {xg / xt:7.2f} {mean([x[2] for x in r]):8.3f} {mean([x[5] for x in r]):8.3f} "
          f"{mean([x[6] for x in r]):6.3f}")

# per arrival at a physical layer
arr = collections.defaultdict(list)  # name -> [(truth, mod, grid)]
cur = None
acc = None
key = None
for e in S:
    if pdg is not None and int(e.pdg) != pdg:
        continue
    tk = (e.eta, e.pt)
    name = phys(int(e.lay1), e.r1)
    if key is not None and (tk != key[0] or name != key[1]):
        arr[key[1]].append(acc)
        acc = None
    if acc is None:
        acc = [0.0, 0.0, 0.0]
    acc[0] += e.x0
    acc[1] += e.rl_mod
    acc[2] += e.rl_grid
    key = (tk, name)
if key is not None:
    arr[key[1]].append(acc)
print("\nper arrival at a physical layer: material since the last hit in an earlier layer (X0)")
print(f"{'layer':>14} {'n':>6} {'true':>7} {'mod':>7} {'grid':>7} {'mod/tr':>7} {'grid/tr':>7} {'rms/mean tr':>11} "
      f"{'rms/mean ratio':>14}")
for name in sorted(arr, key=order):
    v = arr[name]
    if len(v) < 20:
        continue
    t, m, g = mean([x[0] for x in v]), mean([x[1] for x in v]), mean([x[2] for x in v])
    ratios = [x[1] / x[0] for x in v if x[0] > 0]
    print(f"{name:>14} {len(v):6d} {t:7.4f} {m:7.4f} {g:7.4f} {m / t:7.2f} {g / t:7.2f} {rms_rel([x[0] for x in v]):11.2f} "
          f"{rms_rel(ratios):14.2f}")
