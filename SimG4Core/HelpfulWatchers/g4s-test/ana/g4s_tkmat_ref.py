#!/usr/bin/env python3
"""g4s_tkmat_ref.py OUT.csv OUT.md FILE...: a reference table of the true tracker material, from g4s_tkmat.py's step
trees: per arrival at a physical layer (the steps whose destination hit is in that layer, summed per track: the
material since the last hit in an earlier layer, the stack partner's sensor and overlaps included) and per |eta|
bin of the track: the number of arrivals, the mean and rms of the true X0, the mean true lambda_I, and mkFit's
applied radl * p/|p.n| as applied (per sensor) and with the PS/2S stack value once. An arrival after a layer the
track crossed without a hit also carries that layer's material (6-10 % of the TEPX arrivals in T2's forward muons).
Both sides of the detector are folded. The CSV has one row per (layer, |eta| bin); the Markdown file has the matrix
of the mean true X0 and of mkFit / true."""
import sys, math, collections, ROOT


def phys(lay):
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


def order(n):
    a = n.split()
    return (["IT", "OT", "TFPX", "TEPX", "TEDD"].index(a[0]), int(a[1][1:]))


EDGES = [0.0, 0.5, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5]
S = ROOT.TChain("S")
for a in sys.argv[3:]:
    S.Add(a)
arr = collections.defaultdict(list)  # (name, ibin) -> [(x0, lam, mod, once)]
acc = None
for e in S:
    if e.lay0 < 0 or e.lay1 < 0:
        continue
    ib = next((i for i in range(len(EDGES) - 1) if EDGES[i] <= abs(e.eta) < EDGES[i + 1]), None)
    if ib is None:
        continue
    key = (e.eta, e.pt)
    nm = name(phys(int(e.lay1)))
    stack = int(e.lay0) != int(e.lay1) and phys(int(e.lay0)) == phys(int(e.lay1))
    if acc is not None and (acc[0] != key or acc[1] != nm):
        arr[(acc[1], acc[2])].append(acc[3])
        acc = None
    if acc is None:
        acc = [key, nm, ib, [0.0, 0.0, 0.0, 0.0]]
    v = acc[3]
    v[0] += e.x0
    v[1] += e.lam
    v[2] += e.rl_mod
    v[3] += 0.0 if stack else e.rl_mod
if acc is not None:
    arr[(acc[1], acc[2])].append(acc[3])

layers = sorted({k[0] for k in arr}, key=order)
with open(sys.argv[1], "w") as f:
    f.write("layer,eta_lo,eta_hi,n,x0_mean,x0_rms,lambda_mean,mkfit_per_sensor,mkfit_once_per_stack\n")
    for nm in layers:
        for ib in range(len(EDGES) - 1):
            v = arr.get((nm, ib))
            if not v:
                continue
            n = len(v)
            m = [sum(x[j] for x in v) / n for j in range(4)]
            rms = math.sqrt(sum((x[0] - m[0]) ** 2 for x in v) / n)
            f.write(f"{nm},{EDGES[ib]},{EDGES[ib + 1]},{n},{m[0]:.5f},{rms:.5f},{m[1]:.5f},{m[2]:.5f},{m[3]:.5f}\n")


def cell(nm, ib, what):
    v = arr.get((nm, ib))
    if not v or len(v) < 20:
        return ""
    n = len(v)
    t = sum(x[0] for x in v) / n
    if what == "x0":
        return f"{t:.4f}"
    return f"{sum(x[2] for x in v) / n / t:.2f}"


hdr = "| layer | " + " | ".join(f"{EDGES[i]:g}-{EDGES[i + 1]:g}" for i in range(len(EDGES) - 1)) + " |\n"
sep = "|---|" + "---|" * (len(EDGES) - 1) + "\n"
with open(sys.argv[2], "w") as f:
    f.write("Mean true X0 per arrival at a physical layer, by the track's |eta| (cells with at least 20 arrivals)\n\n")
    f.write(hdr + sep)
    for nm in layers:
        f.write(f"| {nm} | " + " | ".join(cell(nm, i, "x0") for i in range(len(EDGES) - 1)) + " |\n")
    f.write("\nmkFit's applied material (per sensor, as in the code) over the true X0, same cells\n\n")
    f.write(hdr + sep)
    for nm in layers:
        f.write(f"| {nm} | " + " | ".join(cell(nm, i, "ratio") for i in range(len(EDGES) - 1)) + " |\n")
print(f"{sum(len(v) for v in arr.values())} arrivals, {len(layers)} layers")
