#!/usr/bin/env python3
"""g4s_tkmat_once.py OUT...: g4s_tkmat.py's per-track material per |eta| bin, with mkFit's per-module material as
applied (each sensor of a PS/2S stack carries the stack's MediumProperties) and as it would be with the stack's
value applied once: the step from one sensor of a stack to the other (two mkFit layers of the same physical layer)
gets no model material."""
import sys, collections, ROOT


def phys(lay):
    """the physical layer of an mkFit layer: stack partners (4/5, ..., 14/15, 28/29, ..., 36/37 and the -z ones) map
    to one number, the pixel layers each to their own"""
    if lay <= 3 or 16 <= lay <= 27 or 38 <= lay <= 49:
        return lay
    if lay <= 15:
        return 100 + (lay - 4) // 2
    return 200 + (lay - 28) // 2 if lay <= 37 else 300 + (lay - 50) // 2


S = ROOT.TChain("S")
for a in sys.argv[1:]:
    S.Add(a)
edges = [0.0, 0.5, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5]
acc = collections.defaultdict(lambda: [0.0, 0.0, 0.0, set()])
for e in S:
    i = next((i for i in range(len(edges) - 1) if edges[i] <= abs(e.eta) < edges[i + 1]), None)
    if i is None or e.lay0 < 0 or e.lay1 < 0:
        continue
    stack = int(e.lay0) != int(e.lay1) and phys(int(e.lay0)) == phys(int(e.lay1))
    a = acc[i]
    a[0] += e.x0
    a[1] += e.rl_mod
    a[2] += 0.0 if stack else e.rl_mod
    a[3].add((e.eta, e.pt))
print("X0 per track between the first and the last hit: true, mkFit per sensor (as applied), mkFit once per stack")
print(f"{'|eta|':>9} {'true':>7} {'mod':>7} {'once':>7} {'mod/tr':>7} {'once/tr':>7}")
for i in sorted(acc):
    t, m, o, n = acc[i][0], acc[i][1], acc[i][2], len(acc[i][3])
    print(f"{edges[i]:.1f}-{edges[i + 1]:.1f} {t / n:7.3f} {m / n:7.3f} {o / n:7.3f} {m / t:7.2f} {o / t:7.2f}")
