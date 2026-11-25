#!/usr/bin/env python3
"""g4s_material.py FILE: for the primary of each event, the material crossed from the origin to |z| = 300 cm (ETL
front), 320 cm (HGCal front face) and 520 cm (HGCal back), in X0 and lambda_I, interpolated linearly in |z| between
path points, in bins of |eta|; and whether the sums along the path never decrease."""
import sys, math, ROOT, collections
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
sf = SnitchFile(sys.argv[1]); f = sf.f
zs = (300.0, 320.0, 520.0)
acc = collections.defaultdict(lambda: [0, [0.0] * 3, [0.0] * 3])
bad = 0
for ev in sf:
    p, s = ev.p, ev.s
    k = p[0].m_daughters_begin
    q = p[k]; path = s[k].m_path
    eta = abs(math.asinh(q.m_p_birth.fZ / math.hypot(q.m_p_birth.fX, q.m_p_birth.fY)))
    b = round(math.floor(eta * 5) / 5, 1)
    a = acc[b]; a[0] += 1
    prev = None
    for j in range(path.size()):
        if prev is not None and (path[j].m_x0 < prev[0] - 1e-6 or path[j].m_lambda < prev[1] - 1e-6): bad += 1
        prev = (path[j].m_x0, path[j].m_lambda)
    for iz, z in enumerate(zs):
        for j in range(1, path.size()):
            z0, z1 = abs(path[j - 1].m_x.fZ), abs(path[j].m_x.fZ)
            if z0 <= z < z1:
                w = (z - z0) / (z1 - z0)
                a[1][iz] += path[j - 1].m_x0 + w * (path[j].m_x0 - path[j - 1].m_x0)
                a[2][iz] += path[j - 1].m_lambda + w * (path[j].m_lambda - path[j - 1].m_lambda)
                break
print(f"{sys.argv[1]}: path points where a sum decreases: {bad}")
print("|eta|    n   X0 to z=300  320   520    lambda to z=300  320   520")
for b in sorted(acc):
    n, x, l = acc[b]
    print(f"{b:.1f}-{b + 0.2:.1f} {n:4d}   " + "  ".join(f"{v / n:5.2f}" for v in x) + "        " + "  ".join(f"{v / n:5.2f}" for v in l))
