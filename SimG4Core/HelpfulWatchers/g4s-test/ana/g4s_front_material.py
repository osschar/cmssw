#!/usr/bin/env python3
"""g4s_front_material.py FILE [ZMIN ZMAX]: the primary's material between |z| ZMIN and ZMAX (default 280-325 cm), per
volume, in X0 per particle: each own stored step contributes the difference of the sums at its start and at the next
own stored step's start (so unstored steps between them are included), filed under the step's volume. For muons."""
import sys, collections, ROOT
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
ROOT.gSystem.Load("libRootG4Snitch.so")
sf = SnitchFile(sys.argv[1])
zmin, zmax = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (280.0, 325.0)
vols = [str(v.path) for v in sf.rows("V")]
acc = collections.Counter(); n = 0
for ev in sf:
    p, s = ev.p, ev.s
    k = p[0].m_daughters_begin; q = p[k]; n += 1
    steps = sorted([st for st in list(s[k].m_steps_sensitive) + list(s[k].m_steps_inert) if st.m_g4_pid == q.m_g4_id],
                   key=lambda st: st.m_x0_beg)
    for a, b in zip(steps, steps[1:]):
        if zmin <= abs(a.m_x_beg.fZ) < zmax:
            name = "/".join(x.split(":")[0] for x in vols[a.m_volume].split("/")[2:5]) + " .. " + vols[a.m_volume].split("/")[-1].split(":")[0]
            acc[name] += b.m_x0_beg - a.m_x0_beg
tot = sum(acc.values())
print(f"{n} particles, X0 between |z| {zmin:.0f} and {zmax:.0f} cm: {tot / n:.3f} per particle")
for name, x in acc.most_common(16): print(f"  {x / n:6.3f}  {name}")
