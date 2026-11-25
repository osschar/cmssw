#!/usr/bin/env python3
"""g4s_primaries.py FILE: per event, the primaries (PDG, pT, eta, phi) and the energy deposited by the sensitive steps
and the inert sums at z > 0 and at z < 0 (GeV), with the event's particle count; reads either format (g4s_read)."""
import sys, math
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile

sf = SnitchFile(sys.argv[1])
print(f"{sys.argv[1]}: {'TTree' if sf.is_ttree else 'RNTuple'}, {len(sf)} events")
for ev in sf:
    p, s = ev.p, ev.s
    prims = []
    for k in range(p[0].m_daughters_begin, p[0].m_daughters_end):
        m = p[k].m_p_birth
        pt = math.hypot(m.fX, m.fY)
        eta = math.asinh(m.fZ / pt) if pt > 0 else float("nan")
        prims.append(f"{p[k].m_pdg:5d} pT {pt:5.1f} eta {eta:+.2f} phi {math.atan2(m.fY, m.fX):+.2f}")
    e = [0.0, 0.0]  # MeV at z > 0, z < 0
    for ps in s:
        for x in ps.m_steps_sensitive:
            e[x.m_x_end.fZ < 0] += x.m_edep.m_total
        for x in ps.m_steps_inert:
            e[x.m_x_end.fZ < 0] += x.m_edep.m_total
        for x in ps.m_inert_sums:
            e[x.m_x_mean.fZ < 0] += x.m_edep.m_total
    ev_no = int(ev.ev) if "ev" in sf.fields else -1
    print(f"ev {ev_no:4d}: {p.size():6d} particles; deposits z>0 {e[0] / 1000:6.1f} GeV, z<0 {e[1] / 1000:6.1f} GeV; "
          + " | ".join(prims))
