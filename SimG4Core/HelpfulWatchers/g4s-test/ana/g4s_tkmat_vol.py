#!/usr/bin/env python3
"""g4s_tkmat_vol.py ETAMIN ETAMAX FILE...: where the primary's X0 between its first and its last tracker hit lies, per
Geant4 volume, for tracks with ETAMIN <= |eta| < ETAMAX (single gun, one primary per event; the outgoing leg in the
tracker as in g4s_tkmat.py). Each own stored step contributes its length over its volume's X0 (table V); the
coverage line compares their sum with the snitch's running X0 sum over the same range (steps not stored are missing
from the first). Volumes inside a tracker module (a path component containing "Module") are grouped per
sub-detector as "<sub-detector> modules: sensors" (sensitive) and "...: other"; the rest by sub-detector and leaf
name (copy numbers dropped)."""
import sys, math, collections, ROOT
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
ROOT.gSystem.Load("libRootG4Snitch.so")

emin, emax = float(sys.argv[1]), float(sys.argv[2])
acc = collections.Counter()
ntrk, x0_run, x0_steps = 0, 0.0, 0.0
for fn in sys.argv[3:]:
    sf = SnitchFile(fn); f = sf.f
    sdcls = [str(getattr(d, "class")) for d in sf.rows("D")]
    trk_sd = {i for i, c in enumerate(sdcls) if "TkAccumulating" in c}
    V = []
    for v in sf.rows("V"):
        comps = [c.split(":")[0] for c in str(v.path).split("/")]
        if "Tracker" in comps:
            sub = comps[comps.index("Tracker") + 1:]
            det = sub[0] if sub else "Tracker"
            if any("Module" in c for c in sub):
                key = f"{det} modules: " + ("sensors" if int(v.sd) in trk_sd else "other")
            else:
                key = f"{det}/{sub[-1]}" if len(sub) > 1 else det
        else:
            key = "/".join(comps[2:4]) + "/" + comps[-1]
        V.append((key, float(v.x0), int(v.sd) in trk_sd))
    for ev in sf:
        p, s = ev.p, ev.s
        k = p[0].m_daughters_begin
        q = p[k]
        pb = q.m_p_birth
        pt = math.hypot(pb.fX, pb.fY)
        if q.m_charge == 0 or not emin <= abs(math.asinh(pb.fZ / pt)) < emax:
            continue
        gid = q.m_g4_id
        tcut, rmax = 1e30, 0.0
        for pp in s[k].m_path:
            r = math.hypot(pp.m_x.fX, pp.m_x.fY)
            rmax = max(rmax, r)
            if r > 118.0 or abs(pp.m_x.fZ) > 282.0 or r < rmax - 1.0:
                tcut = pp.m_x.fT
                break
        own = sorted([st for st in list(s[k].m_steps_sensitive) + list(s[k].m_steps_inert)
                      if st.m_g4_pid == gid and st.m_x_beg.fT < tcut], key=lambda st: st.m_x_beg.fT)
        hit_i = [i for i, st in enumerate(own) if st.m_det_id != 0 and 0 <= st.m_volume < len(V) and V[st.m_volume][2]]
        if len(hit_i) < 2:
            continue
        ntrk += 1
        a, b = hit_i[0], hit_i[-1]
        x0_run += own[b].m_x0_beg - own[a].m_x0_beg
        for st in own[a:b]:
            L = math.sqrt((st.m_x_end.fX - st.m_x_beg.fX) ** 2 + (st.m_x_end.fY - st.m_x_beg.fY) ** 2 +
                          (st.m_x_end.fZ - st.m_x_beg.fZ) ** 2)
            key, x0, _ = V[st.m_volume]
            if x0 > 0:
                acc[key] += L / x0
                x0_steps += L / x0
    f.Close()
print(f"{ntrk} tracks at |eta| {emin}-{emax}; X0 per track between the first and the last hit: running sum "
      f"{x0_run / ntrk:.4f}, stored steps {x0_steps / ntrk:.4f} (coverage {x0_steps / x0_run:.3f})")
for key, x in acc.most_common(25):
    print(f"  {x / ntrk:7.4f}  {100 * x / x0_steps:5.1f} %  {key}")
