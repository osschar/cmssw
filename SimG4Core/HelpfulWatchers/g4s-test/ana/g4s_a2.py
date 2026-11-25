#!/usr/bin/env python3
"""g4s_a2.py FILE...: the MIP segment of single particles from the origin (G4SNITCH_GUN=single samples).

For each event, the primary's leading chain: the primary, continued by a daughter of the same kind born where it
died with more than half of its last energy (Geant4's eBrem can end an electron and start a new one). The shower
starts, for a hadron, at its first inelastic interaction (the chain dies by a process named *Inelastic); for a
photon at its conversion (conv); for an electron where the chain's energy fell below half of the initial (and, as
the first sign, below 90 %). Reported: where (|z| < 300 cm: tracker; 300-320: ETL and HGCal front; then CE-E and CE-H,
split at the largest |z| of an EE cell's step; or nowhere), the material from the origin in X0 and lambda_I at the
start (medians), and the number of distinct HGCal layers the chain crossed with its own sensitive steps before it."""
import sys, math, statistics, collections, ROOT
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
ROOT.gSystem.Load("libRootG4Snitch.so")

def hgc_layer(d):
    det = d >> 28
    if det == 8: return (0, (d >> 20) & 0x1F)    # CE-E silicon
    if det == 9: return (1, (d >> 20) & 0x1F)    # CE-H silicon
    if det == 10: return (1, (d >> 17) & 0x1F)   # CE-H scintillator
    return None

def analyse(fn):
    sf = SnitchFile(fn)
    procs = {int(p.idx): str(p.name) for p in sf.rows("P")}
    # The CE-E / CE-H boundary from the geometry: halfway between the last CE-E and the first CE-H sensitive volume.
    sdname = {int(d.idx): str(d.name) for d in sf.rows("D")}
    ee_z = [abs(v.z) for v in sf.rows("V") if sdname.get(int(v.sd)) == "HGCHitsEE"]
    he_z = [abs(v.z) for v in sf.rows("V") if sdname.get(int(v.sd)) in ("HGCHitsHEfront", "HGCHitsHEback")]
    ee_zmax = 0.5 * (max(ee_z) + min(he_z)) if ee_z and he_z else 362.0
    rows = []
    for ev in sf:
        p, s = ev.p, ev.s
        k = p[0].m_daughters_begin
        q0 = p[k]; E0 = q0.m_p_birth.fT; kind = abs(q0.m_pdg)
        pts = []; layers_seq = []; off_x0 = off_l = 0.0; cur = k; death = None
        while True:
            q = p[cur]; ps = s[cur]
            for pp in ps.m_path:
                pts.append((abs(pp.m_x.fZ), pp.m_p.fT, off_x0 + pp.m_x0, off_l + pp.m_lambda))
            for st in ps.m_steps_sensitive:
                if st.m_g4_pid != q.m_g4_id or st.m_det_id == 0: continue
                hl = hgc_layer(st.m_det_id)
                if hl:
                    layers_seq.append((abs(st.m_x_beg.fZ), hl))
            if ps.m_path.size():
                off_x0 += ps.m_path[ps.m_path.size() - 1].m_x0; off_l += ps.m_path[ps.m_path.size() - 1].m_lambda
            death = (q, procs.get(q.m_death_proc, "?"))
            if kind not in (11, 13, 22) and "Inelastic" in death[1]:
                break  # a hadron's shower starts at its first inelastic interaction, even if a leading hadron goes on
            e_last = ps.m_path[ps.m_path.size() - 1].m_p.fT if ps.m_path.size() else q.m_p_birth.fT
            nxt = None; best = 0.0
            for j in range(q.m_daughters_begin, q.m_daughters_end):
                d = p[j]
                if abs(d.m_pdg) != kind: continue
                dx = math.dist((d.m_x_birth.fX, d.m_x_birth.fY, d.m_x_birth.fZ), (q.m_x_death.fX, q.m_x_death.fY, q.m_x_death.fZ))
                if dx < 1e-3 and d.m_p_birth.fT > 0.5 * e_last and d.m_p_birth.fT > best: nxt, best = j, d.m_p_birth.fT
            if nxt is None: break
            cur = nxt
        start = None; first = None
        if kind == 11:
            for z, e, x0, l in pts:
                if first is None and e < 0.9 * E0: first = (z, x0, l)
                if e < 0.5 * E0: start = (z, x0, l); break
        else:
            dq, dproc = death
            if (kind == 22 and dproc == "conv") or (kind not in (11, 13, 22) and "Inelastic" in dproc):
                start = (abs(dq.m_x_death.fZ), pts[-1][2] if pts else 0.0, pts[-1][3] if pts else 0.0)
        nlay = len({hl for z, hl in layers_seq if start is None or z < start[0]})
        rows.append((start, first, nlay, death[1]))
    return rows, ee_zmax

def kind_of(name): return name.split("_")[0]
for fn in sys.argv[1:]:
    rows, ee_zmax = analyse(fn)
    n = len(rows); reg = collections.Counter(); x0s = []; ls = []; nl = []; firsts = []; none_end = collections.Counter()
    for start, first, nlay, dproc in rows:
        nl.append(nlay)
        if first: firsts.append(first[1])
        if start is None: reg["none"] += 1; none_end[dproc] += 1; continue
        z = start[0]
        r = "trk" if z < 300 else ("front" if z < 320 else ("CE-E" if z <= ee_zmax + 0.5 else "CE-H"))
        reg[r] += 1; x0s.append(start[1]); ls.append(start[2])
    med = lambda v: f"{statistics.median(v):6.2f}" if v else "     -"
    name = fn.rstrip("/").split("/")[-2]
    print(f"{name:12s} n {n:4d}  start: trk {reg['trk']:4d} front {reg['front']:4d} CE-E {reg['CE-E']:4d} CE-H {reg['CE-H']:4d} none {reg['none']:4d}"
          f"  | median X0 {med(x0s)} lambda {med(ls)}" + (f" (90% E: X0 {med(firsts)})" if firsts else "")
          + f"  | MIP layers median {statistics.median(nl):4.1f}"
          + (f"  | none end by {dict(none_end.most_common(3))}" if none_end and kind_of(name) != "mu" else ""))
