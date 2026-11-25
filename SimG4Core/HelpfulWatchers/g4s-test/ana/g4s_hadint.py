#!/usr/bin/env python3
"""g4s_hadint.py OUT.root FILE...: nuclear interactions, decays and kinks of single charged hadrons in the tracker,
for mkFit's stopping and hole policy (single gun; the primary's hits and outgoing leg as in g4s_tkmat.py).

Per primary (tree "H"):
  eta, pt, pdg, nhits (own tracker hits, distinct DetIds in order), fate (0 left the tracker or stopped after its
  last hit outside, 1 inelastic interaction inside r < 118 cm, |z| < 282 cm, 2 decay inside, 3 other end inside),
  death_r, death_z, p_before (|p| before the last step);
  for fate 1 and 2, the leading charged daughter (largest |p|, born in the primary's last step, any charged
  particle kept by the snitch, i.e. above 10 MeV kinetic, nuclei and nuclear fragments excluded): lead_pdg, lead_x = |p_lead| / p_before, lead_ang (rad,
  to the primary's direction in its last step), lead_hits (its own tracker hits), lead_pt; n_ch03 (charged
  daughters above 0.3 GeV);
  kinks along the primary between its first and its last hit, from steps ended by hadElastic or CoulombScat: the
  angle between the step's own direction and its end momentum: nk5, nk20, nk100 (above 5, 20, 100 mrad), kmax
  (largest), kmax_r (its radius), kmax_dp (relative momentum change)."""
import sys, math, array, ROOT
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
ROOT.gSystem.Load("libRootG4Snitch.so")


def unit(v):
    n = math.sqrt(sum(x * x for x in v))
    return [x / n for x in v] if n > 0 else None


def angle(a, b):
    ua, ub = unit(a), unit(b)
    if ua is None or ub is None:
        return -1.0
    return math.acos(max(-1.0, min(1.0, sum(x * y for x, y in zip(ua, ub)))))


def main():
    out = ROOT.TFile(sys.argv[1], "RECREATE")
    cols = ["eta", "pt", "pdg", "nhits", "fate", "death_r", "death_z", "p_before", "lead_pdg", "lead_x", "lead_ang",
            "lead_hits", "lead_pt", "n_ch03", "nk5", "nk20", "nk100", "kmax", "kmax_r", "kmax_dp"]
    H = ROOT.TTree("H", "hadrons in the tracker")
    b = {c: array.array("f", [0.0]) for c in cols}
    for c in cols:
        H.Branch(c, b[c], c + "/F")
    for fn in sys.argv[2:]:
        sf = SnitchFile(fn); f = sf.f
        sdcls = [str(getattr(d, "class")) for d in sf.rows("D")]
        trk_sd = {i for i, c in enumerate(sdcls) if "TkAccumulating" in c}
        V = [int(v.sd) in trk_sd for v in sf.rows("V")]
        P = [str(x.name) for x in sf.rows("P")]
        for ev in sf:
            p, s = ev.p, ev.s
            k = p[0].m_daughters_begin
            q = p[k]
            if q.m_charge == 0:
                continue

            def own_steps(j):
                g = p[j].m_g4_id
                return sorted([st for st in list(s[j].m_steps_sensitive) + list(s[j].m_steps_inert)
                               if st.m_g4_pid == g], key=lambda st: st.m_x_beg.fT)

            def tracker_hits(steps, tcut=1e30):
                ids = []
                for st in steps:
                    if st.m_x_beg.fT >= tcut:
                        break
                    if st.m_det_id != 0 and 0 <= st.m_volume < len(V) and V[st.m_volume]:
                        if not ids or ids[-1] != st.m_det_id:
                            ids.append(st.m_det_id)
                return ids

            tcut, rmax = 1e30, 0.0
            for pp in s[k].m_path:
                r = math.hypot(pp.m_x.fX, pp.m_x.fY)
                rmax = max(rmax, r)
                if r > 118.0 or abs(pp.m_x.fZ) > 282.0 or r < rmax - 1.0:
                    tcut = pp.m_x.fT
                    break
            own = own_steps(k)
            ids = tracker_hits(own, tcut)
            pb = q.m_p_birth
            pt = math.hypot(pb.fX, pb.fY)
            vals = {c: 0.0 for c in cols}
            vals.update(eta=math.asinh(pb.fZ / pt), pt=pt, pdg=q.m_pdg, nhits=len(ids), lead_pdg=0, lead_x=-1,
                        lead_ang=-1, lead_hits=-1, lead_pt=-1, kmax=0, kmax_r=-1, kmax_dp=0)
            xd = q.m_x_death
            dr, dz = math.hypot(xd.fX, xd.fY), xd.fZ
            vals["death_r"], vals["death_z"] = dr, dz
            dproc = P[q.m_death_proc] if 0 <= q.m_death_proc < len(P) else ""
            inside = dr < 118.0 and abs(dz) < 282.0 and xd.fT < tcut
            if inside and "Inelastic" in dproc:
                vals["fate"] = 1
            elif inside and dproc == "Decay":
                vals["fate"] = 2
            elif inside:
                vals["fate"] = 3
            # the primary's direction and momentum before its last step
            if own:
                last = own[-1]
                ddir = [last.m_x_end.fX - last.m_x_beg.fX, last.m_x_end.fY - last.m_x_beg.fY,
                        last.m_x_end.fZ - last.m_x_beg.fZ]
                prev = own[-2] if len(own) > 1 else None
                pbef = math.sqrt(prev.m_p_end.fX ** 2 + prev.m_p_end.fY ** 2 + prev.m_p_end.fZ ** 2) if prev else 0.0
                vals["p_before"] = pbef
                if vals["fate"] in (1, 2) and pbef > 0:
                    lead, lp, nch = None, 0.0, 0
                    for j in range(q.m_daughters_begin, q.m_daughters_end):
                        d = p[j]
                        if d.m_charge == 0 or abs(d.m_pdg) >= 1000000000:
                            continue
                        db = d.m_x_birth
                        if abs(db.fT - xd.fT) > 1e-3:
                            continue
                        pd = d.m_p_birth
                        pm = math.sqrt(pd.fX ** 2 + pd.fY ** 2 + pd.fZ ** 2)
                        nch += pm > 0.3
                        if pm > lp:
                            lead, lp = j, pm
                    vals["n_ch03"] = nch
                    if lead is not None:
                        d = p[lead]
                        pd = d.m_p_birth
                        vals.update(lead_pdg=d.m_pdg, lead_x=lp / pbef, lead_ang=angle(ddir, [pd.fX, pd.fY, pd.fZ]),
                                    lead_hits=len(tracker_hits(own_steps(lead))), lead_pt=math.hypot(pd.fX, pd.fY))
            # kinks along the primary between its first and its last hit
            hit_t = [st.m_x_beg.fT for st in own if st.m_det_id != 0 and 0 <= st.m_volume < len(V) and V[st.m_volume]
                     and st.m_x_beg.fT < tcut]
            if len(hit_t) >= 2:
                t0, t1 = hit_t[0], hit_t[-1]
                prevp = None
                for st in own:
                    if t0 <= st.m_x_beg.fT <= t1 and 0 <= st.m_proc < len(P) and P[st.m_proc] in ("hadElastic",
                                                                                                    "CoulombScat"):
                        d = [st.m_x_end.fX - st.m_x_beg.fX, st.m_x_end.fY - st.m_x_beg.fY, st.m_x_end.fZ - st.m_x_beg.fZ]
                        pe = [st.m_p_end.fX, st.m_p_end.fY, st.m_p_end.fZ]
                        a = angle(d, pe)
                        if a < 0:
                            continue
                        vals["nk5"] += a > 0.005
                        vals["nk20"] += a > 0.020
                        vals["nk100"] += a > 0.100
                        if a > vals["kmax"]:
                            pa = math.sqrt(sum(x * x for x in pe))
                            vals.update(kmax=a, kmax_r=math.hypot(st.m_x_end.fX, st.m_x_end.fY),
                                        kmax_dp=(pa - prevp) / prevp if prevp else 0.0)
                    prevp = math.sqrt(st.m_p_end.fX ** 2 + st.m_p_end.fY ** 2 + st.m_p_end.fZ ** 2)
            for c in cols:
                b[c][0] = vals[c]
            H.Fill()
        f.Close()
    print(f"{H.GetEntries()} primaries")
    out.cd()
    H.Write()
    out.Close()


main()
