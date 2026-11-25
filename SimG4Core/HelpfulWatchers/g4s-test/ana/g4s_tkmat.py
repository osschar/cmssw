#!/usr/bin/env python3
"""g4s_tkmat.py GEO.root OUT.root FILE...: mkFit's tracker material against Geant4's, per step between consecutive
hits of the primary (one charged primary per event, the single gun).

GEO.root is mkFit's geometry dumped by dump_mkfit_geo.C (modules by DetId with their normal and MediumProperties
radl, the (|z|, r) grid). A hit is the primary's own sensitive steps in one tracker DetId (consecutive steps merged),
on the outgoing leg in the tracker (hits after the first path point outside r < 118 cm, |z| < 282 cm, or after the
path's radius has fallen 1 cm below its maximum, are dropped; a hadron's path ends at its inelastic interaction).
Truth: X0 and lambda_I along the path, from the mid point of one hit to the mid point of the next: the snitch's
running sums at the first step of a hit plus half of the hit's own steps (each step's length over its volume's X0
or lambda_I from table V). Model, as mkFit applies material on arriving at a hit (applyMaterialEffects): the
destination's radl times p/|p.n|, n the module normal from mkFit's geometry, p the step direction in the sensor:
  rl_mod  : the module's own radl (MediumProperties; the final fit, refitMaterialPerModule)
  rl_grid : the (|z|, r) grid at the hit, 1 cm x 1 cm (track finding)
Tree "S", one entry per step: eta, pt, pdg, istep, lay0, lay1 (mkFit layers of the two hits, -1 if the DetId is not
in mkFit's geometry), r0, z0, r1, z1, x0, lam (truth), rl_mod, rl_grid, invcos.
Tree "K", one entry per track: eta, pt, pdg, nhits, x0_first (truth from the origin to the first hit's mid point),
x0, lam (truth between the first and the last hit's mid points), rl_mod, rl_grid (sums over hits 2..N), inel (1 if
the primary ended by an inelastic interaction inside r < 118, |z| < 282), death_r, death_z."""
import sys, math, array, collections, ROOT
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
ROOT.gSystem.Load("libRootG4Snitch.so")


def main():
    geo = ROOT.TFile.Open(sys.argv[1])
    mods = {}
    for m in geo.Get("M"):
        mods[int(m.detid)] = (int(m.layer), (m.nx, m.ny, m.nz), m.radl)
    hgrid = geo.Get("radl")
    nzb, nrb = hgrid.GetNbinsX(), hgrid.GetNbinsY()
    grid = [[hgrid.GetBinContent(i + 1, j + 1) for j in range(nrb)] for i in range(nzb)]

    def grid_rl(z, r):
        i, j = int(abs(z)), int(r)
        return grid[i][j] if 0 <= i < nzb and 0 <= j < nrb else 0.0

    out = ROOT.TFile(sys.argv[2], "RECREATE")
    scols = ["eta", "pt", "pdg", "istep", "lay0", "lay1", "r0", "z0", "r1", "z1", "x0", "lam", "rl_mod", "rl_grid",
             "invcos"]
    kcols = ["eta", "pt", "pdg", "nhits", "x0_first", "x0", "lam", "rl_mod", "rl_grid", "inel", "death_r", "death_z",
             "nmiss"]
    S = ROOT.TTree("S", "steps between consecutive hits")
    Kt = ROOT.TTree("K", "tracks")
    sb = {c: array.array("f", [0.0]) for c in scols}
    kb = {c: array.array("f", [0.0]) for c in kcols}
    for c in scols:
        S.Branch(c, sb[c], c + "/F")
    for c in kcols:
        Kt.Branch(c, kb[c], c + "/F")
    nomod = collections.Counter()
    for fn in sys.argv[3:]:
        sf = SnitchFile(fn); f = sf.f
        sdcls = [str(getattr(d, "class")) for d in sf.rows("D")]
        V = [(int(v.sd), float(v.x0), float(getattr(v, "lambda"))) for v in sf.rows("V")]
        procs = [str(pp.name) for pp in sf.rows("P")]
        trk_sd = {i for i, c in enumerate(sdcls) if "TkAccumulating" in c}
        for ev in sf:
            p, s = ev.p, ev.s
            k = p[0].m_daughters_begin
            q = p[k]
            if q.m_charge == 0:
                continue
            gid = q.m_g4_id
            own = [st for st in list(s[k].m_steps_sensitive) + list(s[k].m_steps_inert) if st.m_g4_pid == gid]
            own.sort(key=lambda st: st.m_x_beg.fT)
            path = s[k].m_path
            tcut, rmax = 1e30, 0.0
            for pp in path:
                r = math.hypot(pp.m_x.fX, pp.m_x.fY)
                rmax = max(rmax, r)
                if r > 118.0 or abs(pp.m_x.fZ) > 282.0 or r < rmax - 1.0:
                    tcut = pp.m_x.fT
                    break
            # hits: consecutive own sensitive tracker steps of one DetId
            hits = []  # [detid, x0_at_entry, lam_at_entry, dx0 (own steps), dlam, sum of positions, n, direction]
            prev_tracker = False
            for st in own:
                vol = st.m_volume
                if st.m_x_beg.fT >= tcut:
                    break
                is_hit = st.m_det_id != 0 and 0 <= vol < len(V) and V[vol][0] in trk_sd
                if not is_hit:
                    prev_tracker = False
                    continue
                d = (st.m_x_end.fX - st.m_x_beg.fX, st.m_x_end.fY - st.m_x_beg.fY, st.m_x_end.fZ - st.m_x_beg.fZ)
                L = math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
                dx0 = L / V[vol][1] if V[vol][1] > 0 else 0.0
                dlam = L / V[vol][2] if V[vol][2] > 0 else 0.0
                mid = (0.5 * (st.m_x_beg.fX + st.m_x_end.fX), 0.5 * (st.m_x_beg.fY + st.m_x_end.fY),
                       0.5 * (st.m_x_beg.fZ + st.m_x_end.fZ))
                if prev_tracker and hits[-1][0] == st.m_det_id:
                    h = hits[-1]
                    h[3] += dx0
                    h[4] += dlam
                    h[5] = [h[5][c] + mid[c] for c in range(3)]
                    h[6] += 1
                    h[7] = [h[7][c] + d[c] for c in range(3)]
                else:
                    hits.append([st.m_det_id, st.m_x0_beg, st.m_lambda_beg, dx0, dlam, list(mid), 1, list(d)])
                prev_tracker = True
            if len(hits) < 2:
                continue
            pb = q.m_p_birth
            pt = math.hypot(pb.fX, pb.fY)
            eta = math.asinh(pb.fZ / pt)
            xd = q.m_x_death
            dr, dz = math.hypot(xd.fX, xd.fY), xd.fZ
            dproc = procs[q.m_death_proc] if 0 <= q.m_death_proc < len(procs) else ""
            inel = 1.0 if ("Inelastic" in dproc and dr < 118.0 and abs(dz) < 282.0) else 0.0
            rows = []
            for h in hits:
                pos = [x / h[6] for x in h[5]]
                x0m = h[1] + 0.5 * h[3]
                lamm = h[2] + 0.5 * h[4]
                dd = h[7]
                dn = math.sqrt(sum(x * x for x in dd))
                m = mods.get(int(h[0]))
                if m is None:
                    nomod[int(h[0]) >> 25] += 1
                    lay, rlm, ic = -1, 0.0, 1.0
                else:
                    lay, n_, rlm = m
                    cosv = abs(sum(dd[c] * n_[c] for c in range(3))) / dn if dn > 0 else 1.0
                    ic = 1.0 / max(cosv, 1e-3)
                rows.append((lay, pos, x0m, lamm, rlm, grid_rl(pos[2], math.hypot(pos[0], pos[1])), ic))
            tot = collections.Counter()
            for i in range(len(rows) - 1):
                a, b = rows[i], rows[i + 1]
                vals = dict(eta=eta, pt=pt, pdg=q.m_pdg, istep=i, lay0=a[0], lay1=b[0],
                            r0=math.hypot(a[1][0], a[1][1]), z0=a[1][2], r1=math.hypot(b[1][0], b[1][1]), z1=b[1][2],
                            x0=b[2] - a[2], lam=b[3] - a[3], rl_mod=b[4] * b[6], rl_grid=b[5] * b[6], invcos=b[6])
                for c in scols:
                    sb[c][0] = vals[c]
                S.Fill()
                for c in ("x0", "lam", "rl_mod", "rl_grid"):
                    tot[c] += vals[c]
            kv = dict(eta=eta, pt=pt, pdg=q.m_pdg, nhits=len(rows), x0_first=rows[0][2], x0=tot["x0"], lam=tot["lam"],
                      rl_mod=tot["rl_mod"], rl_grid=tot["rl_grid"], inel=inel, death_r=dr, death_z=dz,
                      nmiss=sum(1 for r_ in rows if r_[0] < 0))
            for c in kcols:
                kb[c][0] = kv[c]
            Kt.Fill()
        f.Close()
    print(f"{Kt.GetEntries()} tracks, {S.GetEntries()} steps; hits whose DetId is not in mkFit's geometry, by "
          f"detector (DetId >> 25): {dict(nomod)}")
    out.cd()
    S.Write()
    Kt.Write()
    out.Close()


main()
