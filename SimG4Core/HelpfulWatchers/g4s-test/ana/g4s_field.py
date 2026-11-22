#!/usr/bin/env python3
"""g4s_field.py OUT.root FILE...: mkFit's parametrised magnetic field against the CMSSW field map, by the bending
integral along the primary's true path between its first and its last tracker hit. Needs the field branch "b"
(G4SNITCH_BFIELD=1). One charged primary per event (the single gun).

The path is the snitch's thinned path (points every <= 5 cm or 0.05 rad) with the tracker hits inserted by time; a
hit is the mid point of the primary's own sensitive steps in one tracker DetId (consecutive steps merged). The map's
field at a hit is interpolated by time between the neighbouring path points. Each segment a -> b contributes
F = t x B dl, t the chord direction, B the mean of its end points (or the model's sampled value, below). Projected:
  G_phi = sum -(F . n_phi), n_phi = z x t_T / s, s = |t_T|: the azimuthal bending, T cm; for a pure Bz it is the
          integral of Bz over the transverse path length, the turn of the transverse direction is q k / pT G_phi, and
          the curvature (so pT) a fit infers scales as 1 / G_phi;
  G_th  = sum s F . n_th, n_th = the unit vector normal to t in the (t, z) plane: the polar bending, q k / pT G_th,
          which only the radial field makes.
Only the outgoing leg in the tracker counts: hits after the first path point outside r < 118 cm, |z| < 282 cm, or
after the path's radius has fallen 1 cm below its maximum, are dropped.
Field models (z, r in cm, tesla; the constants of cmssw_20_v2p2 MkFitCore/src/Config.cc and MkFitGeometryESProducer):
  const : Bz = 3.8112 (Config::Bfield)
  build : Bz = (b0 z^2 + b1 z + c1)(a r^2 + 1), ParabolicMf's constants (Config::mag_*), no Br: track building
  fit   : the same form with bFieldParams (Config::mag_fit_*), plus Br = -(r/2) dBz/dz: the final fit
and as mkFit samples them per step between consecutive hits:
  build_step : build's Bz at the step's start hit, constant over the step (finding: propagation, mini-propagators)
  fit_mid    : fit's Bz at the chord mid point of the step, constant over it (refitBFieldAtMid, without
               refitRadialFieldCorr)
  fit_s<n>   : the final fit's scheme with n sub-steps of equal path length per step (refitBkwSubSteps: 2 in the
               backward pass, 1 in the forward one): per sub-step Bz at its chord mid point (refitBFieldAtMid) plus
               the radial-field correction as FinalFit.cc applies it (refitRadialFieldCorr): the change of r p_phi that
               the constant-Bc helix misses against the parametrised flux, D / (q k) = dPsi/2pi - Bc d(r^2/2) between
               the sub-step's end points, half as a turn D / (2 r0 pT) at its start, half as D / (2 r1 pT) at its
               end. The sub-step ends are taken on the true path.
Per track: dphi_<m> = G_phi(m) / G_phi(map) - 1, the relative pT bias the model gives; dth_<m> = (G_th(m) -
G_th(map)) / G_phi(map), the polar-angle error relative to the azimuthal bending; res_<m>: the rms over the hits of
the azimuthal displacement left when the map's displacement profile D(L) = (k / pT) sum g (L - l) is fitted by
a + b L + c D_m(L) (offset, angle and curvature free), in micron at pT = 1 GeV (it scales as 1 / pT); dpt_<m> = 1 / c - 1,
the relative pT bias of that fit (unweighted, no scattering): the measure closest to what a track fit infers.
kmax*: the step of fit_s1 with the largest kick turn (over G_phi(map)), its radii and z."""
import sys, math, array, ROOT
ROOT.gSystem.Load("libRootG4Snitch.so")

K = 0.299792458e-2  # GeV / (T cm)
C1, B0, B1, A = 3.8114, -3.94991e-06, 7.53701e-06, 2.43878e-11  # ParabolicMf
FC1, FB0, FB1, FA = 3.81036, -2.03767e-06, 7.34495e-06, 3.01291e-07  # bFieldParams


def bz_par(z, r, c1, b0, b1, a):
    return (b0 * z * z + b1 * z + c1) * (a * r * r + 1.0)


def b_build(x, y, z):
    return (0.0, 0.0, bz_par(z, math.hypot(x, y), C1, B0, B1, A))


def b_fit(x, y, z):
    r2 = x * x + y * y
    bz = (FB0 * z * z + FB1 * z + FC1) * (FA * r2 + 1.0)
    dbzdz = (2 * FB0 * z + FB1) * (FA * r2 + 1.0)
    # Br = -(r/2) dBz/dz, so (Bx, By) = -(1/2) dBz/dz (x, y)
    return (-0.5 * dbzdz * x, -0.5 * dbzdz * y, bz)


def flux_fit(x, y, z):
    """the flux function of FinalFit.cc, Psi / 2 pi = Z(z) (a r^4 / 4 + r^2 / 2), with the fit's constants"""
    r2 = x * x + y * y
    return (FB0 * z * z + FB1 * z + FC1) * (0.25 * FA * r2 * r2 + 0.5 * r2)


def b_const(x, y, z):
    return (0.0, 0.0, 3.8112)


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def proj(t, F):
    s = math.hypot(t[0], t[1])
    nphi = (-t[1] / s, t[0] / s, 0.0)
    gphi = -(F[0] * nphi[0] + F[1] * nphi[1])
    tz = t[2]
    nth = (-tz * t[0], -tz * t[1], 1.0 - tz * tz)
    nn = math.sqrt(nth[0] ** 2 + nth[1] ** 2 + nth[2] ** 2)
    gth = s * (F[0] * nth[0] + F[1] * nth[1] + F[2] * nth[2]) / nn if nn > 1e-9 else 0.0
    return gphi, gth


NSUB = (1, 2, 4)
MODELS = ["const", "build", "fit", "build_step", "fit_mid"] + [f"fit_s{n}" for n in NSUB]


def shape_res(Lh, gseg, lseg, gseg_m):
    """rms residual (cm, per k/pT = 1) of D_map(L_i) fitted by a + b L + c D_m(L) over the hits, and c."""
    def D(g, L):
        return sum(gi * (L - li) for gi, li in zip(g, lseg) if li < L)
    dm = [D(gseg, L) for L in Lh]
    dx = [D(gseg_m, L) for L in Lh]
    n = len(Lh)
    # least squares on [1, L, dx]
    import numpy as np
    M = np.array([[1.0, L, x] for L, x in zip(Lh, dx)])
    y = np.array(dm)
    coef, *_ = np.linalg.lstsq(M, y, rcond=None)
    r = y - M @ coef
    return math.sqrt(float((r * r).sum()) / n), float(coef[2])


def main():
    out = ROOT.TFile(sys.argv[1], "RECREATE")
    nt = ROOT.TTree("F", "field integrals per track")
    cols = ["pt", "eta", "phi", "nhits", "rlast", "zlast", "ltr", "gphi_map", "gth_map", "bzmean_map"]
    for m in MODELS:
        cols += ["dphi_" + m, "dth_" + m, "res_" + m, "dpt_" + m]
    cols += ["kmax", "kmax_r0", "kmax_r1", "kmax_z0", "kmax_z1"]
    buf = {c: array.array("f", [0.0]) for c in cols}
    for c in cols:
        nt.Branch(c, buf[c], c + "/F")
    nskip_ins = [0]
    ndrop = [0]
    for fn in sys.argv[2:]:
        f = ROOT.TFile.Open(fn)
        t = f.Get("T")
        sdcls = [str(getattr(d, "class")) for d in f.Get("D")]
        vsd = [int(v.sd) for v in f.Get("V")]
        trk_sd = {i for i, c in enumerate(sdcls) if "TkAccumulating" in c}
        trk_vol = {i for i, s in enumerate(vsd) if s in trk_sd}
        for ev in t:
            p, s, bb = ev.p, ev.s, ev.b
            k = p[0].m_daughters_begin
            q = p[k]
            if q.m_charge == 0:
                continue
            gid = q.m_g4_id
            # hits: own sensitive tracker steps, merged per consecutive DetId
            hits = []
            for st in s[k].m_steps_sensitive:
                if st.m_g4_pid != gid or st.m_volume not in trk_vol or st.m_det_id == 0:
                    continue
                xm = (0.5 * (st.m_x_beg.fX + st.m_x_end.fX), 0.5 * (st.m_x_beg.fY + st.m_x_end.fY),
                      0.5 * (st.m_x_beg.fZ + st.m_x_end.fZ), 0.5 * (st.m_x_beg.fT + st.m_x_end.fT))
                if hits and hits[-1][0] == st.m_det_id:
                    continue
                hits.append((st.m_det_id, xm))
            hits.sort(key=lambda h: h[1][3])
            path = s[k].m_path
            tcut, rmax = 1e30, 0.0
            for pp in path:
                r = math.hypot(pp.m_x.fX, pp.m_x.fY)
                rmax = max(rmax, r)
                if r > 118.0 or abs(pp.m_x.fZ) > 282.0 or r < rmax - 1.0:
                    tcut = pp.m_x.fT
                    break
            nh_all = len(hits)
            hits = [h for h in hits if h[1][3] < tcut]
            ndrop[0] += nh_all - len(hits)
            if len(hits) < 4:
                continue
            bv = bb[k].m_b
            pts = [((pp.m_x.fX, pp.m_x.fY, pp.m_x.fZ), pp.m_x.fT, (bv[i].fX, bv[i].fY, bv[i].fZ), -1)
                   for i, pp in enumerate(path)]
            t0, t1 = hits[0][1][3], hits[-1][1][3]
            # insert hits, with the map's B interpolated in time
            ok = True
            for ih, (_, xm) in enumerate(hits):
                th = xm[3]
                j = next((j for j in range(len(pts) - 1) if pts[j][1] <= th <= pts[j + 1][1]), None)
                if j is None:
                    ok = False
                    break
                ta, tb = pts[j][1], pts[j + 1][1]
                w = (th - ta) / (tb - ta) if tb > ta else 0.0
                B = tuple(pts[j][2][c] + w * (pts[j + 1][2][c] - pts[j][2][c]) for c in range(3))
                pts.insert(j + 1, (xm[:3], th, B, ih))
            if not ok:
                nskip_ins[0] += 1
                continue
            seq = [pp for pp in pts if t0 <= pp[1] <= t1]
            hpos = [h[1][:3] for h in hits]
            # cumulative 3D and transverse path length at each point of seq, and the hits' places in seq
            cl, cT = [0.0], [0.0]
            for u, v in zip(seq, seq[1:]):
                d = [v[0][c] - u[0][c] for c in range(3)]
                cl.append(cl[-1] + math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2))
                cT.append(cT[-1] + math.hypot(d[0], d[1]))
            hidx = [j for j, pp in enumerate(seq) if pp[3] >= 0]

            def at_len(x):
                """position and transverse length at 3D path length x along seq"""
                j = max(0, min(len(cl) - 2, next((j for j in range(len(cl) - 1) if cl[j + 1] >= x), len(cl) - 2)))
                w = (x - cl[j]) / (cl[j + 1] - cl[j]) if cl[j + 1] > cl[j] else 0.0
                pos = [seq[j][0][c] + w * (seq[j + 1][0][c] - seq[j][0][c]) for c in range(3)]
                return pos, cT[j] + w * (cT[j + 1] - cT[j])

            # the final fit's scheme with n sub-steps of equal path length per hit step: for each sub-step its
            # length range, Bc at its chord mid point, and the two half-kicks (transverse position, turn)
            subs = {n: [] for n in NSUB}
            kicks = {n: [] for n in NSUB}
            kmax = [0.0, 0.0, 0.0, 0.0, 0.0]
            for n in NSUB:
                for i in range(len(hidx) - 1):
                    la, lb = cl[hidx[i]], cl[hidx[i + 1]]
                    ends = [at_len(la + (lb - la) * jj / n) for jj in range(n + 1)]
                    for jj in range(n):
                        (p0, lt0), (p1, lt1) = ends[jj], ends[jj + 1]
                        bc = b_fit(*[0.5 * (p0[c] + p1[c]) for c in range(3)])[2]
                        subs[n].append((la + (lb - la) * jj / n, la + (lb - la) * (jj + 1) / n, bc))
                        r0, r1 = math.hypot(p0[0], p0[1]), math.hypot(p1[0], p1[1])
                        kk = (flux_fit(*p1) - flux_fit(*p0)) - bc * 0.5 * (r1 * r1 - r0 * r0)
                        kicks[n].append((lt0, 0.5 * kk / r0))
                        kicks[n].append((lt1, 0.5 * kk / r1))
                        if n == 1 and abs(0.5 * kk * (1 / r0 + 1 / r1)) > abs(kmax[0]):
                            kmax[:] = [0.5 * kk * (1 / r0 + 1 / r1), r0, r1, p0[2], p1[2]]

            def bc_at(n, x):
                for (x0, x1, bc) in subs[n]:
                    if x0 <= x <= x1:
                        return bc
                return subs[n][-1][2]

            G = {m: [0.0, 0.0] for m in ["map"] + MODELS}
            gseg = {m: [] for m in ["map"] + MODELS}
            lseg = []
            Lh = [cT[j] for j in hidx]
            bzsum = 0.0
            cur = 0
            for j, (a, b) in enumerate(zip(seq, seq[1:])):
                if a[3] >= 0:
                    cur = a[3]
                h0 = hpos[cur]
                d = [b[0][c] - a[0][c] for c in range(3)]
                dl = cl[j + 1] - cl[j]
                dT = cT[j + 1] - cT[j]
                if dl < 1e-6 or dT < 1e-6:
                    continue
                tt = [x / dl for x in d]
                xm = 0.5 * (cl[j] + cl[j + 1])
                Bm = {}
                Bm["map"] = [0.5 * (a[2][c] + b[2][c]) for c in range(3)]
                for name, fn_ in (("const", b_const), ("build", b_build), ("fit", b_fit)):
                    ba, bb_ = fn_(*a[0]), fn_(*b[0])
                    Bm[name] = [0.5 * (ba[c] + bb_[c]) for c in range(3)]
                Bm["build_step"] = b_build(*h0)
                Bm["fit_mid"] = (0.0, 0.0, bc_at(1, xm))
                for n in NSUB:
                    Bm[f"fit_s{n}"] = (0.0, 0.0, bc_at(n, xm))
                for m, B in Bm.items():
                    F = [x * dl for x in cross(tt, B)]
                    gp, gt = proj(tt, F)
                    G[m][0] += gp
                    G[m][1] += gt
                    gseg[m].append(gp)
                bzsum += Bm["map"][2] * dT
                lseg.append(0.5 * (cT[j] + cT[j + 1]))
            # the kicks, as zero-length entries
            for n in NSUB:
                m = f"fit_s{n}"
                for (lt, g) in kicks[n]:
                    G[m][0] += g
                    for mm in gseg:
                        gseg[mm].append(g if mm == m else 0.0)
                    lseg.append(lt)
            lT = cT[-1]
            gm = G["map"][0]
            if abs(gm) < 1e-6:
                continue
            pb = q.m_p_birth
            pt = math.hypot(pb.fX, pb.fY)
            buf["pt"][0] = pt
            buf["eta"][0] = math.asinh(pb.fZ / pt)
            buf["phi"][0] = math.atan2(pb.fY, pb.fX)
            buf["nhits"][0] = len(hits)
            buf["rlast"][0] = math.hypot(hpos[-1][0], hpos[-1][1])
            buf["zlast"][0] = hpos[-1][2]
            buf["ltr"][0] = lT
            buf["gphi_map"][0] = gm
            buf["gth_map"][0] = G["map"][1]
            buf["bzmean_map"][0] = bzsum / lT if lT > 0 else 0
            for m in MODELS:
                buf["dphi_" + m][0] = G[m][0] / gm - 1.0
                buf["dth_" + m][0] = (G[m][1] - G["map"][1]) / gm
                rr, cc = shape_res(Lh, gseg["map"], lseg, gseg[m])
                buf["res_" + m][0] = 1e4 * K * rr
                buf["dpt_" + m][0] = 1.0 / cc - 1.0 if cc else -9
            buf["kmax"][0] = kmax[0] / gm
            for j, c in enumerate(["kmax_r0", "kmax_r1", "kmax_z0", "kmax_z1"]):
                buf[c][0] = kmax[1 + j]
            nt.Fill()
        f.Close()
    print(f"{nt.GetEntries()} tracks; skipped for a hit outside the path's time range: {nskip_ins[0]}; hits dropped "
          f"after the outgoing leg: {ndrop[0]}")
    out.cd()
    nt.Write()
    out.Close()


main()
