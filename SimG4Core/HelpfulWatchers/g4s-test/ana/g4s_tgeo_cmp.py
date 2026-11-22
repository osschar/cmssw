#!/usr/bin/env python3
"""g4s_tgeo_cmp.py ENVELOPES DIFFUSE FILE...: the TGeo material session's maps against Geant4 along true tracks.

ENVELOPES: the output of that session's mb_ladder_layers(), mb_discs("ALL") and mb_envelopes() (lines "L name rmin
rmax zhalf", "D name rmin rmax z halfz", "E r0 r1 z0 z1"), printed from /baz/matevz/root-dev/geo-stuff/mat_budget.C
at mb_disc_zpad = 0.5 (`envs.C` beside this script's notes). DIFFUSE: its diffuse table, matdiffuse_rz_zpad05.root
(TH2D diffuse_x0 / diffuse_lam, density per cm over (z, r), material outside every envelope, averaged over phi).
Both from cmsSimGeom-Run4D123.root (tracker T35, as D128 here).

Per track (one charged primary per event, the outgoing leg in the tracker, between the start of its first tracker
hit and the end of its last, as in g4s_tkmat.py), each own Geant4 step is cut into pieces of at most 0.25 cm, and
each piece is classified by its mid point (r, z) against the envelopes:
  g4_in, g4_out : Geant4 X0 (and lambda_I) inside any envelope / outside all of them
  tg_out        : the diffuse table integrated over the same pieces (its phi average along this track's (r, z))
and per disc crossing (the pieces inside one disc envelope, anywhere on the outgoing leg, so that the discs before
the first and after the last hit are crossed completely; before 2026-11-22 15:00 only the pieces between the first
and the last hit were used, which cut those crossings short), Geant4's X0 inside it times |t_z| (normal incidence),
by family, as a plain mean over crossings and area-weighted over r (1 cm rings, weight r), the way the TGeo session
averages its disc budgets (a track hits a random point of the disc)."""
import sys, math, collections
import numpy as np
import ROOT
ROOT.gSystem.Load("libRootG4Snitch.so")

envf, difff = sys.argv[1], sys.argv[2]
L, D, E = [], [], []
for line in open(envf):
    a = line.split()
    if a[0] == "L":
        L.append((a[1], float(a[2]), float(a[3]), float(a[4])))
    elif a[0] == "D":
        D.append((a[1], float(a[2]), float(a[3]), float(a[4]), float(a[5])))
    elif a[0] == "E":
        E.append(tuple(float(x) for x in a[1:5]))
names = [l[0] for l in L]
for d in D:
    names += [d[0] + "+", d[0] + "-"]
assert len(names) == len(E), (len(names), len(E))
E = np.array(E)
disc_idx = [i for i, n in enumerate(names) if not any(n == l[0] for l in L)]

fd = ROOT.TFile.Open(difff)
zpad = fd.Get("zpad")
print(f"diffuse table {difff}, zpad tag {zpad.GetTitle() if zpad else 'missing'}; {len(E)} envelopes from {envf}")
hx, hl = fd.Get("diffuse_x0"), fd.Get("diffuse_lam")
nz, nr = hx.GetNbinsX(), hx.GetNbinsY()
zmin, zmax = hx.GetXaxis().GetXmin(), hx.GetXaxis().GetXmax()
rmin, rmax = hx.GetYaxis().GetXmin(), hx.GetYaxis().GetXmax()
TX = np.array([[hx.GetBinContent(i + 1, j + 1) for j in range(nr)] for i in range(nz)])
TL = np.array([[hl.GetBinContent(i + 1, j + 1) for j in range(nr)] for i in range(nz)])

edges = [0.0, 0.5, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5]
acc = collections.defaultdict(list)
cross = collections.defaultdict(list)  # family -> [(r, x0 * |t_z|)]
for fn in sys.argv[3:]:
    f = ROOT.TFile.Open(fn)
    sdcls = [str(getattr(d, "class")) for d in f.Get("D")]
    trk_sd = {i for i, c in enumerate(sdcls) if "TkAccumulating" in c}
    V = [(float(v.x0), float(getattr(v, "lambda")), int(v.sd) in trk_sd) for v in f.Get("V")]
    for ev in f.Get("T"):
        p, s = ev.p, ev.s
        k = p[0].m_daughters_begin
        q = p[k]
        pb = q.m_p_birth
        pt = math.hypot(pb.fX, pb.fY)
        if q.m_charge == 0:
            continue
        eta = math.asinh(pb.fZ / pt)
        gid = q.m_g4_id
        tcut, rmx = 1e30, 0.0
        for pp in s[k].m_path:
            r = math.hypot(pp.m_x.fX, pp.m_x.fY)
            rmx = max(rmx, r)
            if r > 118.0 or abs(pp.m_x.fZ) > 282.0 or r < rmx - 1.0:
                tcut = pp.m_x.fT
                break
        own = sorted([st for st in list(s[k].m_steps_sensitive) + list(s[k].m_steps_inert)
                      if st.m_g4_pid == gid and st.m_x_beg.fT < tcut], key=lambda st: st.m_x_beg.fT)
        hit_i = [i for i, st in enumerate(own) if st.m_det_id != 0 and 0 <= st.m_volume < len(V) and V[st.m_volume][2]]
        if len(hit_i) < 2:
            continue
        # pieces: mid point, length, 1/X0, 1/lambda, |t_z|, between the first and the last hit (1) or not (0); the
        # whole outgoing leg, so that disc crossings before the first and after the last hit are complete
        P = []
        for ist, st in enumerate(own):
            a = (st.m_x_beg.fX, st.m_x_beg.fY, st.m_x_beg.fZ)
            b = (st.m_x_end.fX, st.m_x_end.fY, st.m_x_end.fZ)
            d = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
            ln = math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
            x0, lam, _ = V[st.m_volume]
            if ln <= 0 or x0 <= 0:
                continue
            n = max(1, int(math.ceil(ln / 0.25)))
            tz = abs(d[2]) / ln
            for j in range(n):
                u = (j + 0.5) / n
                P.append((a[0] + u * d[0], a[1] + u * d[1], a[2] + u * d[2], ln / n, 1.0 / x0,
                          1.0 / lam if lam > 0 else 0.0, tz, 1.0 if hit_i[0] <= ist <= hit_i[-1] else 0.0))
        P = np.array(P)
        r = np.hypot(P[:, 0], P[:, 1])
        z = P[:, 2]
        dl, ix0, ilam, tz = P[:, 3], P[:, 4], P[:, 5], P[:, 6]
        btw = P[:, 7] > 0
        inside = np.zeros(len(P), dtype=bool)
        which = np.full(len(P), -1)
        for ie, e in enumerate(E):
            m = (r >= e[0]) & (r <= e[1]) & (z >= e[2]) & (z <= e[3])
            inside |= m
            which[m & (which < 0)] = ie
        iz = np.clip(((z - zmin) / (zmax - zmin) * nz).astype(int), 0, nz - 1)
        ir = np.clip(((r - rmin) / (rmax - rmin) * nr).astype(int), 0, nr - 1)
        g4x = dl * ix0
        g4l = dl * ilam
        tgx, tgl = dl * TX[iz, ir], dl * TL[iz, ir]
        row = dict(g4_all=g4x[btw].sum(), g4_in=g4x[btw & inside].sum(), g4_out=g4x[btw & ~inside].sum(),
                   tg_out=tgx[btw].sum(), g4l_out=g4l[btw & ~inside].sum(), tgl_out=tgl[btw].sum(),
                   tg_on_in=tgx[btw & inside].sum())
        for i in range(len(edges) - 1):
            if edges[i] <= abs(eta) < edges[i + 1]:
                acc[i].append(row)
        for ie in disc_idx:
            m = which == ie
            if m.any():
                fam = "".join(c for c in names[ie] if c.isalpha())
                cross[fam].append((float(np.average(r[m], weights=dl[m])), float((g4x[m] * tz[m]).sum())))
    f.Close()


def mean(v):
    return sum(v) / len(v) if v else float("nan")


print("X0 per track between the first and the last hit: Geant4 inside / outside the envelopes, and the TGeo diffuse "
      "table along the same path")
print(f"{'|eta|':>9} {'n':>5} {'g4 all':>7} {'g4 in':>7} {'g4 out':>7} {'TGeo out':>8} {'TG/g4':>6} {'in %':>5} "
      f"{'lam g4':>7} {'lam TG':>7} {'TG/g4':>6} {'TG on in':>8}")
for i in sorted(acc):
    a = acc[i]
    m = {c: mean([x[c] for x in a]) for c in a[0]}
    print(f"{edges[i]:.1f}-{edges[i + 1]:.1f} {len(a):5d} {m['g4_all']:7.3f} {m['g4_in']:7.3f} {m['g4_out']:7.3f} "
          f"{m['tg_out']:8.3f} {m['tg_out'] / m['g4_out']:6.2f} {100 * m['g4_in'] / m['g4_all']:5.1f} "
          f"{m['g4l_out']:7.4f} {m['tgl_out']:7.4f} {m['tgl_out'] / m['g4l_out']:6.2f} {m['tg_on_in']:8.4f}")
print("\nper disc crossing, Geant4 X0 inside the disc envelope (zpad 0.5) times |t_z|: plain mean, and area-weighted "
      "over 1 cm rings (weight r)")
for fam in ("TFPX", "TEPX", "TEDD"):
    v = cross.get(fam, [])
    if not v:
        continue
    rings = collections.defaultdict(list)
    for rr, x in v:
        rings[int(rr)].append(x)
    wsum = sum(rr + 0.5 for rr in rings)
    aw = sum((rr + 0.5) * mean(xs) for rr, xs in rings.items()) / wsum
    print(f"  {fam}: {len(v)} crossings, r {min(x[0] for x in v):.1f}-{max(x[0] for x in v):.1f} cm, mean "
          f"{mean([x[1] for x in v]):.4f}, area-weighted {aw:.4f}, rings with crossings {len(rings)}")
