#!/usr/bin/env python3
"""g4s_towers.py FIELDTABLE FILE...: how much of a single particle's visible HGCal energy lands in the tower its track
points to, on the HGCal tower grid of the pointing-tower demonstrator (|eta| from 1.305; EM towers 0.0174 x 2pi/360,
hadronic 0.087 x 2pi/72; phi from -pi; towers pointing to the origin, where these samples' vertex is).
Visible energy: the event's sensitive steps in HGCal cells (DetId detector 8, 9, 10), deposit at the step midpoint.
Pointing: (a) the track at the HGCal face, |z| = 320 cm, from the primary's vertex momentum integrated through the
field table (RK4, 1 cm steps, azimuthally symmetric Br, Bz of grid 0); neutral particles go straight; (b) the vertex
direction. Reported per sample: the fraction of visible energy in the pointed hadronic tower and in its 3x3 (median,
mean, and the share of events with at least 90 % in the 3x3), for (a) and (b); for CE-E deposits on the EM grid, the
fraction in 1, 3x3, 5x5 EM towers around (a). Pions are split by a first inelastic interaction before the face."""
import sys, math, statistics, collections, ROOT
import numpy as np
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
ROOT.gSystem.Load("libRootG4Snitch.so")

ETA0, DETA_EM, NPHI_EM = 1.305, (3.045 - 1.305) / 100, 360
DETA_HAD, NPHI_HAD = DETA_EM * 5, 72
ZFACE = 320.0

# Field: grid 0 of the table, r 0..300, z -700..700 cm in 1 cm steps.
df = ROOT.RDataFrame("field", sys.argv[1]).Filter("grid == 0")
a = df.AsNumpy(["r", "z", "br", "bz"])
R = np.round(a["r"]).astype(int); Z = np.round(a["z"]).astype(int) + 700
BR = np.zeros((301, 1401)); BZ = np.zeros((301, 1401))
BR[R, Z] = a["br"]; BZ[R, Z] = a["bz"]
def field(x, y, z):
    r = math.hypot(x, y)
    if r >= 300 or abs(z) >= 700: return 0.0, 0.0, 0.0
    i, j = int(r), int(z + 700); fr, fz = r - i, (z + 700) - j
    def bil(A): return (A[i, j] * (1 - fr) * (1 - fz) + A[i + 1, j] * fr * (1 - fz) + A[i, j + 1] * (1 - fr) * fz + A[i + 1, j + 1] * fr * fz)
    br, bz = bil(BR), bil(BZ)
    return (br * x / r, br * y / r, bz) if r > 0 else (0.0, 0.0, bz)

def face_point(px, py, pz, q):
    """Where the trajectory from the origin crosses |z| = ZFACE, or None."""
    p = math.sqrt(px * px + py * py + pz * pz)
    x = np.zeros(3); u = np.array([px, py, pz]) / p
    if q == 0:
        if abs(u[2]) < 1e-9: return None
        return u * (ZFACE / abs(u[2]))
    k = 0.0029979 * q / p  # du/ds = k (u x B), per cm, B in T
    def deriv(xx, uu):
        b = field(*xx); return uu, k * np.cross(uu, b)
    h = 1.0
    for _ in range(3000):
        x1, u1 = deriv(x, u)
        x2, u2 = deriv(x + 0.5 * h * x1, u + 0.5 * h * u1)
        x3, u3 = deriv(x + 0.5 * h * x2, u + 0.5 * h * u2)
        x4, u4 = deriv(x + h * x3, u + h * u3)
        xn = x + h / 6 * (x1 + 2 * x2 + 2 * x3 + x4); un = u + h / 6 * (u1 + 2 * u2 + 2 * u3 + u4); un /= np.linalg.norm(un)
        if abs(xn[2]) >= ZFACE:
            w = (ZFACE - abs(x[2])) / (abs(xn[2]) - abs(x[2])); return x + w * (xn - x)
        if math.hypot(xn[0], xn[1]) > 300: return None
        x, u = xn, un
    return None

def eta_phi(x, y, z):
    r = math.hypot(x, y); return (math.asinh(z / r) if r > 0 else 99.0), math.atan2(y, x)
def tower(eta, phi, deta, nphi):
    return int(math.floor((abs(eta) - ETA0) / deta)), int(math.floor((phi + math.pi) / (2 * math.pi / nphi))) % nphi
def within(t, t0, n, nphi):
    dphi = (t[1] - t0[1]) % nphi; dphi = min(dphi, nphi - dphi)
    return abs(t[0] - t0[0]) <= n and dphi <= n

CHARGE = {11: -1, -11: 1, 13: -1, -13: 1, 211: 1, -211: -1, 22: 0, 2112: 0}
for fn in sys.argv[2:]:
    sf = SnitchFile(fn)
    procs = {int(x.idx): str(x.name) for x in sf.rows("P")}
    res = collections.defaultdict(list)
    for ev in sf:
        p, s = ev.p, ev.s
        k = p[0].m_daughters_begin; q0 = p[k]; pdg = q0.m_pdg
        cls = "all"
        if abs(pdg) == 211:
            early = "Inelastic" in procs.get(q0.m_death_proc, "") and abs(q0.m_x_death.fZ) < ZFACE
            cls = "pi interacted before HGCal" if early else "pi reaching HGCal"
        fp = face_point(q0.m_p_birth.fX, q0.m_p_birth.fY, q0.m_p_birth.fZ, CHARGE.get(pdg, 0))
        if fp is None: continue
        eta_a, phi_a = eta_phi(*fp)
        eta_b, phi_b = eta_phi(q0.m_p_birth.fX, q0.m_p_birth.fY, q0.m_p_birth.fZ)
        ta, tb = tower(eta_a, phi_a, DETA_HAD, NPHI_HAD), tower(eta_b, phi_b, DETA_HAD, NPHI_HAD)
        tem = tower(eta_a, phi_a, DETA_EM, NPHI_EM)
        etot = e1 = e3 = e3b = 0.0; eem = eem1 = eem3 = eem5 = 0.0
        for ps in s:
            for st in ps.m_steps_sensitive:
                d = st.m_det_id; det = d >> 28
                if det not in (8, 9, 10) or st.m_edep.m_total <= 0: continue
                e = st.m_edep.m_total
                eta, phi = eta_phi(0.5 * (st.m_x_beg.fX + st.m_x_end.fX), 0.5 * (st.m_x_beg.fY + st.m_x_end.fY), 0.5 * (st.m_x_beg.fZ + st.m_x_end.fZ))
                th = tower(eta, phi, DETA_HAD, NPHI_HAD)
                etot += e
                if th == ta: e1 += e
                if within(th, ta, 1, NPHI_HAD): e3 += e
                if within(th, tb, 1, NPHI_HAD): e3b += e
                if det == 8:
                    tm = tower(eta, phi, DETA_EM, NPHI_EM); eem += e
                    if tm == tem: eem1 += e
                    if within(tm, tem, 1, NPHI_EM): eem3 += e
                    if within(tm, tem, 2, NPHI_EM): eem5 += e
        if etot <= 0: continue
        r = (e1 / etot, e3 / etot, e3b / etot, (eem1 / eem, eem3 / eem, eem5 / eem) if eem > 0 else None, etot)
        res[cls].append(r)
        if cls != "all": res["pi all"].append(r)
    name = fn.rstrip("/").split("/")[-2]
    for cls, rows in sorted(res.items()):
        n = len(rows)
        med = lambda i: statistics.median(r[i] for r in rows); mean = lambda i: statistics.mean(r[i] for r in rows)
        ge90 = sum(1 for r in rows if r[1] >= 0.9) / n; ge90b = sum(1 for r in rows if r[2] >= 0.9) / n
        em = [r[3] for r in rows if r[3]]
        ems = ("  EM grid (CE-E) median 1/3x3/5x5: " + " ".join(f"{statistics.median(x[i] for x in em):.2f}" for i in range(3))) if em else ""
        print(f"{name:12s} {cls:28s} n {n:4d}  E_vis med {statistics.median(r[4] for r in rows):8.1f} MeV | track: 1 tower med {med(0):.2f} mean {mean(0):.2f}, 3x3 med {med(1):.2f} mean {mean(1):.2f}, >=90% in 3x3: {100 * ge90:3.0f}% | vertex dir 3x3 med {med(2):.2f}, >=90%: {100 * ge90b:3.0f}%" + ems)
