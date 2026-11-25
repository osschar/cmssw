#!/usr/bin/env python3
"""g4s_lowpt.py [--field FIELDTABLE] FILE...: low-pT single particles from the origin (G4SNITCH_GUN=single): does the
primary reach the HGCal face, and how much of the event's visible HGCal energy a track follower could associate with
it, against pointing towers. Truth only: Geant4 deposits, no thresholds, no noise, no pileup.

The primary's path is its leading chain, as in g4s_a2.py: the primary, continued by a daughter of the same kind born
where it died with more than half of its last energy; a hadron's chain ends at its first inelastic interaction. The
path is the chain's thinned G4Snitch path (m_path), interpolated linearly between points.

Per event:
- the face: the chain's first crossing of |z| = 320 cm, its radius, |eta| and pT there. Classes: "reach" (crossing
  at |eta| <= 3.045, the tower grid's edge, r >= 30.3 cm), "hole" (crossing at smaller radius, into the beam-pipe
  hole), and no crossing; pions without a crossing split into "pi inelastic before face" (the primary's death process
  contains "Inelastic", at |z| < 320) and "stopped before face" (anything else). For no crossing: the death process
  and position of the chain's last particle.
- visible HGCal energy: the event's sensitive steps in HGCal cells (DetId detector 8, 9, 10) with a deposit above 0,
  each at its step midpoint (as g4s_towers.py).
- towers: hadronic pointing towers, |eta| from 1.305, 0.087 x 2pi/72, pointing at the origin (as g4s_towers.py).
  (a)  the 3x3 around the tower of the chain's true crossing of the face;
  (ar) the 3x3 around the RK4 crossing of the face from the vertex momentum, without energy loss (g4s_towers.py's
       pointing; needs --field), the only face pointing for events without a true crossing;
  (av) the 3x3 around the vertex direction.
- follower ceiling: for each deposit at depth z, the path's position at that z (if the path crosses that z more than
  once, the crossing nearest to the deposit in x-y; deposits beyond the path's z range use its last point).
  (b)  the 3x3 towers around the tower of that position;
  (c3, c5, c10) the deposit within 3, 5, 10 cm of that position in x-y;
  (bx, x5, x10) as (b), c5, c10, but deposits beyond the path's z range use the straight continuation of the path's
       last segment (its direction where the chain ended, e.g. at a pion's inelastic interaction) to the deposit's
       z, in place of the last point: a follower that carries the track's direction into the shower.
  (beyond) the share of the visible energy deeper (in |z|) than the path reaches.
Each is the fraction of the event's visible HGCal energy. Reported per sample and class: medians and means over the
events with visible energy, and the share of events with at least 90 % in the 3x3."""
import sys, math, statistics, collections, argparse, time, ROOT
import numpy as np
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
ROOT.gSystem.Load("libRootG4Snitch.so")
ROOT.gInterpreter.Declare(r'''
#include "SimG4Core/HelpfulWatchers/interface/G4SnitchDataFormat.h"
// The event's sensitive HGCal steps with a deposit, 5 numbers each: the step midpoint x, y, z (cm), the deposit (MeV),
// the DetId detector (8 CE-E silicon, 9 CE-H silicon, 10 CE-H scintillator).
std::vector<double> g4s_lowpt_hgcal_deposits(const std::vector<G4S_ParticleSteps>& s) {
  std::vector<double> out;
  for (auto &ps : s)
    for (auto &st : ps.m_steps_sensitive) {
      unsigned det = st.m_det_id >> 28;
      if ((det == 8 || det == 9 || det == 10) && st.m_edep.m_total > 0) {
        out.push_back(0.5 * (st.m_x_beg.fX + st.m_x_end.fX));
        out.push_back(0.5 * (st.m_x_beg.fY + st.m_x_end.fY));
        out.push_back(0.5 * (st.m_x_beg.fZ + st.m_x_end.fZ));
        out.push_back(st.m_edep.m_total);
        out.push_back(det);
      }
    }
  return out;
}
''')

ETA0, DETA_HAD, NPHI_HAD = 1.305, (3.045 - 1.305) / 20, 72
ETA_EDGE = 3.045
ZFACE = 320.0

ap = argparse.ArgumentParser()
ap.add_argument("--field", help="field table (fieldtable/magFieldTable_D128.root) for the RK4 face pointing (ar)")
ap.add_argument("files", nargs="+")
args = ap.parse_args()

# ---- RK4 face pointing, copied from g4s_towers.py (grid 0 of the table, r 0..300, z -700..700 cm, 1 cm steps).
if args.field:
    df = ROOT.RDataFrame("field", args.field).Filter("grid == 0")
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
def face_point_rk4(px, py, pz, q):
    p = math.sqrt(px * px + py * py + pz * pz)
    x = np.zeros(3); u = np.array([px, py, pz]) / p
    if q == 0:
        if abs(u[2]) < 1e-9: return None
        return u * (ZFACE / abs(u[2]))
    k = 0.0029979 * q / p
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

# ---- towers
def towers(x, y, z):
    """Hadronic tower indices (eta, phi) and z side of points, numpy arrays."""
    r = np.maximum(np.hypot(x, y), 1e-9)
    eta = np.arcsinh(z / r); phi = np.arctan2(y, x)
    ie = np.floor((np.abs(eta) - ETA0) / DETA_HAD).astype(np.int64)
    ip = np.floor((phi + math.pi) / (2 * math.pi / NPHI_HAD)).astype(np.int64) % NPHI_HAD
    return ie, ip, np.sign(z)
def in3x3(t, t0):
    dphi = (t[1] - t0[1]) % NPHI_HAD; dphi = np.minimum(dphi, NPHI_HAD - dphi)
    return (np.abs(t[0] - t0[0]) <= 1) & (dphi <= 1) & (t[2] == t0[2])
def pt_tower(x, y, z):
    ie, ip, sg = towers(np.array([x]), np.array([y]), np.array([z])); return ie[0], ip[0], sg[0]

CHARGE = {11: -1, -11: 1, 13: -1, -13: 1, 211: 1, -211: -1, 22: 0, 2112: 0}

def chain(p, s, k, procs):
    """The leading chain of particle k: path points (n x 3), momenta (n x 4), the last particle and its death process."""
    kind = abs(p[k].m_pdg); X = []; P = []; cur = k
    while True:
        q = p[cur]; path = s[cur].m_path
        for pp in path:
            X.append((pp.m_x.fX, pp.m_x.fY, pp.m_x.fZ)); P.append((pp.m_p.fX, pp.m_p.fY, pp.m_p.fZ, pp.m_p.fT))
        dproc = procs.get(q.m_death_proc, "?")
        if kind not in (11, 13, 22) and "Inelastic" in dproc: break
        e_last = path[path.size() - 1].m_p.fT if path.size() else q.m_p_birth.fT
        nxt = None; best = 0.0
        for j in range(q.m_daughters_begin, q.m_daughters_end):
            d = p[j]
            if abs(d.m_pdg) != kind: continue
            dx = math.dist((d.m_x_birth.fX, d.m_x_birth.fY, d.m_x_birth.fZ), (q.m_x_death.fX, q.m_x_death.fY, q.m_x_death.fZ))
            if dx < 1e-3 and d.m_p_birth.fT > 0.5 * e_last and d.m_p_birth.fT > best: nxt, best = j, d.m_p_birth.fT
        if nxt is None: break
        cur = nxt
    return np.array(X, dtype=float).reshape(-1, 3), np.array(P, dtype=float).reshape(-1, 4), p[cur], dproc

def face_crossing(X, P):
    """The first crossing of |z| = ZFACE: position, momentum (linear interpolation), or None."""
    az = np.abs(X[:, 2])
    idx = np.nonzero((az[:-1] < ZFACE) & (az[1:] >= ZFACE))[0]
    if idx.size == 0: return None
    i = idx[0]; w = (ZFACE - az[i]) / (az[i + 1] - az[i])
    return X[i] + w * (X[i + 1] - X[i]), P[i] + w * (P[i + 1] - P[i])

def path_at_deposits(X, dep):
    """For each deposit (N x 3), the path position at the deposit's z: among the path segments that bracket that z,
    the crossing nearest to the deposit in x-y; the path's last point where no segment brackets it. Returns the
    positions (N x 3) and the mask of the deposits with a bracketing segment."""
    N = dep.shape[0]
    out = np.repeat(X[-1:], N, axis=0); has = np.zeros(N, dtype=bool)
    if X.shape[0] < 2: return out, has
    zr = np.round(dep[:, 2] / 0.05) * 0.05  # deposits share z planes; 0.5 mm changes the crossing by < 0.5 mm
    zu, inv = np.unique(zr, return_inverse=True)
    z1, z2 = X[:-1, 2], X[1:, 2]
    br = ((z1[None, :] - zu[:, None]) * (z2[None, :] - zu[:, None]) <= 0) & (z1 != z2)[None, :]
    ui, si = np.nonzero(br)
    if ui.size == 0: return out, has
    w = (zu[ui] - z1[si]) / (z2[si] - z1[si])
    pts = X[si] + w[:, None] * (X[si + 1] - X[si])
    cnt = np.bincount(ui, minlength=zu.size); K = cnt.max()
    cand = np.full((zu.size, K, 3), np.nan)
    order = np.argsort(ui, kind="stable"); ui_s = ui[order]
    start = np.concatenate(([0], np.cumsum(cnt)[:-1]))
    slot = np.arange(ui_s.size) - start[ui_s]
    cand[ui_s, slot] = pts[order]
    c = cand[inv]  # N x K x 3
    d2 = (c[:, :, 0] - dep[:, None, 0]) ** 2 + (c[:, :, 1] - dep[:, None, 1]) ** 2
    has = ~np.all(np.isnan(d2), axis=1)
    best = np.argmin(np.where(np.isnan(d2), np.inf, d2), axis=1)
    out[has] = c[np.arange(N), best][has]
    return out, has

def extrapolate_last(X, dep, pos, has):
    """pos with the deposits that no segment brackets moved onto the straight continuation of the path's last segment
    of non-zero length in z, at the deposit's z (forward only; otherwise the last point stays)."""
    out = pos.copy()
    dz = np.diff(X[:, 2]); ok = np.nonzero(np.abs(dz) > 1e-6)[0]
    if ok.size == 0: return out
    i = ok[-1]; d = X[i + 1] - X[i]; last = X[-1]
    t = (dep[:, 2] - last[2]) / d[2]
    m = (~has) & (t > 0)
    out[m] = last[None, :] + t[m, None] * d[None, :]
    return out

def eta_of(x, y, z):
    r = math.hypot(x, y); return math.asinh(z / r) if r > 0 else math.copysign(99.0, z)

def analyse(fn):
    sf = SnitchFile(fn)
    procs = {int(x.idx): str(x.name) for x in sf.rows("P")}
    rows = []
    for ev in sf:
        p, s = ev.p, ev.s
        k = p[0].m_daughters_begin; q0 = p[k]; pdg = q0.m_pdg; kind = abs(pdg)
        pt0 = math.hypot(q0.m_p_birth.fX, q0.m_p_birth.fY)
        eta0 = math.asinh(q0.m_p_birth.fZ / pt0) if pt0 > 0 else 99.0
        X, P, ql, dproc = chain(p, s, k, procs)
        fc = face_crossing(X, P)
        r = dict(eta0=eta0, pt0=pt0)
        if fc is not None:
            xf, pf = fc
            r["r_face"] = math.hypot(xf[0], xf[1]); r["eta_face"] = abs(eta_of(*xf)); r["pt_face"] = math.hypot(pf[0], pf[1])
            r["cls"] = "reach" if r["eta_face"] <= ETA_EDGE else "hole"
        else:
            early = kind == 211 and "Inelastic" in procs.get(q0.m_death_proc, "") and abs(q0.m_x_death.fZ) < ZFACE
            r["cls"] = "pi inelastic before face" if early else "stopped before face"
            r["death"] = dproc; r["z_death"] = abs(ql.m_x_death.fZ); r["r_death"] = math.hypot(ql.m_x_death.fX, ql.m_x_death.fY)
            r["ke_death"] = ql.m_p_death.fT - ql.m_mass
        dv = np.asarray(ROOT.g4s_lowpt_hgcal_deposits(s), dtype=float).reshape(-1, 5)
        etot = dv[:, 3].sum() if dv.size else 0.0
        r["evis"] = etot
        if etot > 0:
            x, y, z, e = dv[:, 0], dv[:, 1], dv[:, 2], dv[:, 3]
            tdep = towers(x, y, z)
            if fc is not None:
                r["a"] = e[in3x3(tdep, pt_tower(*fc[0]))].sum() / etot
            if args.field:
                fr = face_point_rk4(q0.m_p_birth.fX, q0.m_p_birth.fY, q0.m_p_birth.fZ, CHARGE.get(pdg, 0))
                if fr is not None: r["ar"] = e[in3x3(tdep, pt_tower(*fr))].sum() / etot
            r["av"] = e[in3x3(tdep, pt_tower(q0.m_p_birth.fX, q0.m_p_birth.fY, q0.m_p_birth.fZ))].sum() / etot
            pa, has = path_at_deposits(X, dv[:, :3])
            tp = towers(pa[:, 0], pa[:, 1], pa[:, 2])
            r["b"] = e[in3x3(tdep, tp)].sum() / etot
            dxy = np.hypot(x - pa[:, 0], y - pa[:, 1])
            for c in (3, 5, 10): r[f"c{c}"] = e[dxy <= c].sum() / etot
            px = extrapolate_last(X, dv[:, :3], pa, has)
            r["bx"] = e[in3x3(tdep, towers(px[:, 0], px[:, 1], px[:, 2]))].sum() / etot
            dxy = np.hypot(x - px[:, 0], y - px[:, 1])
            for c in (5, 10): r[f"x{c}"] = e[dxy <= c].sum() / etot
            r["beyond"] = e[np.abs(z) > np.abs(X[:, 2]).max()].sum() / etot  # energy deeper than the path reaches
        rows.append(r)
    return rows

def fmt(v, w=5): return f"{v:{w}.2f}" if v is not None else " " * (w - 1) + "-"
def med(rows, key):
    v = [r[key] for r in rows if key in r]; return statistics.median(v) if v else None
def mean(rows, key):
    v = [r[key] for r in rows if key in r]; return statistics.mean(v) if v else None
def ge90(rows, key):
    v = [r[key] for r in rows if key in r]; return sum(1 for x in v if x >= 0.9) / len(v) if v else None

KEYS = ["a", "ar", "av", "b", "c3", "c5", "c10", "bx", "x5", "x10"]
CLASSES = ["reach", "hole", "pi inelastic before face", "stopped before face"]
for fn in args.files:
    t0 = time.time(); rows = analyse(fn); dt = time.time() - t0
    parts = fn.rstrip("/").split("/"); name = parts[-2]; n = len(rows)
    cc = collections.Counter(r["cls"] for r in rows)
    reach = [r for r in rows if r["cls"] == "reach"]
    print(f"=== {parts[-3]}/{name}: {n} events, {dt:.0f} s;  " + ", ".join(f"{c} {cc[c]}" for c in CLASSES if cc[c]))
    if reach:
        print(f"  reach: face r med {med(reach, 'r_face'):.1f} cm (min {min(r['r_face'] for r in reach):.1f}, max {max(r['r_face'] for r in reach):.1f}),"
              f" |eta| face med {med(reach, 'eta_face'):.2f} (vertex eta med {med(reach, 'eta0'):.2f}),"
              f" pT face med {med(reach, 'pt_face'):.3f} GeV = {statistics.median(r['pt_face'] / r['pt0'] for r in reach):.2f} of pT0")
    hole = [r for r in rows if r["cls"] == "hole"]
    if hole:
        print(f"  hole:  face r med {med(hole, 'r_face'):.1f} cm, vertex eta med {med(hole, 'eta0'):.2f}, with E_vis > 0: {sum(1 for r in hole if r['evis'] > 0)}")
    for c in ("pi inelastic before face", "stopped before face"):
        rr = [r for r in rows if r["cls"] == c]
        if not rr: continue
        dc = collections.Counter(r["death"] for r in rr)
        zb = collections.Counter("|z|<100" if r["z_death"] < 100 else ("100-270" if r["z_death"] < 270 else ("270-320" if r["z_death"] < 320 else ">=320")) for r in rr)
        print(f"  {c}: death {dict(dc.most_common(4))}; |z| {dict(zb)}; median |z| {med(rr, 'z_death'):.0f} r {med(rr, 'r_death'):.0f} cm,"
              f" vertex eta med {med(rr, 'eta0'):.2f}, KE at death med {1000 * med(rr, 'ke_death'):.1f} MeV, with E_vis > 0: {sum(1 for r in rr if r['evis'] > 0)}")
    groups = [(c, [r for r in rows if r["cls"] == c and r["evis"] > 0]) for c in CLASSES]
    groups.append(("all with E_vis", [r for r in rows if r["evis"] > 0]))
    print(f"  {'class':26s} {'n':>4s} {'Evis med':>8s} | median: " + " ".join(f"{k:>5s}" for k in KEYS) + " | mean: " + " ".join(f"{k:>5s}" for k in KEYS)
          + " | >=90%: a, ar, b, bx | beyond med/mean")
    for c, rr in groups:
        if not rr: continue
        print(f"  {c:26s} {len(rr):4d} {med(rr, 'evis'):8.1f} | median: " + " ".join(fmt(med(rr, k)) for k in KEYS)
              + " | mean: " + " ".join(fmt(mean(rr, k)) for k in KEYS)
              + " | " + " ".join(fmt(ge90(rr, k), 4) for k in ("a", "ar", "b", "bx")) + f" | {fmt(med(rr, 'beyond'), 4)} {fmt(mean(rr, 'beyond'), 4)}")
    sys.stdout.flush()
