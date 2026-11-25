#!/usr/bin/env python3
"""g4s_eloss.py GEO.root OUT.root FILE...: mkFit's energy loss against Geant4's, per step between consecutive tracker
hits of the primary (single gun; hits and the outgoing leg as in g4s_tkmat.py).

Truth: the primary's energy leaving each hit (the end of the hit's last step; the step before a hit is often not
stored, as thin air gaps deposit less than the 10 eV storage cut, so the entry energy is not used); the loss of a
step is the difference between consecutive hits: the material after one sensor up to the end of the next, the
delta rays the primary emits included.
Model: mkFit's applyMaterialEffects() on arriving at the destination hit, transcribed: the mean Bethe-Bloch loss
  dE = 2 xi p/|p.n| (0.5 ln(2 me wmax / I^2) - beta^2 - deltahalf) / beta^2,
  wmax = 2 me beta^2 gamma^2 / (1 + 2 gamma me / mpi + (me / mpi)^2), I = 16 eV * 10.75,
  deltahalf = ln(28.816 eV sqrt(2.33 * 0.498) / I) - 0.5, beta and gamma with the PION mass (mkFit's), me = 0.5 MeV,
with xi = the destination module's bbxi (MediumProperties, xi_mod) or the (|z|, r) grid's (xi_grid), and p the true
momentum leaving the destination hit. dp = dE / beta is what mkFit subtracts from p; here energies are compared.
Tree "S": eta, pt, pdg, lay0, lay1, r1, z1, p1 (GeV), de (truth, GeV), de_mod, de_grid.
Tree "K": eta, pt, pdg, nhits, p_first, de, de_mod, de_grid (sums between the first and the last hit)."""
import sys, math, array, collections, ROOT
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
ROOT.gSystem.Load("libRootG4Snitch.so")

ME, MPI = 0.0005, 0.140
I = 16.0e-9 * 10.75
DELTAHALF = math.log(28.816e-9 * math.sqrt(2.33 * 0.498) / I) - 0.5


def mkfit_de(xi, invcos, p):
    p2 = p * p
    beta2 = p2 / (p2 + MPI * MPI)
    beta = math.sqrt(beta2)
    gamma2 = (p2 + MPI * MPI) / (MPI * MPI)
    gamma = math.sqrt(gamma2)
    wmax = 2 * ME * beta2 * gamma2 / (1 + 2 * gamma * ME / MPI + ME * ME / (MPI * MPI))
    if beta >= 1:
        return 0.0
    return 2 * (xi * invcos * (0.5 * math.log(2 * ME * wmax / (I * I)) - beta2 - DELTAHALF) / beta2)


def main():
    geo = ROOT.TFile.Open(sys.argv[1])
    mods = {int(m.detid): (int(m.layer), (m.nx, m.ny, m.nz), m.bbxi) for m in geo.Get("M")}
    hg = geo.Get("bbxi")
    nzb, nrb = hg.GetNbinsX(), hg.GetNbinsY()
    grid = [[hg.GetBinContent(i + 1, j + 1) for j in range(nrb)] for i in range(nzb)]

    def grid_xi(z, r):
        i, j = int(abs(z)), int(r)
        return grid[i][j] if 0 <= i < nzb and 0 <= j < nrb else 0.0

    out = ROOT.TFile(sys.argv[2], "RECREATE")
    scols = ["eta", "pt", "pdg", "lay0", "lay1", "r1", "z1", "p1", "de", "de_mod", "de_grid"]
    kcols = ["eta", "pt", "pdg", "nhits", "p_first", "de", "de_mod", "de_grid"]
    S, Kt = ROOT.TTree("S", "steps"), ROOT.TTree("K", "tracks")
    sb = {c: array.array("f", [0.0]) for c in scols}
    kb = {c: array.array("f", [0.0]) for c in kcols}
    for c in scols:
        S.Branch(c, sb[c], c + "/F")
    for c in kcols:
        Kt.Branch(c, kb[c], c + "/F")
    nbreak = collections.Counter()
    for fn in sys.argv[3:]:
        sf = SnitchFile(fn); f = sf.f
        sdcls = [str(getattr(d, "class")) for d in sf.rows("D")]
        trk_sd = {i for i, c in enumerate(sdcls) if "TkAccumulating" in c}
        V = [int(v.sd) in trk_sd for v in sf.rows("V")]
        for ev in sf:
            p, s = ev.p, ev.s
            k = p[0].m_daughters_begin
            q = p[k]
            if q.m_charge == 0:
                continue
            m = q.m_mass
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
            hits = []  # [detid, E_out (the end of the hit's last step), mid, dir]
            prev = None
            for st in own:
                is_hit = st.m_det_id != 0 and 0 <= st.m_volume < len(V) and V[st.m_volume]
                if is_hit:
                    if hits and prev is not None and prev.m_det_id == st.m_det_id and hits[-1][0] == st.m_det_id:
                        hits[-1][1] = st.m_p_end.fT
                    else:
                        mid = (0.5 * (st.m_x_beg.fX + st.m_x_end.fX), 0.5 * (st.m_x_beg.fY + st.m_x_end.fY),
                               0.5 * (st.m_x_beg.fZ + st.m_x_end.fZ))
                        d = (st.m_x_end.fX - st.m_x_beg.fX, st.m_x_end.fY - st.m_x_beg.fY,
                             st.m_x_end.fZ - st.m_x_beg.fZ)
                        hits.append([st.m_det_id, st.m_p_end.fT, mid, d])
                prev = st
            if len(hits) < 2:
                continue
            pb = q.m_p_birth
            pt = math.hypot(pb.fX, pb.fY)
            eta = math.asinh(pb.fZ / pt)
            rows = []
            for h in hits:
                em = h[1]
                pm = math.sqrt(max(em * em - m * m, 0.0))
                mo = mods.get(int(h[0]))
                if mo is None:
                    nbreak["DetId not in mkFit geometry"] += 1
                    rows.append((-1, h[2], em, pm, 0.0, 0.0))
                    continue
                lay, n_, xi = mo
                dn = math.sqrt(sum(x * x for x in h[3]))
                cosv = abs(sum(h[3][c] * n_[c] for c in range(3))) / dn if dn > 0 else 1.0
                ic = 1.0 / max(cosv, 1e-3)
                r = math.hypot(h[2][0], h[2][1])
                rows.append((lay, h[2], em, pm, mkfit_de(xi, ic, pm), mkfit_de(grid_xi(h[2][2], r), ic, pm)))
            tot = collections.Counter()
            for a, b in zip(rows, rows[1:]):
                vals = dict(eta=eta, pt=pt, pdg=q.m_pdg, lay0=a[0], lay1=b[0], r1=math.hypot(b[1][0], b[1][1]),
                            z1=b[1][2], p1=b[3], de=a[2] - b[2], de_mod=b[4], de_grid=b[5])
                for c in scols:
                    sb[c][0] = vals[c]
                S.Fill()
                for c in ("de", "de_mod", "de_grid"):
                    tot[c] += vals[c]
            kv = dict(eta=eta, pt=pt, pdg=q.m_pdg, nhits=len(rows), p_first=rows[0][3], **tot)
            for c in kcols:
                kb[c][0] = kv[c]
            Kt.Fill()
        f.Close()
    print(f"{Kt.GetEntries()} tracks, {S.GetEntries()} steps; {dict(nbreak)}")
    out.cd()
    S.Write()
    Kt.Write()
    out.Close()


main()
