#!/usr/bin/env python3
"""g4s_minbias.py FILE...: what one inelastic pp interaction (G4SNITCH_GUN=minbias) puts into HGCal, by the particle
that brought it in, with Geant4 truth. For the charged-track follower: how much of the visible HGCal energy comes in
on charged particles of low pT, and how much of it the towers around the particle's own crossing of the HGCal face
would miss.

Visible energy: sensitive steps with an HGCal DetId (detector = id >> 28: 8 CE-E silicon, 9 CE-H silicon, 10 CE-H
scintillator) and a deposit above zero, at the step midpoint, in MeV. No threshold, noise or digitisation.
Entrant of a deposit: the nearest ancestor-or-self of the particle the step is filed under (a tracked particle; the
steps of untracked descendants, stable particles created with less than 10 MeV of kinetic energy, are filed under
their nearest tracked ancestor) whose birth point has |z| < 320 cm.
Classes: charged entrants by pT at birth (< 0.5, 0.5-1, 1-2, 2-5, > 5 GeV), photons, neutral hadrons (neutral,
100 < |pdg| < 1e9: n, K0L, K0S, Lambda, ...), other. Origin of an entrant: a Geant4 primary (from the generator);
else by its birth point: in the tracker volume (|z| < 270, r < 120 cm) by decay (birth process Decay) or by an
interaction in the material (any other); in front of HGCal (270 <= |z| < 320 cm: the tracker endcap's services, ETL,
HGCal's neutron moderator and thermal screen; HGCal's first CE-E cell is at |z| = 322.2 cm in D128); or elsewhere
(|z| < 270, r >= 120 cm). "Track-like": charged, a primary or a decay product in the tracker volume.
Towers: the hadronic grid of the pointing-tower demonstrator, |eta| from 1.305 in steps of 0.087, 72 in phi from -pi,
pointing at the origin, one grid per endcap. A charged entrant's tower is where its own path (m_path, linear between
points) first crosses |z| = 320 cm; its visible energy is split into that tower's 3x3 (and 1x1, 5x5) and the rest.
Printed: means per interaction over the events with their standard errors, and times 200 for PU200."""
import sys, math, collections, ROOT
import numpy as np
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile
ROOT.gSystem.Load("libRootG4Snitch.so")

ROOT.gInterpreter.Declare(r'''
#include "SimG4Core/HelpfulWatchers/interface/G4SnitchReader.h"
#include <cmath>
#include <vector>
namespace g4smb {
constexpr double ETA0 = 1.305, DETA = (3.045 - 1.305) / 100 * 5, ZFACE = 320.0;
constexpr double T_PROMPT = 25.0;  // ns, global time of the step midpoint
constexpr int NPHI = 72;
struct Tower { int side, ieta, iphi; };
inline Tower tower_of(double x, double y, double z) {
  const double r = std::hypot(x, y);
  const double eta = r > 0 ? std::asinh(z / r) : (z >= 0 ? 99.0 : -99.0);
  const double phi = std::atan2(y, x);
  Tower t;
  t.side = z >= 0 ? 1 : -1;
  t.ieta = (int)std::floor((std::abs(eta) - ETA0) / DETA);
  t.iphi = (((int)std::floor((phi + M_PI) / (2 * M_PI / NPHI))) % NPHI + NPHI) % NPHI;
  return t;
}
// Distance in towers (the larger of |d ieta| and |d iphi|, phi wrapped); 1000 for the other endcap.
inline int tdist(const Tower &a, const Tower &b) {
  if (a.side != b.side) return 1000;
  int dp = std::abs(a.iphi - b.iphi) % NPHI;
  dp = std::min(dp, NPHI - dp);
  return std::max(std::abs(a.ieta - b.ieta), dp);
}
// Flat per-entrant output over all events of a file, plus per-event totals.
struct Out {
  // per entrant with visible HGCal energy
  std::vector<int> ev, pdg, charge, birth_proc, is_primary, has_cross;
  std::vector<double> pt, eta_birth, z_birth, r_birth, e_vis, e_ee, e_hsi, e_hsc, e_own, e_1, e_3, e_5, e_3ee, e_pt, e_3pt, cross_r, cross_eta;
  // per event
  std::vector<double> ev_e_total, ev_e_unattributed, ev_z_vtx;
  std::vector<int> ev_n_part;
};
void process_file(const char *fn, Out &o, long max_events = -1) {
  G4SnitchFile sf(fn);  // a TTree or an RNTuple file; reads every field of T (no pruning)
  const long n_ev = max_events >= 0 ? std::min(max_events, sf.n_events()) : sf.n_events();
  const int ev_off = o.ev_e_total.size();
  for (long ev = 0; ev < n_ev; ++ev) {
    sf.load(ev);
    const std::vector<G4S_Particle> *p = sf.p;
    const std::vector<G4S_ParticleSteps> *s = sf.s;
    const int n = p->size();
    // The entrant of each particle: the nearest ancestor-or-self with |z_birth| < ZFACE; -1 if the chain reaches the
    // administrative mother (slot 0) first.
    std::vector<int> ent(n, -2);
    auto entrant = [&](int i) {
      std::vector<int> chain;
      int k = i, e = -1;
      while (true) {
        if (k <= 0) { e = -1; break; }
        if (ent[k] != -2) { e = ent[k]; break; }
        chain.push_back(k);
        if (std::abs((*p)[k].m_x_birth.fZ) < ZFACE) { e = k; break; }
        k = (*p)[k].m_parent;
      }
      for (int c : chain) ent[c] = e;
      return e;
    };
    // Accumulate per entrant.
    std::vector<double> E(n, 0), Eee(n, 0), Ehsi(n, 0), Ehsc(n, 0), Eown(n, 0), E1(n, 0), E3(n, 0), E5(n, 0), E3ee(n, 0), Ept(n, 0), E3pt(n, 0);
    std::vector<int> cross_state(n, -1);  // -1 not computed, 0 no crossing, 1 crossing
    std::vector<Tower> ctow(n);
    std::vector<double> cr(n, 0), ceta(n, 0);
    auto crossing = [&](int e) {
      if (cross_state[e] >= 0) return;
      cross_state[e] = 0;
      const auto &path = (*s)[e].m_path;
      for (size_t k = 1; k < path.size(); ++k) {
        const auto &a = path[k - 1].m_x, &b = path[k].m_x;
        if (std::abs(a.fZ) < ZFACE && std::abs(b.fZ) >= ZFACE) {
          const double zt = b.fZ > 0 ? ZFACE : -ZFACE;
          const double w = (zt - a.fZ) / (b.fZ - a.fZ);
          const double x = a.fX + w * (b.fX - a.fX), y = a.fY + w * (b.fY - a.fY);
          ctow[e] = tower_of(x, y, zt);
          cr[e] = std::hypot(x, y);
          ceta[e] = std::asinh(zt / cr[e]);
          cross_state[e] = 1;
          return;
        }
      }
    };
    double e_tot = 0, e_unatt = 0;
    std::vector<int> touched;
    for (int i = 1; i < n; ++i) {
      const auto &ps = (*s)[i];
      for (const G4S_Step &st : ps.m_steps_sensitive) {
        const unsigned det = st.m_det_id >> 28;
        if ((det != 8 && det != 9 && det != 10) || st.m_edep.m_total <= 0) continue;
        const double de = st.m_edep.m_total;
        e_tot += de;
        const int e = entrant(i);
        if (e < 0) { e_unatt += de; continue; }
        if (E[e] == 0) touched.push_back(e);
        E[e] += de;
        (det == 8 ? Eee[e] : det == 9 ? Ehsi[e] : Ehsc[e]) += de;
        const bool prompt = 0.5 * (st.m_x_beg.fT + st.m_x_end.fT) < T_PROMPT;
        if (prompt) Ept[e] += de;
        if (i == e && st.m_g4_pid == (*p)[e].m_g4_id) Eown[e] += de;
        if ((*p)[e].m_charge != 0) {
          crossing(e);
          if (cross_state[e] == 1) {
            const Tower tw = tower_of(0.5 * (st.m_x_beg.fX + st.m_x_end.fX), 0.5 * (st.m_x_beg.fY + st.m_x_end.fY),
                                      0.5 * (st.m_x_beg.fZ + st.m_x_end.fZ));
            const int d = tdist(tw, ctow[e]);
            if (d == 0) E1[e] += de;
            if (d <= 1) E3[e] += de;
            if (d <= 1 && det == 8) E3ee[e] += de;
            if (d <= 1 && prompt) E3pt[e] += de;
            if (d <= 2) E5[e] += de;
          }
        }
      }
    }
    for (int e : touched) {
      const G4S_Particle &q = (*p)[e];
      o.ev.push_back(ev_off + ev);
      o.pdg.push_back(q.m_pdg);
      o.charge.push_back(q.m_charge);
      o.birth_proc.push_back(q.m_birth_proc);
      o.is_primary.push_back(q.m_parent == 0);
      o.has_cross.push_back(q.m_charge != 0 ? cross_state[e] : -1);
      o.pt.push_back(std::hypot(q.m_p_birth.fX, q.m_p_birth.fY));
      const double pz = q.m_p_birth.fZ, ptt = std::hypot(q.m_p_birth.fX, q.m_p_birth.fY);
      o.eta_birth.push_back(ptt > 0 ? std::asinh(pz / ptt) : (pz >= 0 ? 99.0 : -99.0));
      o.z_birth.push_back(q.m_x_birth.fZ);
      o.r_birth.push_back(std::hypot(q.m_x_birth.fX, q.m_x_birth.fY));
      o.e_vis.push_back(E[e]);
      o.e_ee.push_back(Eee[e]);
      o.e_hsi.push_back(Ehsi[e]);
      o.e_hsc.push_back(Ehsc[e]);
      o.e_own.push_back(Eown[e]);
      o.e_1.push_back(E1[e]);
      o.e_3.push_back(E3[e]);
      o.e_3ee.push_back(E3ee[e]);
      o.e_pt.push_back(Ept[e]);
      o.e_3pt.push_back(E3pt[e]);
      o.e_5.push_back(E5[e]);
      o.cross_r.push_back(cr[e]);
      o.cross_eta.push_back(ceta[e]);
    }
    o.ev_e_total.push_back(e_tot);
    o.ev_e_unattributed.push_back(e_unatt);
    o.ev_n_part.push_back(n);
    // The vertex: the birth point of the first primary.
    o.ev_z_vtx.push_back(n > 1 ? (*p)[1].m_x_birth.fZ : 0.0);
  }
}
}
''')

PT_EDGES = [0.5, 1.0, 2.0, 5.0]
CH_CLASSES = ["ch pT<0.5", "ch 0.5-1", "ch 1-2", "ch 2-5", "ch pT>5"]
CLASSES = CH_CLASSES + ["photon", "neutral hadron", "other"]
# Origin of an entrant: a Geant4 primary; else by its birth point: in the tracker volume (|z| < 270, r < 120 cm) by
# decay or by an interaction in the material; in front of HGCal (270 <= |z| < 320 cm: tracker endcap services, ETL,
# HGCal's neutron moderator and thermal screen); or elsewhere (|z| < 270, r >= 120 cm: the barrel calorimeters).
ORIGINS = ["primary", "decay trk", "matter trk", "front", "barrel"]

def main(files):
    o = ROOT.g4smb.Out()
    procs = None
    for fn in files:
        ROOT.g4smb.process_file(fn, o)
        if procs is None:
            sf = SnitchFile(fn); procs = {int(x.idx): str(x.name) for x in sf.rows("P")}; sf.f.Close()
    A = {k: np.array(list(getattr(o, k))) for k in ("ev", "pdg", "charge", "birth_proc", "is_primary", "has_cross", "pt",
                                                     "eta_birth", "z_birth", "r_birth", "e_vis", "e_ee", "e_hsi", "e_hsc",
                                                     "e_own", "e_1", "e_3", "e_5", "e_3ee", "e_pt", "e_3pt", "cross_r", "cross_eta")}
    ev_tot = np.array(list(o.ev_e_total)); ev_un = np.array(list(o.ev_e_unattributed)); nev = len(ev_tot)
    zv = np.array(list(o.ev_z_vtx))
    ne = len(A["pdg"]); one = np.ones(ne); allm = np.ones(ne, bool)
    ch = A["charge"] != 0
    cls = np.full(ne, CLASSES.index("other"))
    cls[ch] = np.digitize(A["pt"][ch], PT_EDGES)
    cls[(~ch) & (A["pdg"] == 22)] = CLASSES.index("photon")
    apdg = np.abs(A["pdg"])
    cls[(~ch) & (apdg > 100) & (apdg < 1000000000)] = CLASSES.index("neutral hadron")
    is_decay = np.array([("Decay" in procs.get(int(b), "")) for b in A["birth_proc"]], dtype=bool)
    zb = np.abs(A["z_birth"]); rb = A["r_birth"]
    origin = np.where(A["is_primary"] == 1, 0, np.where(zb >= 270, 3, np.where(rb >= 120, 4, np.where(is_decay, 1, 2))))
    trk = ch & (origin <= 1)  # charged, from the generator or a decay in the tracker volume: track-like
    xs = A["has_cross"] == 1
    e_out = np.where(xs, A["e_vis"] - A["e_3"], 0.0)
    grid = xs & (np.abs(A["cross_eta"]) >= 1.305) & (np.abs(A["cross_eta"]) < 3.045)  # crossing on the tower grid

    def per_event(mask, w):
        return np.bincount(A["ev"][mask], weights=w[mask], minlength=nev)
    def ms(x):
        return x.mean(), (x.std(ddof=1) / math.sqrt(len(x)) if len(x) > 1 else 0.0)
    tot_m, tot_e = ms(ev_tot)

    print(f"{len(files)} files, {nev} events, one inelastic pp interaction each (14 TeV); vertex z mean {zv.mean():.2f} cm, rms {zv.std():.2f} cm")
    print(f"Visible HGCal energy per interaction: {tot_m:.2f} +- {tot_e:.2f} MeV (x200: {200 * tot_m:.0f} MeV); median {np.median(ev_tot):.2f} MeV;"
          f" events with none: {np.sum(ev_tot == 0)}; not attributed to an entrant: {ev_un.sum():.3f} MeV in all")
    for name, k in (("CE-E Si", "e_ee"), ("CE-H Si", "e_hsi"), ("CE-H scint", "e_hsc")):
        m, e = ms(per_event(allm, A[k]))
        print(f"   {name:10s} {m:8.2f} +- {e:6.2f} MeV  ({100 * m / tot_m:.1f}%)")
    print(f"   at global time t < 25 ns: {ms(per_event(allm, A['e_pt']))[0] / tot_m:.4f} of it")
    print("Means per interaction +- standard error over events; x200 = PU200 (in-time only). N: entrants with visible energy > 0.\n")

    def class_table(title, rows):
        print(title)
        print(f"{'':22s} {'N':>15s} {'N x200':>9s} {'E_vis (MeV)':>17s} {'E x200':>9s} {'share':>7s} {'t<25ns':>7s}")
        for name, sel in rows:
            nm, nerr = ms(per_event(sel, one)); em, eerr = ms(per_event(sel, A["e_vis"]))
            ept = A["e_pt"][sel].sum() / nev
            print(f"{name:22s} {nm:8.3f} +-{nerr:6.3f} {200 * nm:9.1f} {em:8.2f} +-{eerr:6.2f} {200 * em:9.0f} {100 * em / tot_m:6.1f}% {ept / em if em > 0 else 0:7.3f}")
        print()
    class_table("1. Entrants by class (charged by pT at birth, GeV):",
                [(n, cls == c) for c, n in enumerate(CLASSES)]
                + [("charged, all", ch), ("charged pT<2", ch & (A["pt"] < 2)), ("charged pT<5", ch & (A["pt"] < 5))])
    class_table("2. Track-like charged entrants only (Geant4 primaries, and decays in the tracker volume):",
                [(n, trk & (cls == c)) for c, n in enumerate(CH_CLASSES)]
                + [("track-like, all", trk), ("track-like pT<2", trk & (A["pt"] < 2)), ("track-like pT<5", trk & (A["pt"] < 5))])

    print("3. Visible energy per interaction (MeV) [entrants per interaction], by class and origin:")
    print(f"{'':16s} " + " ".join(f"{x:>18s}" for x in ORIGINS))
    for c, name in enumerate(CLASSES):
        cells = []
        for oi in range(len(ORIGINS)):
            sel = (cls == c) & (origin == oi)
            cells.append(f"{A['e_vis'][sel].sum() / nev:8.2f} [{sel.sum() / nev:6.2f}]")
        print(f"{name:16s} " + " ".join(f"{x:>18s}" for x in cells))
    print()

    def contain_table(title, rows):
        print(title)
        print(f"{'':22s} {'E_vis':>8s} {'in 1x1':>7s} {'in 3x3':>7s} {'in 5x5':>7s} {'out of 3x3 (MeV)':>17s} {'x200':>7s} {'CE-E: all':>9s} {'out':>6s} {'out t<25ns':>10s} {'no crossing':>11s}")
        for name, sel in rows:
            E = A["e_vis"][sel].sum() / nev
            if E <= 0:
                print(f"{name:22s} {0:8.2f}"); continue
            e1 = A["e_1"][sel & xs].sum() / nev; e3 = A["e_3"][sel & xs].sum() / nev; e5 = A["e_5"][sel & xs].sum() / nev
            om, oe = ms(per_event(sel, e_out)); enox = A["e_vis"][sel & (A["has_cross"] == 0)].sum() / nev
            oee = (A["e_ee"] - A["e_3ee"])[sel & xs].sum() / nev
            opt = (A["e_pt"] - A["e_3pt"])[sel & xs].sum() / nev
            print(f"{name:22s} {E:8.2f} {e1 / E:7.3f} {e3 / E:7.3f} {e5 / E:7.3f} {om:8.2f} +-{oe:6.2f} {200 * om:7.0f} {A["e_ee"][sel].sum() / nev / E:9.3f} {oee / om if om > 0 else 0:6.3f} {opt / om if om > 0 else 0:10.3f} {enox / E:11.4f}")
        print()
    contain_table("4. Charged entrants: visible energy in the hadronic towers (0.087 x 5 deg, pointing at the origin) around the"
                  " entrant's own crossing of |z| = 320 cm\n   (fractions of the row's visible energy), the energy outside its 3x3 (MeV per interaction), and the CE-E share of all and of the outside:",
                  [(n, ch & (cls == c)) for c, n in enumerate(CH_CLASSES)]
                  + [("charged pT<2", ch & (A["pt"] < 2)), ("charged pT<5", ch & (A["pt"] < 5)), ("charged, all", ch),
                     ("pT<2, on the grid", ch & (A["pt"] < 2) & grid), ("pT<5, on the grid", ch & (A["pt"] < 5) & grid), ("all, on the grid", ch & grid)])
    contain_table("5. The same for track-like charged entrants (Geant4 primaries, decays in the tracker volume):",
                  [(n, trk & (cls == c)) for c, n in enumerate(CH_CLASSES)]
                  + [("track-like pT<2", trk & (A["pt"] < 2)), ("track-like pT<5", trk & (A["pt"] < 5)), ("track-like, all", trk),
                     ("pT<2, on the grid", trk & (A["pt"] < 2) & grid), ("pT<5, on the grid", trk & (A["pt"] < 5) & grid), ("all, on the grid", trk & grid)])
    print("   (on the grid: the crossing of |z| = 320 cm at 1.305 <= |eta| < 3.045, where the 3x3 lies on the tower grid)\n")

    print("6. Charged entrants by |eta| of the crossing of |z| = 320 cm: visible energy and energy outside the pointed 3x3 (MeV per interaction):")
    eb = [1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 99]
    ae = np.abs(A["cross_eta"])
    print(f"{'|eta| at crossing':18s} " + " ".join(f"{'%.1f-%s' % (eb[i], ('%.1f' % eb[i + 1]) if eb[i + 1] < 90 else 'inf'):>14s}" for i in range(len(eb) - 1)) + f" {'<1.3':>14s}")
    for name, sel in (("charged, all", ch), ("charged pT<2", ch & (A["pt"] < 2)), ("track-like pT<2", trk & (A["pt"] < 2))):
        cells = []
        for i in range(len(eb) - 1):
            m = sel & xs & (ae >= eb[i]) & (ae < eb[i + 1])
            cells.append(f"{A['e_vis'][m].sum() / nev:6.2f} ({e_out[m].sum() / nev:5.2f})")
        m = sel & xs & (ae < eb[0]); cells.append(f"{A['e_vis'][m].sum() / nev:6.2f} ({e_out[m].sum() / nev:5.2f})")
        print(f"{name:18s} " + " ".join(f"{x:>14s}" for x in cells))
    print("   (E_vis, in brackets the part outside the pointed 3x3)\n")

    print("7. Charged entrants by species and origin: visible energy per interaction (MeV) (fraction outside the pointed 3x3), and entrants per interaction:")
    species = [("e+-", apdg == 11), ("mu+-", apdg == 13), ("pi+-", apdg == 211), ("K+-", apdg == 321), ("p, pbar", apdg == 2212),
               ("other charged", ~np.isin(apdg, [11, 13, 211, 321, 2212]))]
    print(f"{'':14s} " + " ".join(f"{x:>17s}" for x in ORIGINS) + f" {'N':>8s}")
    for name, sp in species:
        sp = sp & ch
        cells = []
        for oi in range(len(ORIGINS)):
            sel = sp & (origin == oi); E = A["e_vis"][sel].sum() / nev; eo = e_out[sel].sum() / nev
            cells.append(f"{E:7.2f} ({eo / E if E > 0 else 0:5.3f})")
        print(f"{name:14s} " + " ".join(f"{x:>17s}" for x in cells) + f" {sp.sum() / nev:8.3f}")
    print()
    for name, sel in (("charged", ch), ("photon", cls == CLASSES.index("photon")), ("neutral hadron", cls == CLASSES.index("neutral hadron"))):
        E = A["e_vis"][sel].sum()
        if E > 0:
            print(f"{name:15s}: own steps of the entrant make {100 * A['e_own'][sel].sum() / E:5.1f}% of its visible energy, its descendants the rest")

if __name__ == "__main__":
    main(sys.argv[1:])
