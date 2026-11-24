#!/usr/bin/env python3
"""g4s_mt_check.py A.root -- B.root...: the snitch's output of one job (A, e.g. one thread) against another split over
several files (B, e.g. G4Snitch_t*.root of a multi-threaded job), event by event, matched by branch "ev" (the CMSSW
event number). Per event it compares the counts of particles, stored sensitive and inert steps, inert sums and path
points, the summed deposits of the sensitive and the inert steps and of the inert sums, the summed path length
of the paths (from the points), and the particles by (PDG code, birth process name, birth position rounded to 1 um,
momentum rounded to 1 keV); process names come from each file's own table P, so files numbered differently compare
equal. Prints the events only in one side and every event that differs, then a summary line."""
import sys, math, collections, ROOT
ROOT.gSystem.Load("libRootG4Snitch.so")


def summarise(files):
    out = {}
    for fn in files:
        f = ROOT.TFile.Open(fn)
        P = [str(x.name) for x in f.Get("P")]
        t = f.Get("T")
        if not t.GetBranch("ev"):
            sys.exit(f"{fn}: no branch ev (written before the event number was stored)")
        for ev in t:
            p, s = ev.p, ev.s
            n_sens = n_inert = n_sums = n_path = 0
            e_sens = e_inert = e_sums = plen = 0.0
            parts = collections.Counter()
            for k in range(p.size()):
                q = p[k]
                b = q.m_x_birth
                pm = q.m_p_birth
                proc = P[q.m_birth_proc] if 0 <= q.m_birth_proc < len(P) else ""
                parts[(q.m_pdg, proc, round(b.fX, 4), round(b.fY, 4), round(b.fZ, 4),
                       round(pm.fX, 6), round(pm.fY, 6), round(pm.fZ, 6))] += 1
                st = s[k]
                n_sens += st.m_steps_sensitive.size()
                n_inert += st.m_steps_inert.size()
                n_sums += st.m_inert_sums.size()
                n_path += st.m_path.size()
                for x in st.m_steps_sensitive:
                    e_sens += x.m_edep.m_total
                for x in st.m_steps_inert:
                    e_inert += x.m_edep.m_total
                for x in st.m_inert_sums:
                    e_sums += x.m_edep.m_total
                path = st.m_path
                for i in range(path.size() - 1):
                    a, c = path[i].m_x, path[i + 1].m_x
                    plen += math.sqrt((c.fX - a.fX) ** 2 + (c.fY - a.fY) ** 2 + (c.fZ - a.fZ) ** 2)
            key = int(ev.ev)
            if key in out:
                print(f"event {key} twice ({fn})")
            out[key] = dict(particles=p.size(), sens=n_sens, inert=n_inert, sums=n_sums, path=n_path,
                            e_sens=round(e_sens, 6), e_inert=round(e_inert, 6), e_sums=round(e_sums, 6),
                            plen=round(plen, 3), parts=parts)
        f.Close()
    return out


args = sys.argv[1:]
i = args.index("--")
A, B = summarise(args[:i]), summarise(args[i + 1:])
only_a, only_b = sorted(set(A) - set(B)), sorted(set(B) - set(A))
if only_a:
    print(f"events only in A: {only_a}")
if only_b:
    print(f"events only in B: {only_b}")
ndiff = 0
for e in sorted(set(A) & set(B)):
    a, b = A[e], B[e]
    d = [k for k in a if a[k] != b[k]]
    if d:
        ndiff += 1
        print(f"event {e}: differs in {d}: " + ", ".join(f"{k} {a[k]} / {b[k]}" for k in d if k != "parts"))
nq = len(next(iter(A.values()))) if A else 0
print(f"{len(set(A) & set(B))} events in both ({len(A)} in A, {len(B)} in B), {nq} quantities each; "
      f"events differing: {ndiff}")
