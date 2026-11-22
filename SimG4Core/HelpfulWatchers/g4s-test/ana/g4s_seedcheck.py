"""g4s_seedcheck.py FILE...: per event, the number of particles, of primaries, the first primary's vertex z and
momentum x, z, and the summed sensitive deposit (MeV); to check that equal G4SNITCH_SEED give equal events."""
import sys, ROOT
ROOT.gSystem.Load("libRootG4Snitch.so")
def summ(fn):
    f = ROOT.TFile.Open(fn); t = f.Get("T"); out = []
    for ev in t:
        p, s = ev.p, ev.s
        q = p[1]
        esens = sum(ps.m_edep_sum_sensitive.m_total for ps in s)
        out.append((p.size(), p[0].m_daughters_end - p[0].m_daughters_begin, round(q.m_x_birth.fZ, 6), round(q.m_p_birth.fX, 6), round(q.m_p_birth.fZ, 6), round(esens, 3)))
    return out
for fn in sys.argv[1:]:
    print(fn.split("/")[-2], summ(fn))
