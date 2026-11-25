#!/usr/bin/env python3
"""g4s_subtree_check.py pick ORIG EVENT           print the Geant4 id of the tracked non-primary with the largest subtree
   g4s_subtree_check.py ORIG EVENT ID SELECTED    compare a snitch file written with select_track_id=ID (one event)
                                                   with the subtree of ID in entry EVENT of the full file ORIG.
The subtree is laid out as G4Snitch::select_output() writes it: slot 0 the administrative mother, then each particle's
daughters as one contiguous block, breadth first. Compared exactly: the particles, and per particle its steps, inert
sums and path points."""
import sys, ROOT
from SimG4Core.HelpfulWatchers.g4s_read import SnitchFile

def subtree(p, root):
    order = [0, root]; parent = {root: 0}; k = 1
    while k < len(order):
        o = order[k]; db, de = p[o].m_daughters_begin, p[o].m_daughters_end
        if db >= 0:
            for j in range(db, de): parent[j] = k; order.append(j)
        k += 1
    return order, parent

# The process and volume tables are numbered in order of first use within a job, so indices are compared by name.
def tables(sf):
    return [str(x.name) for x in sf.rows("P")], [str(x.path) for x in sf.rows("V")]

def pkey(q, P, V):
    v = lambda a: (a.fX, a.fY, a.fZ, a.fT)
    nm = lambda t, i: t[i] if i >= 0 else None
    return (q.m_pdg, q.m_g4_id, q.m_g4_level, q.m_was_tracked, v(q.m_x_birth), v(q.m_p_birth), v(q.m_x_last),
            v(q.m_p_last), v(q.m_x_death), v(q.m_p_death), nm(P, q.m_birth_proc), nm(P, q.m_death_proc),
            nm(V, q.m_death_volume), q.m_n_daughters_skipped, q.m_n_sub_daughters_skipped, q.n_daughters())

def skey(ps, P, V):
    nm = lambda t, i: t[i] if i >= 0 else None
    st = lambda v: [(x.m_g4_pid, x.m_pdg, nm(P, x.m_proc), nm(V, x.m_volume), x.m_det_id, x.m_edep.m_total, x.m_x_beg.fX,
                     x.m_x_end.fZ, x.m_x0_beg, x.m_lambda_beg) for x in v]
    return (st(ps.m_steps_sensitive), st(ps.m_steps_inert),
            [(nm(V, s.m_volume), s.m_n_steps, s.m_edep.m_total, s.m_x_mean.fZ, tuple(s.m_cov[i] for i in range(6))) for s in ps.m_inert_sums],
            [(x.m_x.fX, x.m_x.fZ, x.m_p.fT, x.m_x0, x.m_lambda) for x in ps.m_path])

if sys.argv[1] == "pick":
    sf = SnitchFile(sys.argv[2]); p = sf.entry(int(sys.argv[3])).p
    nprim = p[0].m_daughters_end - p[0].m_daughters_begin
    best, bid = 0, 0
    for i in range(1, p.size()):
        if p[i].m_was_tracked and p[i].m_g4_id > nprim:
            n = len(subtree(p, i)[0])
            if n > best: best, bid = n, p[i].m_g4_id
    print(bid)
    sys.exit(0)

sfo = SnitchFile(sys.argv[1]); fo = sfo.f; to = sfo.entry(int(sys.argv[2])); po, so = to.p, to.s
gid = int(sys.argv[3])
sfs = SnitchFile(sys.argv[4]); fs = sfs.f; ts = sfs.entry(0); ps_, ss = ts.p, ts.s
Po, Vo = tables(sfo); Ps, Vs = tables(sfs)
root = next(i for i in range(1, po.size()) if po[i].m_was_tracked and po[i].m_g4_id == gid)
order, parent = subtree(po, root)
bad = []
if ps_.size() != len(order): bad.append(f"particles {ps_.size()} against {len(order)}")
for k in range(1, min(len(order), ps_.size())):
    o = order[k]
    if pkey(po[o], Po, Vo) != pkey(ps_[k], Ps, Vs) or ps_[k].m_parent != parent[o]: bad.append(f"particle {k}")
    if skey(so[o], Po, Vo) != skey(ss[k], Ps, Vs): bad.append(f"steps of particle {k}")
n_steps = sum(ss[k].m_steps_sensitive.size() + ss[k].m_steps_inert.size() for k in range(ss.size()))
print(f"subtree of Geant4 id {gid}: {len(order) - 1} particles, {n_steps} steps; differences: {len(bad)}" + (f" {bad[:5]}" if bad else ""))
