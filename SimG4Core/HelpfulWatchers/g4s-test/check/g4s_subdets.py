#!/usr/bin/env python3
"""g4s_subdets.py FILE: sensitive steps with a deposit and a DetId, by detector class and DetId subdetector, with
their |eta| range. Tracker subdetectors in phase 2: 1 PixelBarrel, 2 PixelEndcap (inner tracker), 5 TOB, 4 TID (outer
tracker); Ecal: 1 EB, 2 EE; Muon: 1 DT, 2 CSC, 3 RPC, 4 GEM, 5 ME0; HCalSD: 1 HB, 2 HE, 3 HO, 4 HF."""
import sys, math, ROOT, collections
ROOT.gSystem.Load("libRootG4Snitch.so")
f = ROOT.TFile.Open(sys.argv[1]); t = f.Get("T")
vol_sd = [int(v.sd) for v in f.Get("V")]
cls = {int(d.idx): str(getattr(d, "class")) for d in f.Get("D")}
n = collections.Counter(); eta = {}
for ev in t:
    for ps in ev.s:
        for st in ps.m_steps_sensitive:
            if st.m_edep.m_total <= 0 or st.m_det_id == 0: continue
            c = cls.get(vol_sd[st.m_volume], "none")
            # HCalSD gives HCAL test-numbering indices: the top four bits are the HCAL part (1 HB, 2 HE, 3 HO, 4 HF).
            # HGCal DetIds carry no subdetector field there: grouped by class only.
            sub = (st.m_det_id >> 28) if c == "HCalSD" else ("-" if c.startswith("HGC") else (st.m_det_id >> 25) & 7)
            k = (c, sub)
            n[k] += 1
            x = st.m_x_beg; r = math.hypot(x.fX, x.fY); e = abs(math.asinh(x.fZ / r)) if r > 0 else 9
            lo, hi = eta.get(k, (9, 0)); eta[k] = (min(lo, e), max(hi, e))
for k in sorted(n, key=str):
    print(f"  {k[0]:34s} subdet {k[1]}  steps {n[k]:6d}  |eta| {eta[k][0]:.2f}-{eta[k][1]:.2f}")
