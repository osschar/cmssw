#!/usr/bin/env python3
"""copyno_offsets.py NUMBERING_XML > TABLE: the copy-number offsets that DDG4Builder adds when it builds the Geant4
geometry from DDD: copy = DDD copy + CopyNoOffset + CopyNoTag, SpecPar parameters of the logical parts (the muon
numbering decodes its levels from them). The TGeo export (TGeoMgrFromDdd) keeps the DDD copy number, so a Geant4 path
maps onto it after subtracting these. Output: one line per PartSelector, "regex<TAB>offset", the regex for the
volume name without the DDD namespace. For D128: Geometry/MuonCommonData/data/muonNumbering/TDR_DeV/v5/muonNumbering.xml."""
import re, sys
import xml.etree.ElementTree as ET

root = ET.parse(sys.argv[1]).getroot()
def strip(tag): return tag.split('}')[-1]
consts = {}
for el in root.iter():
    if strip(el.tag) == "Constant":
        consts[el.get("name")] = el.get("value")
def ev(expr, depth=0):
    if depth > 20:
        raise ValueError("constant recursion: " + expr)
    expr = re.sub(r"\[(?:[\w]+:)?(\w+)\]", lambda m: "(" + str(ev(consts[m.group(1)], depth + 1)) + ")", expr)
    if not re.fullmatch(r"[0-9+\-*/() .]+", expr):
        raise ValueError("unexpected expression: " + expr)
    return int(eval(expr))  # digits and arithmetic only, checked above
for sp in root.iter():
    if strip(sp.tag) != "SpecPar":
        continue
    sels = [c.get("path") for c in sp if strip(c.tag) == "PartSelector"]
    pars = {c.get("name"): c.get("value") for c in sp if strip(c.tag) == "Parameter"}
    if "CopyNoOffset" not in pars and "CopyNoTag" not in pars:
        continue
    off = ev(pars.get("CopyNoOffset", "0")) + ev(pars.get("CopyNoTag", "0"))
    for s in sels:
        name = s.split("//")[-1].split(":")[-1]  # //GE0L.* -> GE0L.*
        print(f"{name}\t{off}")
