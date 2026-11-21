# Mapping the snitch's volumes onto the TGeo export

A G4Snitch file's volume table (tree `V`) names each placement by its Geant4 path, `name:copy` per level from the
world down. The TGeo export of the same geometry (`Fireworks/Geometry`, `TGeoMgrFromDdd`; for D128
`/foo/matevz/cms-dev/g4-snitch/fwgeom/cmsSimGeom-Run4_D128.root`) has the same names (with the DDD namespace) and the
DDD copy numbers, with two differences:
- the export has no `CMSE` level (an identity container): skip it;
- Geant4's copy number is DDD's plus `CopyNoOffset + CopyNoTag`, SpecPars that `DDG4Builder` adds for the muon numbering
  (`MEP` 1 -> 1201, `GE0L` 6 -> 3307). `copyno_offsets.py` reads them from the numbering XML into a table;
  `copyno_offsets_D128.txt` is the D128 one (`Geometry/MuonCommonData/data/muonNumbering/TDR_DeV/v5/muonNumbering.xml`).

`v2tgeo.C+(snitch, tgeo, offsets[, nmax])` (ROOT, ACLiC) maps every volume and checks that the mapped node's global
position equals the volume's origin as Geant4 saw it (`V.x/y/z`). D128, `run20cms` (gun into most of CMS): 120209 of
120210 volumes mapped (the one left is `CMSE` itself), 20992 path levels with an offset, largest position difference
2.3e-7 cm. ROOT gives each physical node a fixed id after `gGeoManager->GetCache()->BuildIdArray()` (or
`CloseGeometry("i")`): `cd(path)` then `GetCurrentNodeId()`, and `CdNode(id)` back.
