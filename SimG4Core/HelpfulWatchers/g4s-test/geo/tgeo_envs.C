// tgeo_envs.C(zpad): print the TGeo material session's envelopes (mb_ladder_layers, mb_discs, mb_envelopes) for
// g4s_tgeo_cmp.py. Run in a scratch directory holding a COPY of /baz/matevz/root-dev/geo-stuff/mat_budget.C and a link
// to its cmsSimGeom-Run4D123.root, so nothing is compiled or written in that session's directory:
//   root -l -b -q -e '.L mat_budget.C+' -e '.x tgeo_envs.C(0.5)' | grep -E '^(L|D|E) ' > envelopes_zpad05.txt
void tgeo_envs(double zpad = 0.5) {
  mb_init();
  mb_disc_zpad = zpad;
  for (auto &l : mb_ladder_layers()) printf("L %s %.4f %.4f %.4f\n", l.name.Data(), l.rmin, l.rmax, l.zhalf);
  for (auto &d : mb_discs("ALL")) printf("D %s %.4f %.4f %.4f %.4f\n", d.name.Data(), d.rmin, d.rmax, d.z, mb_disc_halfz(d));
  for (auto &e : mb_envelopes()) printf("E %.4f %.4f %.4f %.4f\n", e[0], e[1], e[2], e[3]);
}
