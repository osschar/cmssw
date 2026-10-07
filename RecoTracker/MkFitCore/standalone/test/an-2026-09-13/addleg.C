// Put the legend INSIDE the plot as a TLegend, written into the data file.
// It was an HTML row above the canvas, which ate iframe height and was the
// reason the plots clipped. Colours are set here too, so the page no longer
// has to know them -- the data file becomes self-describing.
#include <TLegend.h>
static void one(const char *fn, std::vector<std::tuple<const char*,int,const char*>> ent,
                double x1, double y1, double x2, double y2, const char *title){
  TFile f(fn,"UPDATE");
  TLegend *L = new TLegend(x1,y1,x2,y2);
  L->SetName("leg"); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.040);
  // A TLegend written to file without ever being PAINTED has fX1NDC..fY2NDC all
  // zero -- the constructor only fills fX1..fY2. Set both, or a reader that
  // trusts GetX1NDC() gets (0,0,0,0) and the legend lands in the corner.
  L->SetX1NDC(x1); L->SetY1NDC(y1); L->SetX2NDC(x2); L->SetY2NDC(y2);
  int n = 0;
  for (auto &e : ent){
    TH1F *h = (TH1F*) f.Get(std::get<0>(e));
    if (!h) { printf("  %s: MISSING %s\n", fn, std::get<0>(e)); continue; }
    h->SetLineColor(std::get<1>(e)); h->SetLineWidth(3);
    h->SetMarkerColor(std::get<1>(e)); h->SetStats(0);
    if (n == 0 && title) h->SetTitle(title);
    L->AddEntry(h, std::get<2>(e), "l");
    h->Write(h->GetName(), TObject::kOverwrite);
    ++n;
  }
  L->Write("leg", TObject::kOverwrite);
  f.Close();
  printf("  %s: legend with %d entries\n", fn, n);
}
void addleg(){
  const char *D = "/foo/matevz/mic-dev/validation-slides/data/";
  one(Form("%s200-bkfit-chi2.root",D),
      {{"chi2_pure",kBlue,"good_frac = 1.0"},{"chi2_mostly",kBlack,"#geq 0.8"},{"chi2_dirty",kRed,"< 0.8"}},
      0.62,0.68,0.89,0.88, "backward-fit per-hit #chi^{2} by seed purity;log_{10} #chi^{2};hits");
  one(Form("%s240-inward-search.root",D),
      {{"chi2_tobps_on",kBlue,"TOB PS 4-9"},{"chi2_pixb_on",kBlack,"PixB 0-3"},{"chi2_disks_on",kRed,"fwd disks"}},
      0.62,0.68,0.89,0.88, "search #chi^{2}, MC-matched hits;log_{10} #chi^{2};hits");
  one(Form("%s250-search-material.root",D),
      {{"chi2_pixb_on",kBlue,"PixB, material ON"},{"chi2_pixb_off",kRed,"PixB, OFF"},
       {"chi2_disks_off",kMagenta+1,"disks, OFF"}},
      0.58,0.68,0.89,0.88, "search #chi^{2}: material ON vs OFF;log_{10} #chi^{2};hits");
  one(Form("%s230-full-chain.root",D),
      {{"ratio",kBlue,"median #chi^{2} / 1.3863"},{"lam",kBlack,"whitened eigenvalues"}},
      0.15,0.72,0.52,0.88, "full chain, all 48 configs;value;count");
  one(Form("%s330-eta.root",D),
      {{"won",kBlue,"WON (best hit)"},{"presel",kRed,"lost to pre-selection"}},
      0.42,0.72,0.89,0.88, "pixel barrel vs |#eta|;|#eta|;fraction");
  one(Form("%s340-preselect-dq.root",D),
      {{"ratio_dq",kRed,"dq / dq-cut"},{"ratio_dphi",kBlue,"d#phi / d#phi-cut"}},
      0.62,0.74,0.89,0.88, "rejected true hits: residual / its own cut;log_{10} ratio;hits");
  one(Form("%s350-q-covariance.root",D),
      {{"track_term",kBlack,"cut's track term (9#sigma_{q})"},
       {"dq_rejected",kRed,"|dq|, rejected true hits"},
       {"dq_passed",kBlue,"|dq|, accepted true hits"}},
      0.35,0.48,0.84,0.68,   // moved: was clipped at the pad edge "pixel barrel: the dq cut's track term vs the residual;|#eta|;cm");
}
