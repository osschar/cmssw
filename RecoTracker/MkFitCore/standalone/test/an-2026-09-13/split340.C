// One data file per slide, never shared -- so the dq/dphi ratios get their own.
void split340(){
  const char *D="/foo/matevz/mic-dev/validation-slides/data/";
  TFile in(Form("%s330-eta.root",D));
  TH1F *rq=(TH1F*)in.Get("ratio_dq"), *rp=(TH1F*)in.Get("ratio_dphi");
  rq->SetDirectory(nullptr); rp->SetDirectory(nullptr);
  rp->SetTitle("rejected true hits: residual / its own cut;log_{10} ratio;hits");
  { TFile o(Form("%s340-preselect-dq.root",D),"RECREATE"); rq->Write(); rp->Write(); }
  // and drop them from 330 so each slide owns exactly its own numbers
  TH1F *w=(TH1F*)in.Get("won"), *p=(TH1F*)in.Get("presel");
  w->SetDirectory(nullptr); p->SetDirectory(nullptr);
  { TFile o(Form("%s330-eta.root",D),"RECREATE"); w->Write(); p->Write(); }
  printf("split: 330-eta.root = won/presel ; 340-preselect-dq.root = ratio_dq/ratio_dphi\n");
}
