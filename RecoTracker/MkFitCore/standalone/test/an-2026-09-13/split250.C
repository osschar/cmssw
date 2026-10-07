// One data file per slide: 240 keeps material-ON (its subject is the region
// split), 250 keeps the ON/OFF pair (its subject is the material lever).
void split250(){
  const char *D="/foo/matevz/mic-dev/validation-slides/data/";
  TFile in(Form("%s240-inward-search.root",D));
  const char *keep240[3]={"chi2_pixb_on","chi2_tobps_on","chi2_disks_on"};
  const char *keep250[6]={"chi2_pixb_on","chi2_pixb_off","chi2_tobps_on",
                          "chi2_tobps_off","chi2_disks_on","chi2_disks_off"};
  std::vector<TH1F*> a,b;
  for(auto n:keep240){ TH1F*h=(TH1F*)in.Get(n); if(h){h=(TH1F*)h->Clone(); h->SetDirectory(nullptr); a.push_back(h);} }
  for(auto n:keep250){ TH1F*h=(TH1F*)in.Get(n); if(h){h=(TH1F*)h->Clone(); h->SetDirectory(nullptr); b.push_back(h);} }
  { TFile o(Form("%s250-search-material.root",D),"RECREATE"); for(auto h:b) h->Write(); }
  { TFile o(Form("%s240-inward-search.root",D),"RECREATE"); for(auto h:a) h->Write(); }
  printf("240 keeps %zu (material ON only); 250 keeps %zu (ON/OFF pair)\n",a.size(),b.size());
}
