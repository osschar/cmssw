#include "mkplots.C"
void runplots(){
  const char *D = "/foo/matevz/mic-dev/validation-slides/data/";
  p200("val-bkfit-t5.root",       Form("%s200-bkfit-chi2.root",D));
  p240("val-search-mat1.root","val-search-mat0.root", Form("%s240-inward-search.root",D));
  p230("val-covxport.root",       Form("%s230-full-chain.root",D));
  p330("val-miss-pt5c.root",      Form("%s330-eta.root",D));
}
