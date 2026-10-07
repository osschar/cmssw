#include "tot.C"
void tots2(){
  printf("\nPER-TRACK chi2/(2*n_hits), PURE tracks, expectation ~0.96\n");
  printf("  %-9s %7s %9s %9s %9s %9s | %9s %9s %9s\n","config","tracks","p25","median","p75","p90","all med","all p90","all max");
  tot("val-bkfit-sc1.root",1); tot("val-bkfit-nomat.root",0);
}
