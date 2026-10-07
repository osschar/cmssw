#include "tot.C"
void tots(){
  printf("\nPER-TRACK chi2/(2*n_hits), expectation ~0.96 at 15 d.o.f.\n");
  printf("                   |------------- PURE tracks -------------| |------ all tracks -----|\n");
  printf("  %-9s %7s %9s %9s %9s %9s | %9s %9s %9s\n","var scale","tracks","p25","median","p75","p90","median","p90","max");
  tot("val-bkfit-sc1.root",1); tot("val-bkfit-sc10.root",10);
  tot("val-bkfit-sc100.root",100); tot("val-bkfit-sc10000.root",10000);
}
