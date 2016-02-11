#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include <rawkfunc.h>
//#include <dynamory.h>

#include <rwk_stats.h>
#include <rwk_popstats.h>

int main(int argc, char **argv) {

  int nsam = 10;
  char factor[] = "chr";
  int winsize = 10000;
  
  // expects stdin in bed format, with the 4th and 5th columns
  // representing the number of ref and alt alleles respectively:
  // chrom    pos-1    pos    nref    nalt

  int ncols = 5;
  int lcols = 1024;
  int lwidth = 2048;
  char delim = '\t';
  char newline = '\n';
  char buffer[lwidth];

  char **array;
  array = calloc(ncols, sizeof (char*));
  
  char chr[lcols];
  long long int startpos;
  long long int stoppos;
  int nref;
  int nalt;
  int coln;
  
  long long int start_region;
  long long int stop_region;

  
  int s,pi;

  double tw_val, pi_val;
  
  struct rwkThetaW thetaW;
  struct rwkThetaPi thetaPi;
  
  rwkThetaWInit(&thetaW, nsam);
  rwkThetaPiInit(&thetaPi, nsam);

  char *tmp;
  int iwin = 1;
  int startindex = 0;
  while (fgets(buffer, sizeof(buffer), stdin)) {
    coln = 0;
    tmp = buffer;
    array[0] = tmp;
    while (*tmp && newline != *tmp) {
      if (delim == *tmp) {
	*tmp = '\0';
	coln++;
	tmp++;
	array[coln] = tmp;

      } else {
	tmp++;
      }
    }
    *tmp = '\0';
    
    startpos = atoll(array[1]);
    stoppos = atoll(array[2]);
    nref = atoi(array[3]);
    nalt = atoi(array[4]);

    if (startindex == 0) {
      strcpy(chr, array[0]);
      start_region = startpos;
      startindex = 1;
    }
    
    if (strcmp(array[0], chr) != 0) {
      printf("%s\t%lld\t%lld\t%d\t%d\t%d", chr, start_region, stop_region, thetaW.nsam, thetaW.nsites, thetaW.s);
      tw_val = thetaW.eval(&thetaW);
      pi_val = thetaPi.eval(&thetaPi);
      printf("\t%lf\t%lf", tw_val / thetaW.nsites, pi_val / thetaPi.nsites);
      printf("\t%lf\n", rwkTajD(thetaW.nsam, thetaW.s, tw_val, pi_val));
      thetaW.reset(&thetaW);
      thetaPi.reset(&thetaPi);
      strcpy(chr, array[0]);
      iwin = 1;
      start_region = 0;
    } else if (iwin >= winsize) {
      printf("%s\t%lld\t%lld\t%d\t%d\t%d", chr, start_region, stop_region, thetaW.nsam, thetaW.nsites, thetaW.s);
      tw_val = thetaW.eval(&thetaW);
      pi_val = thetaPi.eval(&thetaPi);
      printf("\t%lf\t%lf", tw_val / thetaW.nsites, pi_val / thetaPi.nsites);
      printf("\t%lf\n", rwkTajD(thetaW.nsam, thetaW.s, tw_val, pi_val));
      thetaW.reset(&thetaW);
      thetaPi.reset(&thetaPi);
      iwin = 1;
      start_region = startpos;
    }
    
    if (nref > 0 && nalt > 0) {
      s = 1;
      pi = nref;
    } else {
      s = 0;
      pi = 0;
    }

    if (nref + nalt == nsam) {
      thetaW.add(&thetaW, s);
      thetaPi.add(&thetaPi, pi);
    }
    iwin++;
    stop_region = stoppos;
  }
    
  printf("%s\t%lld\t%lld\t%d\t%d\t%d", chr, start_region, stop_region, thetaW.nsam, thetaW.nsites, thetaW.s);
  tw_val = thetaW.eval(&thetaW);
  pi_val = thetaPi.eval(&thetaPi);
  printf("\t%lf\t%lf", tw_val / thetaW.nsites, pi_val / thetaPi.nsites);
  printf("\t%lf\n", rwkTajD(thetaW.nsam, thetaW.s, tw_val, pi_val));

    

  // Free memory
  free(array);
  
  return 0; }
