#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_stats.h>
#include <rwk_parse.h>
#include <rwk_popstats.h>

int main(int argc, char **argv) {

  //... | popstats nsam fcol [-w <lwin>]
  
  int nsam = atoi(argv[1]);
  int fcol = atoi(argv[2]) - 1;
  
  //char factor[1000];
  int winsize = 10000;
  
  // expects stdin in bed format, with the 4th and 5th columns
  // representing the number of ref and alt alleles respectively:
  // chrom    pos-1    pos    nref    nalt

  int ncols = 6;
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
  
  long long int start_region;
  long long int stop_region;

  
  int s,pi;

  double tw_val, pi_val;
  
  struct rwkThetaW thetaW;
  struct rwkThetaPi thetaPi;
  
  rwkThetaWInit(&thetaW, nsam);
  rwkThetaPiInit(&thetaPi, nsam);
  
  int coln;
  char *tmp;
  int iwin = 1;
  int startindex = 0;
  while (fgets(buffer, sizeof(buffer), stdin)) {
    rwkStrtoArray(array, buffer, &delim);
    
    startpos = atoll(array[1]);
    stoppos = atoll(array[2]);
    nref = atoi(array[4]);
    nalt = atoi(array[5]);

    if (startindex == 0) {
      strcpy(chr, array[fcol]);
      start_region = startpos;
      startindex = 1;
    }
    
    if (strcmp(array[fcol], chr) != 0) {
      tw_val = thetaW.eval(&thetaW);
      pi_val = thetaPi.eval(&thetaPi);
      printf("%s\t%lld\t%lld\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",
	     chr, start_region, stop_region,
	     thetaW.nsam, thetaW.nsites, thetaW.s,
	     tw_val / thetaW.nsites, pi_val / thetaPi.nsites,
	     rwkTajD(thetaW.nsam, thetaW.s, tw_val, pi_val));
      
      thetaW.reset(&thetaW);
      thetaPi.reset(&thetaPi);
      strcpy(chr, array[fcol]);
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
    thetaW.add(&thetaW, s);
    thetaPi.add(&thetaPi, pi);
    iwin++;
    stop_region = stoppos;
  }
  
  tw_val = thetaW.eval(&thetaW);
  pi_val = thetaPi.eval(&thetaPi);
  printf("%s\t%lld\t%lld\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",
	 chr, start_region, stop_region,
	 thetaW.nsam, thetaW.nsites, thetaW.s,
	 tw_val / thetaW.nsites, pi_val / thetaPi.nsites,
	 rwkTajD(thetaW.nsam, thetaW.s, tw_val, pi_val));

  // Free memory
  free(array);
  
  return 0; }
