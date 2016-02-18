#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <evo_vcf.h>

#include <rwk_parse.h>
#include <rwk_htable.h>

int main(int argc, char **argv) {

  char f0[] = "GT:AD:DP";
  char f1[] = "GT:AD:DP:GQ:PL";
  char f2[] = "GT:DP";

  int minDP = 8;

  
  int dpi;
  char smpl[1024];
  char **samparray;
  samparray = calloc(24, sizeof (char*));
  char smpdel = ':';
  
  char delim = '\t';
  char newline = '\n';
  int lwidth = 100000;
  char buffer[lwidth];

  char **array;
  int coln, ncols;
  int nref, nalt;
  
  struct VCFsample SMP;
  struct VariantCallFormat VCF;
  VCF.attach = evoVcfAttach;

  while (fgets(buffer, sizeof(buffer), stdin)) {
    
    if (buffer[0] == '#' && buffer[1] != '#') {      
      ncols = rwkCountCols(buffer, delim);
      array = calloc(ncols, sizeof (char*));

    } else if (buffer[0] != '#') {
      rwkStrtoArray(array, buffer, &delim);
      VCF.attach(&VCF, array, ncols);

      if (strcmp(f0, VCF.FORMAT) == 0) {
	dpi = 2;
      } else if (strcmp(f1, VCF.FORMAT) == 0) {
	dpi = 2;
      } else if (strcmp(f2, VCF.FORMAT) == 0) {
	dpi = 1;
      }
      
      nref = 0;
      nalt = 0;
      if (strlen(VCF.REF) == 1 && strlen(VCF.ALT) == 1) {
	
	for (int i = 0; i < VCF.nsamples; i++) {
	  strcpy(smpl, VCF.SAMPLES[i]);
	  rwkStrtoArray(samparray, smpl, &smpdel);
	  if (atoi(samparray[dpi]) >= 1) {
	    SMP.GT[0] = '\0';
	    getGT(&SMP, VCF.SAMPLES[i]);
	    if (strcmp("0/1", SMP.GT) == 0) {
	      nref++;
	      nalt++;
	    } else if (strcmp("1/1", SMP.GT) == 0) {
	      nalt += 2;
	    } else if (strcmp("0/0", SMP.GT) == 0) {
	      nref += 2;
	    }
	  }
	}
      }
      printf("%s\t%d\t%d\t%s\t%d\t%d\n", VCF.CHROM, VCF.POS - 1, VCF.POS,
	     "nalleles", nref, nalt);
    }
  }
  
  free(array);
  
  return 0; }
