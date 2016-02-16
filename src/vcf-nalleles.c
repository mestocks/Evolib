#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <evo_vcf.h>

#include <rwk_parse.h>
#include <rwk_htable.h>

int count_columns(char *buffer, char delim) {
  char *tmp = buffer;
  int count = 0;
  while (*tmp) {
    if (delim == *tmp) { count++; }
    tmp++; }
  return count + 1; }

int main(int argc, char **argv) {
  
  struct VCFsample SMP;
  struct VariantCallFormat VCF;
  VCF.attach = evoVcfAttach;
  int coln;
  int ncols;
  int lwidth = 100000;
  char delim = '\t';
  char newline = '\n';
  char buffer[lwidth];

  char *tmp;
  char **array;
  
  while (fgets(buffer, sizeof(buffer), stdin)) {
    
    if (buffer[0] == '#' && buffer[1] != '#') {
      ncols = count_columns(buffer, delim);
      array = calloc(ncols, sizeof (char*));
    } else if (buffer[0] != '#') {
      rwkStrtoArray(array, buffer, &delim);
      VCF.attach(&VCF, array, ncols);
    
      int ref = 0;
      int alt = 0;
      for (int i = 0; i < VCF.nsamples; i++) {
	SMP.GT[0] = '\0';
	getGT(&SMP, VCF.SAMPLES[i]);
	if (strcmp("0/1", SMP.GT) == 0) {
	  ref++;
	  alt++;
	} else if (strcmp("1/1", SMP.GT) == 0) {
	  alt += 2;
	} else if (strcmp("0/0", SMP.GT) == 0) {
	  ref += 2;
	}
      }
      if (strlen(VCF.REF) == 1 && strlen(VCF.ALT) == 1) {
	printf("%s\t%d\t%d\t%s\t%d\t%d\n", VCF.CHROM, VCF.POS - 1, VCF.POS,
	       "nalleles", ref, alt);
      } else {
	printf("%s\t%d\t%d\t%s\t%d\t%d\n", VCF.CHROM, VCF.POS - 1, VCF.POS,
	       "nalleles", 0, 0);
      }
    }
  }
  
  free(array);
  
  return 0; }
