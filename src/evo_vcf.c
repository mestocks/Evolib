#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "evo_vcf.h"


void evoVcfAttach(struct VariantCallFormat *pvcf, char **array, int ncols) {
  pvcf->CHROM = array[0];
  pvcf->POS = atoi(array[1]);
  pvcf->ID = array[2];
  pvcf->REF = array[3];
  pvcf->ALT = array[4];
  pvcf->QUAL = array[5];
  pvcf->FILTER = array[6];
  pvcf->INFO = array[7];
  pvcf->FORMAT = array[8];
  pvcf->SAMPLES = &array[9];
  pvcf->nsamples = ncols - 9;  
}
