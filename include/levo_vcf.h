
struct VariantCallFormat {
  char *CHROM;
  int POS;
  char *ID;
  char *REF;
  char *ALT;
  char *QUAL;
  char *FILTER;
  char *INFO;
  char *FORMAT;
  char **SAMPLES;

  int nsamples;

  void (*attach)(struct VariantCallFormat *, char **, int);
};

void _vcf_attach(struct VariantCallFormat *pvcf, char **array, int ncols) {

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

char *getGT(char *GT, char *ptr) {
  char *tmp;
  tmp = ptr;
  //char GT[10];
  int c = 0;
  while (*tmp) {
    if (':' == *tmp) {
      GT[c] = '\0';
      break;
    } else {
      GT[c] = *tmp;
    }
    c++;
    tmp++;
  }
  //return GT;
}
