
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
};

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
