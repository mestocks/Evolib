#ifndef evo_vcf_h__
#define evo_vcf_h__

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



//extern void evoVcfAttach(struct VariantCallFormat *pvcf, char **array, int ncols);

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

struct VCFsample {
  int DP;
  char GT[10];
};

/*

  GT:DP

  

 */




void parse_sample(struct VCFsample *SMP, char *format, char *sample) {

  int fi, si;
  char *ftmp, *stmp;
  int findex, sindex;

  char fchar[1000];
  char schar[1000];

  char newline = '\n';

  fi = 0;
  si = 0;
  findex = 0;
  sindex = 0;
  ftmp = format;
  stmp = sample;

  while (*stmp && *ftmp) {
    if (findex == sindex) {
      if (':' == *ftmp) {
	fchar[fi] = '\0';
	fi = 0;
	findex++;
      } else {
	fchar[fi] = *ftmp;
	fi++;
      }
      ftmp++;
    } else {
      if (':' == *stmp) {
	schar[si] = '\0';
	si = 0;
	sindex++;
	if (strcmp("DP", fchar) == 0) {
	  SMP->DP = atoi(schar);
	  printf("\n - %s %s\n", schar, fchar);
	}
	if (strcmp("GT", fchar) == 0) {
	  strcpy(SMP->GT, schar);
	  printf("\n - %s %s\n", schar, fchar);
	}
      } else {
	schar[si] = *stmp;
	si++;
      }
      stmp++;
    }
    printf("\n%s %s\n", ftmp, stmp);
  }
  //fchar[fi] = '\0';
  //fi = 0;
  //findex++;
  //schar[si] = '\0';
  //si = 0;
  //sindex++;
  //printf("::%s\n", fchar);
  if (strcmp("DP", fchar) == 0) {
    SMP->DP = atoi(schar);
    printf("\n -- %s %s\n", schar, fchar);
  }
  if (strcmp("GT", fchar) == 0) {
    strcpy(SMP->GT, schar);
    printf("\n -- %s %s\n", schar, fchar);
  }
}

char *getGT(struct VCFsample *SMP, char *ptr) {
  char *tmp;
  tmp = ptr;
  //char GT[10];
  char dp[10];
  int c = 0;
  int f = 0;
  while (*tmp) {
    if (':' == *tmp) {
      if (f == 0) {
	SMP->GT[c] = '\0';
      } else if (f == 1) {
	dp[c] = '\0';
	SMP->DP = atoi(dp);
      }
      c = 0;
      f++;
    } else {
      if (f == 0) {
	SMP->GT[c] = *tmp;
      } else if (f == 1) {
	dp[c] = *tmp;
      }
      c++;
    }
    tmp++;
  }
  if (f == 0) {
    SMP->GT[c] = '\0';
  } else if (f == 1) {
    dp[c] = '\0';
    SMP->DP = atoi(dp);
  }
  //return GT;
}

char *old_getGT(char *GT, char *ptr) {
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

#endif
