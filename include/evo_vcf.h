
#ifndef evo_vcf_h__
#define evo_vcf_h__

// Contig0    1    .    A    .    .    .    NCC=17    GT:DP    0/0:0
// 

struct VCFInfo {
  char *id;
  char *value;
  struct VCFInfo *next;
};

struct VCFFormat {
  char *id;
  struct VCFFormat *next;
};

struct VCFSample {
  char *id;
  char *value;
  struct VCFSample *next;
};

struct VCFSample *init_record() {
  struct VCFSample *ptr = calloc(1, sizeof (struct VCFSample));
  ptr->id = NULL;
  ptr->value = NULL;
  ptr->next = NULL;
  return ptr; }

struct VCFSample *create_pool(int poolsize) {
  struct VCFSample *curr;
  struct VCFSample *prev;
  struct VCFSample *head;
  
  head = init_record();
  curr = head;
  prev = head;
  
  for (int i = 1; i < poolsize; i++) {
    curr = init_record();
    prev->next = curr;
    prev = curr;
  }
  return head; }

void print_pool(struct VCFSample *head) {
  struct VCFSample *curr = head;
  while (curr->next != NULL) {
    printf("%p %p\n", curr, curr->next);
    curr = curr->next;
  }
  printf("%p %p\n", curr, curr->next);
}



void free_pool(struct VCFSample *head) {
  struct VCFSample *curr;
  struct VCFSample *prev;
  curr = head;
  while (curr->next != NULL) {
    prev = curr;
    curr = curr->next;
    free(prev);
  }
  free(head);
}

// Create pool of VCFSample records
//
// *-*-*-*
// ^
// pool_head
//
// struct VCFSample *prec;
// record_pool = calloc(poolsize, sizeof (struct VCFSample*));
// struct VCFSample *pool_head = record_pool;
// for (i = 0; i < poolsize; i++) {
//     prec = calloc(1, sizeof (struct VCFSample));
//     
// }


struct VariantCallFormat {
  char *CHROM;
  int POS;
  char *ID;
  char *REF;
  char *ALT;
  char *QUAL;
  char *FILTER;
  struct VCFInfo *INFO;
  struct VCFFormat *FORMAT;
  struct VCFSample **SAMPLES;
  //char **SAMPLES;

  int nsamples;

  void (*attach)(struct VariantCallFormat *, char *, int);
};

void new_attach(struct VariantCallFormat *pvcf, char *buffer, int ncols) {

  int coln;
  char *tmp;
  int nsam;
  char delim = '\t';
  char equal = '=';
  char semi = ';';
  char colon = ':';
  char newline = '\n';

  nsam = ncols - 9;
  
  struct VCFInfo info;
  pvcf->INFO = &info;
  struct VCFInfo *curr = &info;
  struct VCFInfo *head = &info;
  struct VCFFormat format;
  pvcf->FORMAT = &format;
  printf("%p ", pvcf->FORMAT);
  struct VCFFormat *Fcurr = &format;
  struct VCFFormat *Fhead = &format;

  struct VCFSample *current_sample;
  pvcf->SAMPLES = calloc(nsam, sizeof (struct VCFSample*));
  for (int k = 0; k < nsam; k++) {
    pvcf->SAMPLES[k] = calloc(1, sizeof (struct VCFSample));
  }
  
  tmp = buffer;
  coln = 0;
  pvcf->CHROM = tmp;
  while (*tmp && newline != *tmp) {
    if (delim == *tmp) {
      *tmp = '\0';
      coln++;
      tmp++;
      if (coln == 1) pvcf->POS = atoi(tmp);
      else if (coln == 2) pvcf->ID = tmp;
      else if (coln == 3) pvcf->REF = tmp;
      else if (coln == 4) pvcf->ALT = tmp;
      else if (coln == 5) pvcf->QUAL = tmp;
      else if (coln == 6) pvcf->FILTER = tmp;
      else if (coln == 7) {
	curr->id = tmp;
	curr->next = NULL;
      } else if (coln == 8) {
	Fcurr->id = tmp;
	Fcurr->next = NULL;
      } else if (coln > 8) {
	Fcurr = pvcf->FORMAT;
	current_sample = pvcf->SAMPLES[coln - 9];
	current_sample->id = Fcurr->id;	
	current_sample->value = tmp;
	current_sample->next = NULL;
      } 
    } else {
      if (coln == 7) {
	if (equal == *tmp) {
	  *tmp = '\0';
	  tmp++;
	  curr->value = tmp;
	} else if (semi == *tmp) {
	  *tmp = '\0';
	  tmp++;
	  curr->next = calloc(1, sizeof (struct VCFInfo));
	  curr = curr->next;
	  curr->id = tmp;
	  curr->next = NULL;
	} else {
	  tmp++;
	}
      } else if (coln == 8) {
	if (colon == *tmp) {
	  *tmp = '\0';
	  tmp++;
	  Fcurr->next = calloc(1, sizeof (struct VCFFormat));
	  Fcurr = Fcurr->next;
	  Fcurr->id = tmp;
	  Fcurr->next = NULL;
	} else {
	  tmp++;
	}
      } else if (coln > 8) {
	if (colon == *tmp) {
	  *tmp = '\0';
	  tmp++;
	  current_sample->next = init_record();
	  current_sample = current_sample->next;
	  current_sample->value = tmp;
	  current_sample->next = NULL;
	  
	  //Fcurr->next = calloc(1, sizeof (struct VCFFormat));
	  //Fcurr = Fcurr->next;
	  //Fcurr->id = tmp;
	  //Fcurr->next = NULL;
	} else {
	  tmp++;
	}
      } else {
	tmp++;
      }
    }
  }
  *tmp = '\0';
  struct VCFInfo *pinfo = head;
  struct VCFInfo *linfo;
  int i = 0;
  while (pinfo != NULL) {
    //printf("%s %s %p\n", pinfo->id, pinfo->value, pinfo->next);
    linfo = pinfo;
    pinfo = pinfo->next;
    if (i > 0) {
      free(linfo);
    }
    i++;
  }
  struct VCFFormat *pformat = Fhead;
  struct VCFFormat *lformat;
  int j = 0;
  while (pformat != NULL) {
    lformat = pformat;
    pformat = pformat->next;
    if (j > 0) {
      free(lformat);
    }
    j++;
  }
  //printf("%s %s %p %s %s %p\n", head->id, head->value, head->next, curr->id, curr->value, curr->next);
  struct VCFSample *lsam;
  struct VCFSample *psam;
  for (int l = 0; l < nsam; l++) {
    psam = pvcf->SAMPLES[l];
    while (psam != NULL) {
      lsam = psam;
      psam = psam->next;
      free(lsam);
    }
  }
  
  //free(curr);
  //free(Fcurr);
}



//extern void evoVcfAttach(struct VariantCallFormat *pvcf, char **array, int ncols);

void evoVcfAttach(struct VariantCallFormat *pvcf, char **array, int ncols) {
  pvcf->CHROM = array[0];
  pvcf->POS = atoi(array[1]);
  pvcf->ID = array[2];
  pvcf->REF = array[3];
  pvcf->ALT = array[4];
  pvcf->QUAL = array[5];
  pvcf->FILTER = array[6];
  //pvcf->INFO = array[7];
  //pvcf->FORMAT = array[8];
  //pvcf->SAMPLES = &array[9];
  pvcf->nsamples = ncols - 9;  
}

struct VCFsample {
  int DP;
  char GT[10];
};


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
