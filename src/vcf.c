#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct FORMAT {
  int nitems;
  char **parray;
};

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

int count_columns(char *buffer, char delim) {
  char *tmp = buffer;
  int count = 0;
  while (*tmp) {
    if (delim == *tmp) { count++; }
    tmp++; }
  return count + 1; }

int main(int argc, char **argv) {

  struct VariantCallFormat VCF;
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
      printf("%d\n", ncols);
      array = calloc(ncols, sizeof (char*));
    } else if (buffer[0] != '#') {
      tmp = buffer;
      array[0] = tmp;
      coln = 0;
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

      VCF.CHROM = array[0];
      VCF.POS = atoi(array[1]);
      VCF.ID = array[2];
      VCF.REF = array[3];
      VCF.ALT = array[4];
      VCF.QUAL = array[5];
      VCF.FILTER = array[6];
      VCF.INFO = array[7];
      VCF.FORMAT = array[8];
      VCF.SAMPLES = &array[9];
      VCF.nsamples = ncols - 9;
      
      printf("%s\t%d", VCF.CHROM, VCF.POS);
      for (int i = 0; i < VCF.nsamples; i++) { printf("\t%s", VCF.SAMPLES[i]); }
      printf("\n");
    }
  }
  //
  
  free(array);
  
  return 0; }

int getDP(char *sample, int findex) {
  int dp;
  int cindex = 0;
  int windex = 0;
  char delim = ':';
  char *tmp = sample;

  char col[20];
  
  while (*tmp) {
    if (delim == *tmp) {
      if (cindex == findex) {
	break;
      }
      cindex++;
      windex = 0;
    }
    else {
      col[windex] = *tmp;
      windex++;
    }
    tmp++;
  }
  col[windex] = '\0';
  if (cindex == findex) {
    dp = atoi(col);
  }
  return dp;
}

int get_index(char *format, char *value, char delim) {
  int cindex = 0;
  int windex = 0;
  char *tmp = format;
  char col[20];
  while (*tmp) {
    if (delim == *tmp) {
      col[windex] = '\0';
      if (strcmp(value, col) == 0) {
	return cindex;
      }
      cindex++;
      windex = 0;
    }
    else {
      col[windex] = *tmp;
      windex++;
    }
    tmp++;
  }
}

void parse_columns(char **array, int ncols, int lwidth,
		   char *buffer, char delim) {
  
  int cindex = 0;
  int windex = 0;
  char *tmp = buffer;
  char *col = calloc(lwidth, sizeof (char));

  while (*tmp) {
    if (delim == *tmp) {
      col[windex] = '\0';
      strcpy(array[cindex], col);
      cindex++;
      windex = 0;
    } else {;
      col[windex] = *tmp;
      windex++;
    }
    tmp++;
  }
  free(col);
}


int main_old(int argc, char **argv) {

  int dp;
  int dpindex;
  int ncols;
  const int lwidth = 10000;
  char delim = '\t';
  char dformat = ':';
  char buffer[lwidth];

  char format[] = "GT";
  
  // pointer to <ncols> pointers to char arrays
  char **array;
  
  while (fgets(buffer, sizeof(buffer), stdin)) {
    
    if (buffer[0] == '#' && buffer[1] != '#') {
      ncols = count_columns(buffer, delim);
      array = calloc(ncols, sizeof (char*));
      for (int i = 0; i < ncols; i++) {
	array[i] = calloc(lwidth, sizeof (char));
      }
    } else if (buffer[0] != '#') {
      parse_columns(array, ncols, lwidth, buffer, delim);
      
      dpindex = get_index(array[8], format, dformat);
      //printf("%s %d\n", buffer, dpindex);
      printf("%s %s", array[0], array[1]);
      for (int k = 9; k < ncols; k++) {
	dp = getDP(array[k], dpindex);
	printf(" %d", dp);
      }
      printf("\n");
    }

  }

  
  for (int j = 0; j < ncols; j++) {
    free(array[j]);
  }
  free(array);
  
  return 0; }