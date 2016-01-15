#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int getDP(char *sample, int lwidth) {
  int dp;
  int cindex = 0;
  int windex = 0;
  char delim = ':';
  char *tmp = sample;

  char col[20];
  
  while (*tmp) {
    if (delim == *tmp) {
      col[windex] = '\0';
      cindex++;
      windex = 0;
    }
    else {
      col[windex] = *tmp;
    }
    if (cindex == 3) {
      dp = atoi(col);
    }
    tmp++;
  }
  if (cindex == 2) {
    dp = atoi(col);
  }
  return dp;
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
    } else {
      col[windex] = *tmp;
      windex++;
    }
    tmp++;
  }
  free(col);
}

int main(int argc, char **argv) {

  int ncols = 20;
  const int lwidth = 10000;
  char delim = '\t';
  char buffer[lwidth];

  // pointer to <ncols> pointers to char arrays
  char **array = calloc(ncols, sizeof (char*));

  for (int i = 0; i < ncols; i++) {
    array[i] = calloc(lwidth, sizeof (char));
  }
  
  while (fgets(buffer, sizeof(buffer), stdin)) {
    if (buffer[0] != '#') {
      parse_columns(array, ncols, lwidth, buffer, delim);
      //printf("%s %s %s %s %s\n", array[0], array[1], array[8], array[9], array[10]);
      printf("%s %s", array[0], array[1]);
      for (int k = 9; k < ncols; k++) {
	int dp = getDP(array[k], lwidth);
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
