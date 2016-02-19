/*

... | dna2codon

chr    pos-1    pos    nuc

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <rwk_parse.h>

int main(int argc, char **argv) {

  int ncols = 7;
  int lcols = 1024;
  int lwidth = 2048;
  char delim = '\t';
  char newline = '\n';
  char buffer[lwidth];
  
  char **array;
  array = calloc(ncols, sizeof (char*));

  int start = 0;
  char codon[4];
  codon[3] = '\0';
  
  while (fgets(buffer, sizeof(buffer), stdin)) {
    rwkStrtoArray(array, buffer, &delim);
    codon[start] = *array[4];
    if (start == 2) {
      printf("%s\t%d\t%s\t%s\n", array[0], atoi(array[1]) - 2, array[2], codon);
      start = 0;
    } else {
      start++;
    }
  }
  
  free(array);
  
  return 0; }
