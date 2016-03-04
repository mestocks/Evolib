#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <evo_vcf.h>

#include <rwk_parse.h>
#include <rwk_htable.h>



int main(int argc, char **argv) {

  int minDP;
  
  if (argc == 1) {
    minDP = 0;
  } else {
    minDP = atoi(argv[1]);
  }
    
  char delim = '\t';
  char newline = '\n';
  int lwidth = 100000;
  char buffer[lwidth];
  
  int coln, ncols;
  int nref, nalt;
  
  struct VCFsample SMP;
  struct VariantCallFormat VCF;
  VCF.attach = new_attach;
  
  while (fgets(buffer, sizeof(buffer), stdin)) {
  
    if (buffer[0] == '#' && buffer[1] != '#') {      
      ncols = rwkCountCols(buffer, delim);

      // Create pool of nodes
      
      struct VCFSample *head;
      struct VCFSample *curr;
      struct VCFSample *new_node;

      int nrecs = 10;      
      head = create_pool(nrecs);

      print_pool(head);

      printf("######\n");

      // Cleave nodes from pool
      
      struct VCFSample *new_head;
      struct VCFSample *new_tail;
      struct VCFSample *new_curr;
      
      int nnodes = 2;
      new_head = head;
      new_curr = new_head;

      while (nnodes > 1) {
	new_curr = new_curr->next;
	new_tail = new_curr;
	head = new_curr->next;
	nnodes--;
      }
      new_curr->next = NULL;

      print_pool(new_head);
      
      printf("######\n");

      print_pool(head);

      printf("######\n");

      // Return nodes to pool
      
      new_tail->next = head;
      head = new_head;
      
      printf("######\n");

      print_pool(head);
      
      free_pool(head);
      //free_pool(new_curr);

      
      //VCF.SAMPLES = calloc(ncols, sizeof (struct VCFSample*));
      
    } else if (buffer[0] != '#') {
      
      
      //printf("%s", buffer);
      //VCF.attach(&VCF, buffer, ncols);
      //printf("%p %p - %s\n", VCF.FORMAT->id, VCF.SAMPLES[0]->id, VCF.SAMPLES[0]->value);
    }
  }
  
  return 0; }
