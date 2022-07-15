#ifndef __rnabp_H__
#define __rnabp_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



#define FALSE (0)
#define TRUE (1)



#define MAX_BP   (5)


struct nucbp{
      int cifid;
      char resname[4];
      char resclass;
      char chain[4];
      char ins;
      char name[MAX_BP][5];  //W:WC
      char type[MAX_BP][3];  // BP, TP, BF
      double eval[MAX_BP];

      int oth_base_index[MAX_BP]; // This is added for helix computation.

      int numbp;
};

int is_canonical(char res1, char res2, char* bpname);
 
void rnabp_scan_out(struct nucbp *self, int numres, char *outfile);


int get_numres(char* outfile);




#endif// __BPNET_RNABP_H__
