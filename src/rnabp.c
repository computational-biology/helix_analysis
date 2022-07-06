/*
 * =====================================================================================
 *
 *       Filename:  rnabp.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  Thursday 31 March 2022 06:25:34  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include "rnabp.h"

int is_canonical(char res1, char res2, char* bpname){
       if(strcmp(bpname, "W:WC") != 0) return FALSE;

       if(res1 == 'G' && res2 == 'C') return TRUE;
       if(res1 == 'C' && res2 == 'G') return TRUE;

       if(res1 == 'G' && res2 == 'U') return TRUE;
       if(res1 == 'U' && res2 == 'G') return TRUE;

       if(res1 == 'A' && res2 == 'U') return TRUE;
       if(res1 == 'U' && res2 == 'A') return TRUE;

       return FALSE;
 }
void rnabp_scan_out(struct nucbp *self, int numres, char *outfile)
{
      FILE* fp = fopen(outfile, "r");
      if(fp == NULL){    /* Exception Handling */
	    fprintf(stderr, "Error in %s at line %d... Unable to open file %s\n", __FILE__, __LINE__, outfile);
	    exit(EXIT_FAILURE);
      }
      char line[256];
      while(fgets(line, 256, fp) != NULL){
	    if(line[0] == '#') continue;
	    break; /* pass the lines */
      }

      int i;
      char sep[] = "\t \n";
      char *token;
      while(line != NULL){
	    if(strlen(line)<=6) break;

	    token = strtok(line, sep);
	    int i = atoi(token);
	    i --;
	    self[i].numbp = 0;

	    token = strtok(NULL, sep);
	    self[i].cifid = atoi(token);

	    token = strtok(NULL, sep);
	    strcpy(self[i].resname, token);
	    self[i].resclass= get_nuc_class(self[i].resname);

	    token = strtok(NULL, sep);
	    if(strlen(token) != 1){
		  fprintf(stderr,"Error... ins-code having more than one characters\n");
		  exit(EXIT_FAILURE);
	    }
	    self[i].ins = token[0];

	    token = strtok(NULL, sep);
	    strcpy(self[i].chain, token);

	    /* For base pairs */
	    while((token = strtok(NULL, sep)) != NULL){

		  int numbp = self[i].numbp;
		  int j = atoi(token);
		  j --;

		  self[i].oth_base_index[numbp] = j;


		  token = strtok(NULL, sep);

		  token = strtok(NULL, sep);

		  token = strtok(NULL, sep);

		  token = strtok(NULL, sep);

		  token = strtok(NULL, sep);
		  strcpy(self[i].name[numbp], token);

		  token = strtok(NULL, sep);
		  strcpy(self[i].type[numbp], token);

		  token = strtok(NULL, sep);
		  self[i].eval[numbp] = atof(token);
		  self[i].numbp++;

	    }
	    fgets(line, 256, fp);

      }
      fclose(fp);
}


int get_numres(char* outfile)
{
      FILE* fp = fopen(outfile, "r");
      if(fp == NULL){    /* Exception Handling */
	    fprintf(stderr, "Error in %s at line %d... Unable to open file %s\n", __FILE__, __LINE__, outfile);
	    exit(EXIT_FAILURE);
      }
      char line[256];
      int flag = 0;
      int nres = -1;
      while(fgets(line, 256, fp) != NULL){
	    if(line[0] != '#') break;
	    if(strncmp(line, "#HEADER   Cleaned number", 24) != 0) continue;
	    nres = atoi(line+38);
	    flag = 1;
	    break;
      }
      if(flag == 0){
	    fprintf(stderr, "Error.... Something wrong in out file. No cleaned residue entry found.\n");
	    fclose(fp);
	    exit(EXIT_FAILURE);
      }
      return nres;
}

