/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  Thursday 30 June 2022 09:32:55  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "helix.h"
#include "rnabp.h"
#include "hlxseq.h"
#include "struct.h"
#include "char_editdist.h"

void needleman_wunsch_sec_seq(char* file1, int numres1, char* file2, int numres2){
      char filepath[512];
      char basename[64];
      char ext[20];
      char infile[512];
      char seqfile[515];
      fname_split(filepath, basename, ext, file1);
      fname_join(infile, filepath, basename, ".dat");
      
      char* s_aligned;
      struct fasta seq1;
      fasta_init(&seq1, infile, numres1);
      

      fname_split(filepath, basename, ext, file2);
      fname_join(infile, filepath, basename, ".dat");
      
      char* t_aligned;
      struct fasta seq2;
      fasta_init(&seq2, infile, numres2);


      char_Needleman_Wunsch_seq_align(seq1.secseq, seq2.secseq, &s_aligned, & t_aligned);

      fname_join(seqfile, filepath, basename, ".ssa");
      
      FILE	*fp;										/* output-file pointer */

      fp	= fopen( seqfile, "w" );
      if ( fp == NULL ) {
	    fprintf ( stderr, "couldn't open file '%s'; %s\n",
			seqfile, strerror(errno) );
	    exit (EXIT_FAILURE);
      }
      fprintf(fp, "%s\n", s_aligned);
      fprintf(fp, "%s\n", t_aligned);
      if( fclose(fp) == EOF ) {			/* close output file   */
	    fprintf ( stderr, "couldn't close file '%s'; %s\n",
			seqfile, strerror(errno) );
	    exit (EXIT_FAILURE);
      }


}


int main(int argc, char* argv[])
{

      
      struct hlxinfo rna1;
      struct hlxinfo rna2;
      
      char rule[6];
      strcpy(rule, argv[1]);

      hlxinfo_create(&rna1, rule, argv[2]);
      hlxinfo_create(&rna2, rule, argv[3]);
      hlxseq_generate(&rna1, &rna2, argv[2], argv[3]);
      hlx_needleman_wunsch_generate(&rna1, &rna2, argv[2], argv[3]);
      needleman_wunsch_sec_seq(argv[2], rna1.numres, argv[3], rna2.numres);



}

