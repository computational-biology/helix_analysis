/*
 * =====================================================================================
 *
 *       Filename:  hlxseq.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  Friday 01 July 2022 09:18:41  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#ifndef  __hlxseq_H__
#define  __hlxseq_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "biodefs.h"
#include "bioio.h"
#include "polymer.h"
#include "rnabp.h"
#include "helix.h"
#include "editdist.h"

struct hlxinfo{
      struct nucbp* rnabp;
      int numres;
      struct helix* helix;
      int hlxcount;

      struct polymer* poly;

};

struct lcsubstr{
      int* s;
      int* t;
      int slen;
      int tlen;
      int maxlen;
      int longest;
      int sindex[40];
      int tindex[40];
      int numstr;
};

void hlxinfo_free(struct hlxinfo* hlx);


void hlxinfo_create(struct hlxinfo* hlx, char* rule, char* infile);
void hlxseq_generate(struct hlxinfo* hlx1, struct hlxinfo* hlx2, char* file1, char* file2);
void hlx_needleman_wunsch_generate(struct hlxinfo* hlx1, struct hlxinfo* hlx2, char* file1, char* file2);



#endif   /* ----- #ifndef __hlxseq_H__  ----- */
