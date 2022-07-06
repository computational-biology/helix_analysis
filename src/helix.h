/*
 * =====================================================================================
 *
 *       Filename:  helix.h
 *
 *    Description:  This header file manipulates hlx file. 
 *
 *        Version:  1.0
 *        Created:  27/06/21 01:34:03 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */



#ifndef  __helix_H__
#define  __helix_H__

#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "rnabp.h"
#include "spgraph.h"





#define HELIX_MAX 30
#define HAIRPIN_LOOP_MAX 8

enum helix_type{INCR_DECR=0, DECR_INCR, INCR_INCR, DECR_DECR};

struct helix{
      int i[HELIX_MAX];
      int j[HELIX_MAX];
      int seq[HELIX_MAX];
      int size;
      int is_cano[HELIX_MAX];
      int is_hairpin;
      int pinsize;
      char loop_res[HAIRPIN_LOOP_MAX];
      char loopname[64];
};				/* ----------  end of struct helix  ---------- */

struct pseudo_helix{
      struct helix hlx;
      struct pseudo_helix* branch[6];
      int count;
};

void helix_get_bpname(struct helix* hlx, struct nucbp* nbp, int i, char* bpname);

int get_pair_index(struct nucbp* nbp, int i, int j);
void pseudo_helix_init(struct pseudo_helix** helix, int numres);
void pseudo_helix_free(struct pseudo_helix* self);
void pseudo_helix_compute(struct pseudo_helix* all_pseudo_helix[], int* counter, struct nucbp* nbp, int nres);
void pseudo_helix_fprint(struct pseudo_helix* all_pseudo_helix[], int counter, struct nucbp* nbp, int nres, FILE* fp);

//void helix_printpdb(struct helix* self, int hlxcount, struct nucbp* nbp, int nres, struct polymer* poly);
void helix_init(struct helix** self, int numres);
void helix_create_seq(struct helix* self, struct nucbp* nbp, int* seq, char* rule);
void helix_compute(struct helix* self, int* hlxcount, struct nucbp* nbp, int nres);
void helix_fprint(struct helix* self, struct nucbp* nbp, int count, FILE* fp);
void helix_free(struct helix* self);

void helix_free_all(struct helix* self);






#endif   /* ----- #ifndef __helix_H__  ----- */
