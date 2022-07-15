/*
 * =====================================================================================
 *
 *       Filename:  polymer.h
 *
 *    Description:  Biopolymer 
 *
 *        Version:  1.0
 *        Created:  05/09/21 08:54:44 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#ifndef  __polymer_H__
#define  __polymer_H__

#include "biodefs.h"
#include "bioio.h"
#define MAX_HYDRO (15)
struct residue{
      struct atom* atoms;
      int size;
      char name[4];
      char cls[4];
      int  serial;
      int O2p;
      int C1p;
      //struct atom H[MAX_HYDRO];
      //int precursor_index[MAX_HYDRO];
      int numh;
};


//void residue_addh(struct residue* res, int precur_indx, Point3d position, char* h_name);
struct atom* residue_get_atom(struct residue* res, char* atom_loc_name);
void residue_printpdb(FILE* fp, struct residue* res);


struct polymer{
      struct atom* atoms;
      int numatom;
      struct residue* residues;
      int numres;
};


void polymer_create(struct polymer* polymer, struct atom* atoms, int numatom);
void polymer_free(struct polymer* polymer);
void polymer_set_residue_predefined_atoms(struct polymer* poly);
struct residue* residue_at(struct polymer*, int resindex);
int polymer_ressize(struct polymer* poly, int resindex);
void polymer_printpdb(FILE* fp, struct polymer* poly);

#endif   /* ----- #ifndef __polymer_H__  ----- */
