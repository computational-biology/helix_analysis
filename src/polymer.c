/*
 * =====================================================================================
 *
 *       Filename:  polymer.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  16/09/21 04:00:40 PM IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  PARTHAJIT ROY (PR), roy.parthajit@gmail.com
 *   Organization:  The University of Burdwan
 *
 * =====================================================================================
 */

#include "polymer.h"
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  partition_into_residues
 *  Description:  This finction returns an array containing the partition 
 *  		  information.
 *  		  Ex:    |0|15|20|30|35|
 *  		  here   num_residues is 4. first residue starts at 0 and goes up to 14
 *  		  i.e. it contains 15 elements. Second starts at 15 and goes up to 19.
 *  		  The last element starts from 30 goes up to 34. 
 * =====================================================================================
 */
// #define UNDEFINED_VAL (-99999998.0)
// void point3d_uninitialize(Point3d* point)
//{
//      point->x = UNDEFINED_VAL-1.0;
//}
//
//int  Point3d_is_uninit(Point3d* point){
//     if(point->x < UNDEFINED_VAL) return TRUE;
//     else return FALSE;
//}
static void partition_into_residues(const struct atom* atom_array, const int num_atoms, const int init_numres_guess, int** partition_array, int* num_residues)
{ 
     int max_size = init_numres_guess+1;

      int* partition = (int*) malloc ( max_size * sizeof(int) );
      if ( partition == NULL ) {
	    fprintf ( stderr, "\ndynamic memory allocation failed\n" );
	    exit (EXIT_FAILURE);
      }

      int index = 1;
      partition[0] = 0;

      for(int i=1; i<num_atoms; ++i){
	    if(atom_array[i].resid != atom_array[i-1].resid ||
			strcmp(atom_array[i].chain, atom_array[i-1].chain) != 0||
			atom_array[i].model != atom_array[i-1].model || strcmp(atom_array[i].ins, atom_array[i-1].ins) != 0){

		  if(index == max_size){
			max_size *= 2;
			partition = (int*) realloc ( partition, max_size * sizeof(int) );
			if ( partition == NULL ) {
			      fprintf ( stderr, "\ndynamic memory reallocation failed\n" );
			      exit (EXIT_FAILURE);
			}
		  }
		  partition[index] = i;
		  index ++;
	    }
      }
      if(index == max_size){
	    max_size += 1;
	    partition = (int*) realloc ( partition, max_size * sizeof(int) );
	    if ( partition == NULL ) {
		  fprintf ( stderr, "\ndynamic memory reallocation failed\n" );
		  exit (EXIT_FAILURE);
	    }
      }
      partition[index] = num_atoms;
      index++;
      *partition_array = partition;
      *num_residues = index-1;
}
void residue_set_predefined_atoms(struct residue* res)
{
      res->O2p = -1;
      res->C1p = -1;
      for(int i=0; i<res->size; ++i){
	    if(strcmp(res->atoms[i].loc, "O2*") == 0) res->O2p = i;
	    if(strcmp(res->atoms[i].loc, "C1*") == 0) res->C1p = i;
      }
}
void residue_printpdb(FILE* fp, struct residue* res)
{
    for(int i=0; i<res->size; ++i){
        print_pdb_line(fp, res->atoms +i);
    }
//    for(int i=0; i<res->numh; ++i){
//        print_pdb_line(fp, res->H +i);
//    }
}

struct atom* residue_get_atom(struct residue* res, char* atom_loc_name)
{

      
      for(int i=res->size -1; i>=0; --i){  // this is reversed because most base atoms comes at the end of the list.
	    if(strcmp(res->atoms[i].loc, atom_loc_name) == 0){
		  return res->atoms + i;
	    }
      }
      return NULL;
}
//void residue_addh(struct residue* res, int precur_indx, Point3d position, char* h_name)
//{
//      if(Point3d_is_uninit(&position) == TRUE){
//          return;
//      }
//      if(res->numh == MAX_HYDRO){    /* Exception Handling */ 
//	    fprintf(stderr, "Error in function %s()... hydrogen max limit encountered.\n", __func__);
//	    exit(EXIT_FAILURE);
//      }
//      res->H[res->numh] = res->atoms[0];
//
//      res->H[res->numh].center = position;
//      strcpy(res->H[res->numh].loc, h_name);
//      strcpy(res->H[res->numh].symbol, "H");
//      res->H[res->numh].id = 0;
//      res->precursor_index[res->numh] = precur_indx;
//      res->numh = res->numh + 1;
//      return;
//      
//}

void polymer_create(struct polymer* polymer, struct atom* atoms, int numatom){

      polymer->atoms = atoms;
      polymer->numatom = numatom;
      int* res_partition;
      int numres;

      
      partition_into_residues(atoms, numatom, (numatom/5), &res_partition, &numres);
      polymer->numres = numres;
      polymer->residues = (struct residue*) malloc (numres * sizeof(struct residue));
      int atmpos;
      int ressize;
      
      for(int i=0; i<numres; ++i){
	    atmpos = res_partition[i];
	    ressize= res_partition[i+1] - atmpos;


	    polymer->residues[i].atoms = polymer->atoms + atmpos;
//	    print_pdb_line(stdout, polymer->residues[i].atoms);
	    polymer->residues[i].size  = ressize;
	    strcpy(polymer->residues[i].name, polymer->residues[i].atoms[0].resname);
	    polymer->residues[i].numh = 0;
	    polymer->residues[i].serial = i+1;
	    char resname[5];
	     strcpy(resname, polymer->residues[i].name);
	     if(is_guavar(resname) == TRUE){
	         strcpy(polymer->residues[i].cls,"G");
	     }else if(is_adevar(resname) == TRUE){
	         strcpy(polymer->residues[i].cls,"A");
	     }else if(is_cytvar(resname) == TRUE){
	         strcpy(polymer->residues[i].cls,"C");
	     }else if(is_uravar(resname) == TRUE){
	         strcpy(polymer->residues[i].cls,"U");
	     }else{
	         fprintf ( stderr, "\ninvalid residue type encountered in during polymer creation.\n" );
		     exit (EXIT_FAILURE);
		  }
      }
      free(res_partition);
      res_partition = NULL;
}

int polymer_ressize(struct polymer* poly, int resindex)
{
      return poly->residues[resindex].size;
} 

struct residue* residue_at(struct polymer* poly, int resindex){
      return poly->residues + resindex;
}

void polymer_set_residue_predefined_atoms(struct polymer* poly){
      for(int i=0; i<poly->numres; ++i){
	    residue_set_predefined_atoms(poly->residues + i);
      }
}

void polymer_printpdb(FILE* fp, struct polymer* poly)
{
    for(int i=0; i<poly->numres; ++i){
        residue_printpdb(fp, poly->residues + i);
    }
}


void polymer_free(struct polymer* polymer)
{
      free(polymer->residues);
      polymer->residues = NULL;
}
