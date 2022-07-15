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



#include "helix.h"


struct stack{
      void** elem;
      int top;
      int smax;
};
void stack_init(struct stack* self, int size){
      self->smax = size;
      self->elem = (void**) malloc (size * sizeof(void*));
      self->top  = -1; 
}

void stack_push(struct stack* self, void* val){
      if(self->top == self->smax-1){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s. (File: %s, Line %d)... Stack overflows.\n", __func__, __FILE__, __LINE__);
	    exit(EXIT_FAILURE);
      }
      self->elem[++self->top] = val;
}
void* stack_pop(struct stack* self){
      return self->elem[self->top--];
}
void stack_free(struct stack* self){
      free(self->elem);
}

void helix_get_bpname(struct helix* hlx, struct nucbp* nbp, int i, char* bpname)
{

      int base1 = hlx->i[i];
      int base2 = hlx->j[i];
      int indx = get_pair_index(nbp, base1, base2);
      bpname[0] = nbp[base1].resclass;
      bpname[1] = '-';
      bpname[2] = nbp[base2].resclass;
      bpname[3] = ' ';
      strcpy(bpname+4, nbp[base1].name[indx]);
}


int convert_res(char res)
{
      switch(res){
	    case 'A': 
	    case 'a': return 1;
	    case 'G': 
	    case 'g': return 2;
	    case 'C': 
	    case 'c': return 3;
	    case 'U': 
	    case 'u': return 4;
	    case '#': return 0; // # means not considered for matching.
	    default:
		         /* Exception Handling */ 
			    fprintf(stderr, "Error in function %s()... \n", __func__);
			    exit(EXIT_FAILURE);
      }
}

int convert_edge(char edge)
{
      switch(edge){
	    case 'W':
	    case 'w':
	    case '+': return 1;
	    case 'S':
	    case 's':
	    case 'z': return 2;
	    case 'H':
	    case 'h':
	    case 'g': return 3;
	    case '#': return 0; // # means not considered for matching.
	    default:
		         /* Exception Handling */ 
			    fprintf(stderr, "Error in function %s()... \n", __func__);
			    exit(EXIT_FAILURE);
      }
}

int create_bp_integer(int res1, int res2, int edge1, int edge2, int ori)
{
      int num = res1;
      num = 10 * num + res2;
      num = 10 * num + edge1;
      num = 10 * num + edge2;
      num = 10 * num + ori;
      return num;

}



int get_seq(char res1, char res2, char edge1, char edge2, char ori, char* rule)
{
      int orival = 0;
      if(rule[0] == 'F'){
	    res1 = '#';
	    res2 = '#';
      }
      if(rule[1] == 'F'){
	    edge1 = '#';
	    edge2 = '#';
      }
      if(rule[2] == 'F'){
	    orival = 0;
      }else{
	    if(ori == 'C'){
		  orival = 1;
	    }else{
		  orival = 2;
	    }
      }
      if(rule[3] == 'F'){
          if(res1 == 'A') res1 = 'G';
          if(res1 == 'U') res1 = 'C';
          if(res2 == 'A') res2 = 'G';
          if(res2 == 'U') res2 = 'C';
      }

      int res1val = convert_res(res1);
      int res2val = convert_res(res2);
      int edge1val = convert_edge(edge1);
      int edge2val = convert_edge(edge2);
      int val = create_bp_integer(res1val, res2val, edge1val, edge2val, orival);
      return val;


}


void helix_create_seq(struct helix* self, struct nucbp* nbp, int* seq, char* rule)
{
      for(int k=0; k<self->size; ++k){
	    seq[k] = get_seq(nbp[self->i[k]].resclass, nbp[self->j[k]].resclass, nbp[self->i[k]].name[0][0], nbp[self->i[k]].name[0][2], nbp[self->i[k]].name[0][3], rule);
      }
}

void helix_printf(struct helix* self, struct nucbp* nbp, FILE* fp){
      if(self->is_hairpin == TRUE){
	    fprintf(fp, "TYPE  STEM-LOOP  with (%s) loop structure.\n", self->loopname);
      }else{
	    fprintf(fp, "TYPE  PURE HELIX structure.\n");
      }
      for(int k=0; k<self->size; ++k){
	    if(self->is_hairpin == TRUE){
		  fprintf(fp, "STEM  ");
	    }else{
		  fprintf(fp, "HELX  ");
	    }
	    fprintf(fp, "%c %c  ", nbp[self->i[k]].resclass, nbp[self->j[k]].resclass);
	    fprintf(fp, "%5d %-4s", nbp[self->i[k]].cifid,nbp[self->i[k]].chain) ;
	    fprintf(fp, "%5d %-4s", nbp[self->j[k]].cifid,nbp[self->j[k]].chain) ;
	    fprintf(fp, "%5d:%-5d  ", self->i[k]+1, self->j[k]+1);
	    fprintf(fp, "\n");
      }
      if(self->is_hairpin == TRUE){
	    fprintf(fp, "LOOP  ");
	    for(int i=0; i<self->pinsize; ++i){
		  fprintf(fp, "%c",self->loop_res[i]);
	    }
	    fprintf(fp, "\n");
      }
      fprintf(fp, "#\n");
}


int get_pair_index(struct nucbp* nbp, int i, int j){ 
      for(int k=0; k<nbp[i].numbp; ++k){
	    if(nbp[i].oth_base_index[k] == j) {
		  return k;
	    }
      }
      return -1;
}
int is_paired(struct nucbp* nbp, int i, int j, int* pair_type){ // pair_type 1 for canonical, 0 for non canonical.
      for(int k=0; k<nbp[i].numbp; ++k){
	    if(nbp[i].oth_base_index[k] == j) {

		  if(is_canonical(nbp[i].resclass, nbp[j].resclass, nbp[i].name[k]) == TRUE){
			*pair_type = 1; 
		  }else{
			*pair_type = 0;
		  }
		  return 1;
	    }
      }
      return 0;
}

int is_next_resi(struct nucbp* nbp, int nres, int i, int nexti){
      if(i<0 || i >=nres || nexti<0 || nexti >=nres) return FALSE;
      if(nexti-i > 1) return FALSE;
      if(strcmp(nbp[i].chain, nbp[nexti].chain) != 0) return FALSE;
      else if(nbp[i].cifid < nbp[nexti].cifid) return TRUE;
      else if(nbp[i].cifid == nbp[nexti].cifid && nbp[i].ins <  nbp[nexti].ins) return TRUE;
      else return FALSE;
}

int is_prev_resi(struct nucbp* nbp, int nres, int i, int previ){
      if(i<0 || i >=nres || previ<0 || previ >=nres) return FALSE;
      if(i-previ >1) return FALSE;
      if(strcmp(nbp[i].chain, nbp[previ].chain) != 0) return FALSE;
      else if(nbp[i].cifid >  nbp[previ].cifid) return TRUE;
      else if(nbp[i].cifid == nbp[previ].cifid && nbp[i].ins > nbp[previ].ins) return TRUE;
      else return FALSE;
}




/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  helix_comp
 *  Description:  This function computes pure helix from
 *  		  rnabp data base.  
 * =====================================================================================
 */

int helix_calc(struct helix* self, struct nucbp* nbp, int nres, struct graph* g, int i, int j, enum helix_type htype){
      /* The base case */
      if(graph_edge_index(g, i, j) != -1) return FALSE;
      int is_cano;
      if(is_paired(nbp, i, j, &is_cano) == FALSE) return FALSE;
      else{
	    graph_set_edge(g, i, j);
	    graph_set_edge(g, j, i);
	    
	    
	    
	    
	    if(self->size == HELIX_MAX){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s. (File: %s, Line %d)... Helix Overflow Encountered.\n", __func__, __FILE__, __LINE__);
		  exit(EXIT_FAILURE);
	    }	       
	    self->i[self->size]=i;
	    self->j[self->size]=j;
	    self->is_cano[self->size]=is_cano;
	    self->size ++;

	    
	    if(htype == INCR_DECR){
		  if((is_next_resi(nbp, nres, i, i+1) == TRUE && is_prev_resi(nbp, nres, j, j-1) == TRUE)){

			helix_calc(self, nbp, nres, g, i+1, j-1, htype);
		  } 
	    }else if(htype == DECR_INCR){
		  if((is_prev_resi(nbp, nres, i, i+1) == TRUE && is_next_resi(nbp, nres, j, j-1) == TRUE)){
			helix_calc(self, nbp, nres, g, i+1, j-1, htype);

			
		  } 
	    }else if(htype == INCR_INCR){
		  if((is_next_resi(nbp, nres, i, i+1) == TRUE && is_next_resi(nbp, nres, j, j-1) == TRUE)){
			
			
			helix_calc(self, nbp, nres, g, i+1, j-1, htype);

			
		  } 
	    }else{// DECR_DECR
		  if((is_prev_resi(nbp, nres, i, i+1) == TRUE && is_prev_resi(nbp, nres, j, j-1) == TRUE)){

			
			helix_calc(self, nbp, nres, g, i+1, j-1, htype);

			
		  } 
	    }
      }

      return TRUE;
}		/* -----  end of function helix_calc  ----- */


struct pseudo_helix* pseudo_helix_getnode(){
      struct pseudo_helix* node;
      node = (struct pseudo_helix*) malloc (sizeof(struct pseudo_helix));
      node->hlx.size = 0;
      node->count = 0;
      for(int i=0; i<6; ++i){
	    node->branch[i] = NULL;
      }
      return node;
}
void pseudo_helix_free(struct pseudo_helix* self){
      for(int i=0; i<self->count; ++i){
	    pseudo_helix_free(self->branch[i]);
      }
      free(self);
      return;
}
struct pseudo_helix* pseudo_helix_calc(struct nucbp* nbp, int nres, struct graph* g, int i, int j, enum helix_type htype){
      struct pseudo_helix * self = pseudo_helix_getnode();
      struct helix* h = &self->hlx;
      int hval;
      hval = helix_calc(h, nbp, nres, g, i, j, htype);
      if(self->hlx.size < 3){
	    self->hlx.size = 0;
	    free(self);
	    return NULL;
      }
      
      int  i_last = h->i[h->size - 1];
      int  j_last = h->j[h->size - 1];
      if((is_next_resi(nbp, nres, i_last, i_last+1) == TRUE) || (is_prev_resi(nbp, nres, i_last, i_last+1) == TRUE)){
	    for(int k=0; k<nbp[i_last + 1].numbp; ++k){
		  struct pseudo_helix* tmphlx = pseudo_helix_calc(nbp, nres, g, i_last+1, nbp[i_last +1].oth_base_index[k], htype);
		  if(tmphlx != NULL){
			self->branch[self->count] = tmphlx;
			self->count ++;
		  }
	    }
      }

      if((is_prev_resi(nbp, nres, j_last, j_last - 1) == TRUE)|| (is_next_resi(nbp, nres, j_last, j_last - 1) == TRUE)){
	    for(int k=0; k<nbp[j_last - 1].numbp; ++k){
		  struct pseudo_helix* tmphlx = pseudo_helix_calc(nbp, nres, g, nbp[j_last - 1].oth_base_index[k], j_last - 1, htype);
		  if(tmphlx != NULL){
			self->branch[self->count] = tmphlx;
			self->count ++;
		  }
	    }
      }
      return self;
}

void helix_pseudo_fprint(struct pseudo_helix* self, struct nucbp* nbp, struct stack* stk, FILE* fp){
      stack_push(stk, self);
      if(self->count == 2){
	    fprintf(fp, "More than one branch found...\n");
      }
      if(self->count == 0){
	    struct pseudo_helix* tmphlx = (struct pseudo_helix*) stk->elem[0];
//            struct helix* f_elem = &tmphlx->hlx;
//
//	    for(int i=0; i<=stk->top; ++i){
//		  tmphlx = (pseudo_helix*) stk->elem[i];
//		  fprintf(fp, "%5d~", tmphlx->hlx.i[0] + 1) ;
//		  for(int k=0; k<tmphlx->hlx.size; ++k){
//			fprintf(fp, " ");
//		  }
//	    }
//	    fprintf(fp,"\n");
	    for(int i=0; i<=stk->top; ++i){
		  tmphlx = (struct pseudo_helix*) stk->elem[i];
		  for(int k=0; k<tmphlx->hlx.size; ++k){
			fprintf(fp,"PHLX  ");
			fprintf(fp, "%5d:%-5d  ", tmphlx->hlx.i[k]+1, tmphlx->hlx.j[k]+1);
			fprintf(fp, "%c=%c", nbp[tmphlx->hlx.i[k]].resclass, nbp[tmphlx->hlx.j[k]].resclass);
			fprintf(fp, "  %5d %-4s", nbp[tmphlx->hlx.i[k]].cifid,nbp[tmphlx->hlx.i[k]].chain);
			fprintf(fp, "  %5d %-4s", nbp[tmphlx->hlx.j[k]].cifid,nbp[tmphlx->hlx.j[k]].chain);
			if(i<stk->top && k == tmphlx->hlx.size -1){
			      fprintf(fp," <---Break\n");
			}else{
			      fprintf(fp,"\n");
			}
		  }
	    }
	    fprintf(fp,"#\n\n");

	    
//	    for(int i=0; i<=stk->top; ++i){
//		  tmphlx = (pseudo_helix*) stk->elem[i];
//		  for(int k=0; k<tmphlx->hlx.size; ++k){
//			fprintf(fp, "|");
//		  }
//		  fprintf(fp,"~");
//	    }
//	    fprintf(fp,"\n");
//
//	    
//	    for(int i=0; i<=stk->top; ++i){
//		  tmphlx = (pseudo_helix*) stk->elem[i];
//		  for(int k=0; k<tmphlx->hlx.size; ++k){
//			fprintf(fp, "%c", nbp[tmphlx->hlx.j[k]].resclass);
//		  }
//		  fprintf(fp,"~");
//	    }
      }

	    
     

      
      for(int k=0; k<self->count; ++k){
	    helix_pseudo_fprint(self->branch[k], nbp, stk, fp);
      }
      
      
      stack_pop(stk);

      
}


void pseudo_helix_compute(struct pseudo_helix* all_pseudo_helix[], int* counter, struct nucbp* nbp, int nres)
{

      int cnt = 0;
      struct graph g;
      graph_init(&g, nres, UNDIRECTED); 
      for(int i=0; i<nres; ++i){
	    for(int k=0; k<nbp[i].numbp; ++k){
		  int j = nbp[i].oth_base_index[k];
		  if(graph_edge_index(&g, i, j) != -1) continue;

		  struct pseudo_helix* tmp = pseudo_helix_calc(nbp, nres, &g, i, j, INCR_DECR);
		  if(tmp== NULL){
			tmp = pseudo_helix_calc(nbp, nres, &g, i, j, DECR_INCR);
		  }
		  if(tmp== NULL){
			tmp = pseudo_helix_calc(nbp, nres, &g, i, j, INCR_INCR);
		  }
		  if(tmp == NULL){
			tmp = pseudo_helix_calc(nbp, nres, &g, i, j, DECR_DECR);
		  }
		  if(tmp != NULL){
			if(tmp->count == 0) {
			      free(tmp);
			}else{
			      all_pseudo_helix[cnt] = tmp;
			      cnt ++;
			}

		  }
	    }
      }
      *counter = cnt;
      graph_free(&g);
}

void pseudo_helix_init(struct pseudo_helix** helix, int numres)
{
      
      *helix	= (struct pseudo_helix*) malloc ( (numres/4) * sizeof(struct pseudo_helix) );
      if ( *helix==NULL ) {
	    fprintf ( stderr, "\ndynamic memory allocation failed in function %s()\n" , __func__);
	    exit (EXIT_FAILURE);
      }
}
void pseudo_helix_fprint(struct pseudo_helix* all_pseudo_helix[], int counter, struct nucbp* nbp, int nres, FILE* fp){
      struct stack stk;
      stack_init(&stk, nres/4);

      
      for(int i=0; i<counter; ++i){
	    helix_pseudo_fprint(all_pseudo_helix[i], nbp, &stk, fp);
      }

      
      stack_free(&stk);
}
void helix_init(struct helix** self, int numres){
      int maxsize = numres /4;
      *self = (struct helix*) malloc ( maxsize * sizeof(struct helix) );
      if ( *self == NULL ) {
	    fprintf ( stderr, "\nDynamic memory allocation failed\n" );
	    exit (EXIT_FAILURE);
      }
}

void get_loop_name(char* seq, int len, char* name){
      if(len == 4){
	    if(toupper(seq[0]) == 'G' 
			&&
	       (toupper(seq[2]) == 'G' || toupper(seq[2]) == 'A') 
	                &&
	        toupper(seq[3]) == 'A'){
		  strcpy(name, "GNRA Tetra");
	    }else if(toupper(seq[0]) == 'U' 
			&&
	             toupper(seq[2]) == 'C'
	                &&
	             toupper(seq[3]) == 'G'){
		  strcpy(name, "UNCG Tetra");
	    }else{
		  strcpy(name,"Tetra");
	    }
      }else{
	    name[0] = '\0';
      }
}

int helix_check_hairpin(struct helix* self, struct nucbp* nbp){
      const int loop_size = HAIRPIN_LOOP_MAX;
      self->is_hairpin = FALSE;
      self->pinsize = 0;
      char* ich = nbp[self->i[self->size-1]].chain;
      char* jch = nbp[self->j[self->size-1]].chain;

      /*
      int icif = nbp[self->i[self->size-1]].cifid;
      int jcif = nbp[self->j[self->size-1]].cifid;

      int iins = nbp[self->i[self->size-1]].ins;
      int jins = nbp[self->j[self->size-1]].ins;
      */
      if(strcmp(ich, jch) != 0) return FALSE;

      int pin_size = abs(self->i[self->size-1] - self->j[self->size-1]) - 1;
      if(pin_size > 0 && pin_size <= loop_size){
	    self->is_hairpin = TRUE;
	    self->pinsize = pin_size;
	    int index = self->i[self->size - 1] + 1; // Starting of loop.
	    if(index > self->j[self->size-1]){
		  index = self->j[self->size - 1] + 1; // Starting of Loop.
	    }
	    for(int i=0; i<self->pinsize; ++i){
		  self->loop_res[i] = nbp[index+i].resclass;
	    }
	    get_loop_name(self->loop_res, self->pinsize, self->loopname);

	    
	    return TRUE;
      }
      return FALSE;
}

void helix_compute(struct helix* self, int* hlxcount, struct nucbp* nbp, int nres){
      struct graph g;
      *hlxcount = 0;
      graph_init(&g, 100*nres, UNDIRECTED); 

      
      
      for(int i=0; i<nres; ++i){
	    
	    
	    
	    for(int k=0; k<nbp[i].numbp; ++k){
		  int j = nbp[i].oth_base_index[k];
		  if(graph_edge_index(&g, i, j) != -1) continue;
		  self[*hlxcount].size = 0;

		  
		  
		  int hval = helix_calc(self+ *hlxcount, nbp, nres, &g, i, j, INCR_DECR);
		  if(hval == FALSE){
			hval = helix_calc(self+ *hlxcount, nbp, nres, &g, i, j, DECR_INCR);
		  }

		  if(hval == FALSE){
			hval = helix_calc(self+ *hlxcount, nbp, nres, &g, i, j, INCR_INCR);
		  }
		  if(hval == FALSE){
			helix_calc(self+ *hlxcount, nbp, nres, &g, i, j, DECR_DECR);
		  }
		  helix_check_hairpin(self + *hlxcount, nbp);

		  *hlxcount = *hlxcount + 1;
		  
//		  if(self[*hlxcount].size >= 3){
//			*hlxcount = *hlxcount + 1;
//		  }
	    }

	    
      }


      
      graph_free(&g);
}
void helix_fprint(struct helix* self, struct nucbp* nbp, int count, FILE* fp){
      for(int i=0; i<count; ++i){
	    if(self[i].size >=3){
		  helix_printf(self+i, nbp, fp);
	    }
      }
}

void helix_free(struct helix* self){
      free ( self );
}



