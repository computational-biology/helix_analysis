//
// Created by parthajit on 23/2/20.
//

#ifndef __spgraph_H__
#define __spgraph_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define EMPTY (-1)
#define MAX_DEGS (10)
enum gtype {DIRECTED =0, UNDIRECTED=1};


struct djset{
      int* leader;
      int* rank;
      int* link;
      int* cardinal;
      int* cycles;
      int size;
      int numset;

};

int intcmp (const void* a, const void* b) ;

void djset_sort(struct djset* self, int elem, int* array);
void djset_init(struct djset* self, int size);

void djset_free(struct djset* self);


int djset_find(struct djset* self, int v);


void djset_link(struct djset* self, int x, int y);

void djset_union(struct djset* self, int u, int v);

int djset_cycles(struct djset* self, int vertex);
int djset_composize(struct djset* self, int vertex);
int djset_next(struct djset* self, int vertex);




struct graph{
      int v;
      int* adj;
      int* deg;
      double* w;
      double memsize;
      enum gtype type;
};

void graph_free(struct graph* self);


void graph_init(struct graph* self, int vnum, enum gtype type);


int graph_edge_at(struct graph* self, int vertex, int index);


int graph_set_edge(struct graph* self, int i, int j);

double graph_get_wt(struct graph* self, int i, int index);
int graph_edge_index(struct graph* self, int i, int j);


void graph_set_wt(struct graph* self, int i, int index, double wt);


int graph_deg(struct graph* self, int vertex);


void graph_kruskal_component(struct graph* self, struct djset* set);


void graph_compo_isomorph_name(struct graph* self, struct djset* set, int vertex, char* name);
#endif //__spgraph_H__
