/**
 * Copyright (C) 2015 Universidad Simón Bolívar
 *
 * @brief Implementation of a adjacent list graph
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <limits.h>
#include <errno.h>
#include <math.h>

#include "dlist.h"
#include "types.h"
#include "util.h"
#include "memory.h"
#include "graph_adj.h"


void init_graph_adj(struct graph_adj *g, int nodes)
{
     g->n_nodes = nodes;
     g->n_arcs = 0;
     g->adj_list = (struct arc **)xcalloc(nodes, sizeof(struct arc *));
     g->degree = (struct adj *)xcalloc(nodes, sizeof(struct adj));
     for (int i = 0; i < nodes; i++) {
	  g->degree[i].din = 0;
	  g->degree[i].dout = 0;
	  g->adj_list[i] = NULL;
     }
}

void add_arc(struct graph_adj *g, int id, int from, int to)
{
     struct arc *a;

     a = (struct arc *)xmalloc(sizeof(struct arc));
     a->id = id;
     a->to = to;
     a->next = NULL;

     a->next = g->adj_list[from];
     g->adj_list[from] = a;

     g->degree[from].dout++;
     g->degree[to].din++;
     g->n_arcs++;
}

void print_graph_adj(struct graph_adj *g)
{
     struct arc *current;
     
     printf("\nNum. of Nodes %d --- Num. of Arcs %d\n",g->n_nodes, g->n_arcs);
     for (int i = 0; i < g->n_nodes; i++) {
	  printf("Node %d - (in out) (%d, %d): ", i,  g->degree[i].din, g->degree[i].dout);
	  current = g->adj_list[i];
	  while (current != NULL) {
		    printf("(id %d f %d t %d) ",  current->id, i, current->to);
		    current = current->next;
	  }
	  printf("\n");
     }
}

void free_graph_adj(struct graph_adj *g)
{
     struct arc *current, *tmp;
     
     free(g->degree);
     for (int i = 0; i < g->n_nodes; i++) {
	  current = g->adj_list[i];
	  while (current != NULL) {
		    tmp = current;
		    current = current->next;
		    free(tmp);
	  }
     }
     free(g->adj_list);
     g->n_nodes = 0;
     g->n_arcs = 0;
}
