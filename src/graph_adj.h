/**
 * Copyright (C) 2015 Universidad Simón Bolívar
 *
 * @brief Implementation of a adjacent list graph
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#ifndef ___GRAPH_ADJ_H
#define ___GRAPH_ADJ_H

struct arc {
     int id;
     int to;
     struct arc *next;
};

struct graph_adj {
     int n_nodes;
     int n_arcs;
     struct adj {
	  int din;
	  int dout;
     } *degree;
     struct arc **adj_list;
};

void init_graph_adj(struct graph_adj *g, int nodes);

void add_arc(struct graph_adj *g, int id, int from, int to);

void print_graph_adj(struct graph_adj *g);

void free_graph_adj(struct graph_adj *g);

#endif
