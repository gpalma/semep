/**
 * Copyright (C) 2012, Universidad Simón Bolívar
 *
 * @brief Implementation of Abstract Data Type Graph how adjacents list
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#ifndef ___GRAPH_H
#define ___GRAPH_H

#include "dlist.h"

struct edge {
  long cost;
  long id;
  long from;
  long to;
};

struct edge_list {
  long n;
  struct edge item;
  struct list_head list;
};

struct graph {
  long n_nodes;
  long n_edges;
  long n_buckets;
  struct valency {
    long din;
    long dout;
  } *degree;
  struct edge_list *adj_list;

  /* Private variable */
  struct edge_list *_tmp;
};

typedef struct graph graph_t;

typedef int (*edge_cost_fn_t)(const struct edge *);

/**
 * Iterate over all the edges of a graph
 * @edge:    struct edge_list *
 * @graph:   struct graph *
 */
#define graph_for_each(edge, graph)					\
  for (long _i=0; (_i < (graph)->n_nodes); ++_i)                        \
    for (edge = list_entry( (&(((graph)->adj_list[_i]).list))->next, struct edge_list, list); \
         &((edge)->list) != &(((graph)->adj_list[_i]).list);            \
         edge = list_entry(((edge)->list).next, struct edge_list, list))
/**
 * Iterate over all the edges of a graph safe against removal of edge in the graph
 * @edge:    struct edge_list *
 * @graph:   struct graph *
 */
#define graph_for_each_safe(edge, graph)				\
  for (long _i=0; (_i < (graph)->n_buckets); ++_i)                      \
    for (edge = list_entry( (&(((graph)->adj_list[_i]).list))->next, struct edge_list, list), \
             (graph)->_tmp = list_entry(((edge)->list).next, struct edge_list, list); \
         &((edge)->list) != &(((graph)->adj_list[_i]).list);            \
         edge = (graph)->_tmp, (graph)->_tmp = list_entry((((graph)->_tmp)->list).next, struct edge_list, list))

/**
 * Iterate over the adjacent list of a graph
 * @edge:       struct edge_list *, Edge to return
 * @adj_head:   struct edge_list *, Head of adjacent list
 */
#define adj_for_each(edge, adj_head)					\
  for (edge = list_entry( (&((adj_head).list))->next, struct edge_list, list); \
       &((edge)->list) != &((adj_head).list);                           \
       edge = list_entry(((edge)->list).next, struct edge_list, list))

/**
 * Iterate over the adjacent list of a graph. Safe against
 *             removal of edge in the adjacent list
 * @edge:       struct edge_list *, Edge to return
 * @adj_head:   struct edge_list *, Head of adjacent list
 * @tmp         struct edge_list *, Temporal structure
 */
#define adj_for_each_safe(edge, adj_head, tmp)				\
  for (edge = list_entry( (&((adj_head).list))->next, struct edge_list, list), \
           tmp = list_entry(((edge)->list).next, struct edge_list, list); \
       &((edge)->list) != &((adj_head).list);                           \
       edge = tmp, tmp = list_entry(((tmp)->list).next, struct edge_list, list))

void set_edge(struct edge *e, long id, long from, long to, long c);

void init_graph(struct graph *g, long nodes);

struct edge_list *get_adjacent_list(struct graph *g, long from);

struct edge_list *new_edge_list(long id, long from, long to, long c);

void add_arc_to_graph(struct graph *g, long id, long from, long to, long c);

void print_graph(const struct graph *g);

void set_edges_cost(const struct graph *g);

void free_edge_item(struct edge_list *el);

void free_graph(struct graph *g);

void dfs_search(const struct graph *g, long s, long *d, long *f, long *pred);

long min_distance(const struct graph *g, long s, long t);

void all_pairs_shortest(const struct graph *g, long **dist);

void graph_inverse(const struct graph *orig, struct graph *inv);

void graph_undirect(const struct graph *orig, struct graph *ud);

long max_distance(const struct graph *g, long s, long t);

long min_path(const struct graph *g, long start, long goal);

struct long_list *topological_sort(const struct graph *g, long s);

long *calculate_depth(const struct graph *g);

void add_reprensentative_ancestor(struct graph *g);

long *calculate_depth_bfs(const struct graph *g);

bool detect_cycle(const struct graph *g);

bool *get_spanning_tree(const struct graph *g);

VEC(long) *get_euler_tour(const struct graph *g);

bool find_edge(struct graph *g, long v1, long v2);

#endif /* ___GRAPH_H */
