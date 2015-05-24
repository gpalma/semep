/**
 * Copyright (C) 2014, 2015 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>

#include "dlist.h"
#include "types.h"
#include "graph.h"
#include "memory.h"
#include "hash_map.h"
#include "util.h"
#include "metric.h"
#include "semEP.h"

#define COST       0
#define NOCOLOR   -1
#define NONODE    -1
#define INFTY      LONG_MAX
#define ROOT       0

/*
 * Macros used for binaries heaps implementation
 */
#define PARENT(i)     ((long)(((i) - 1) / 2))
#define LEFT(i)       (((i) * 2) + 1)
#define RIGHT(i)      (((i) * 2) + 2)

typedef struct color {
     long id;
     double sim_entity_1;
     double sim_entity_2;
     double sim_between;
     VEC(long) id_nodes;
     VEC(long) entities1;
     VEC(long) entities2;
} color_t;

typedef color_t *color_p;
DEFINE_VEC(color_p);

typedef struct node {
     long id;
     long pos1;
     long pos2;
     double sim;
     color_p cp;
} node_t;

typedef node_t *node_p;
DEFINE_VEC(node_p);

typedef struct prediction {
     long entity1;
     long entity2;
     double prob;
} prediction_t;

DEFINE_VEC(prediction_t);

typedef struct prediction_array {
     unsigned nr;
     VEC(prediction_t) *pred;
} prediction_array_t;

typedef struct item_entry {
     double value;
     struct hash_entry entry;
} item_entry_t;

typedef struct saturation {
     long node;
     long n_adj_colors;
     bool *color_used;
} saturation_t;

typedef struct pqueue {
     long size;
     saturation_t **heap;
} pqueue_t; 

static struct hash_map sim_map; /* Hash table with all similitudes calculated */
static double (*metricPtr)(const struct graph *g, long x, long y);

/*************************************
 *************************************
 **
 **  Utilities
 **
 ************************************
 ************************************/

static double get_similarity(const struct graph *g, long x, long y)
{
     int l;
     char *buf;
     item_entry_t *item;
     struct hash_entry *hentry;
     double sim;

     if (x > y) {
	  l = asprintf(&buf, "%ld-%ld", x, y);
     } else {
	  l = asprintf(&buf, "%ld-%ld", y, x);
     }

     if (l == -1)
	  fatal("Error in output directory");

     hentry = hmap_find_member(&sim_map, buf, l);
     if (hentry == NULL) {
	  sim = (*metricPtr)(g, x, y);
	  //sim = sim_dtax(g, x, y);
	  item = (item_entry_t *)xmalloc(sizeof(item_entry_t));
	  item->value = sim;
	  if (hmap_add(&sim_map, &item->entry, buf, l) != 0)
	       fatal("Error in similarity hash table");
     } else {
	  item = hash_entry(hentry, item_entry_t, entry);
	  sim = item->value;
     }
     free(buf);
     return sim;
}

static struct lpairs get_max_group_depth(const VEC(long) *v1)
{
     long max_depth = 0; /* ROOT depth */
     long max_node = ROOT;
     long i, n, node;
     const long *depth;
     struct lpairs p;

     depth = get_nodes_depth();
     n = VEC_SIZE(*v1);
     for (i = 0; i < n; i++) {
	  node = VEC_GET(*v1, i);
	  if (max_depth <  depth[node]) {
	       max_depth = depth[node];
	       max_node = node;
	  }
     }
     p.x = max_node;
     p.y = max_depth;

     return p;
}

/*************************************
 *************************************
 **
 **  Make the graph to coloring
 **
 ************************************
 ************************************/

static bool *nodes_in_the_term(const VEC(long) *terms, long max)
{
     unsigned i;
     bool *in_term;

     in_term = (bool *)xcalloc(max, sizeof(bool));
     memset(in_term, false, max*sizeof(bool));

     for (i = 0; i < VEC_SIZE(*terms); i++)
	  in_term[VEC_GET(*terms, i)] = true;

     return in_term;
}

static void compute_color_nodes(void *object, const VEC(long) *v1, const VEC(long) *v2,
				VEC(node_p) *color_nodes, double threshold,
				const bool* in_terms1, const bool* in_terms2, bool matrix)
{
     long i, j, n, m, x, y, cont;
     double sim;
     node_p new;
     double **M = NULL;
     graph_t *g;

     if (matrix) {
	  M = (double **)object;
     } else {
	  g = (graph_t *)object;
     }

     cont = 0;
     n = VEC_SIZE(*v1);
     m = VEC_SIZE(*v2);
     for (i = 0; i < n; i++) {
	  x = VEC_GET(*v1, i);
	  for (j = 0; j < m; j++) {
	       y = VEC_GET(*v2, j);
	       if (matrix)
		    sim = M[x][y];
	       else
		    sim = get_similarity(g,x,y);
	       if  (sim > threshold) {
		    if (!((in_terms2[x] && in_terms1[y]) && (x != y))) {
			 new = (node_p)xcalloc(1, sizeof(node_t));
			 new->id = cont;
			 new->pos1 = i;
			 new->pos2 = j;
			 new->sim = sim;
			 new->cp = NULL;
			 VEC_PUSH(node_p, *color_nodes, new);
			 cont++;
		    }
	       }
	  }
     }
}

static void build_graph_to_coloring_matrix(struct graph *gc, VEC(node_p) *vn,
					   double **M, const VEC(long) *v1, const VEC(long) *v2,
					   double threshold1, double threshold2)
{
     long i, j, n, a, b, c, d, cont;
     node_p x, y;

     n = VEC_SIZE(*vn);
     cont = 0;
     for (i = 0; i < n-1; i++) {
	  x = VEC_GET(*vn, i);
	  a = x->pos1;
	  b = x->pos2;
	  for (j = i+1; j < n; j++) {
	       y = VEC_GET(*vn, j);
	       c = y->pos1;
	       d = y->pos2;
	       assert(a <= c);
	       if ( (M[VEC_GET(*v1, a)][VEC_GET(*v1, c)] <= threshold1) ||
		    (M[VEC_GET(*v2, b)][VEC_GET(*v2, d)] <= threshold2)) {
		    add_arc_to_graph(gc, cont, i, j, COST);
		    add_arc_to_graph(gc, cont, j, i, COST);
		    cont++;
	       }
	  }
     }
     
}

static double **similarity_between_all(const struct graph *g, const VEC(long) *v)
{
     long i, j, n;
     double **sim_matrix;

     n = VEC_SIZE(*v);
     sim_matrix = double_matrix(0, n, 0, n);
     for (i = 0; i < n-1; i++) {
	  for (j = i+1; j < n; j++) {
	       sim_matrix[i][j] = get_similarity(g, VEC_GET(*v, i), VEC_GET(*v, j));
	       DEBUG("annt1 %ld annt2 %ld sim %.3f\n", VEC_GET(*v, i), VEC_GET(*v, j), sim_matrix[i][j]);
	  }
     }
     return sim_matrix;
}

static void build_graph_to_coloring_graph(struct graph *gc, VEC(node_p) *vn,
					  double **simM1, double **simM2,
					  double threshold1, double threshold2)
{
     long i, j, n, a, b, c, d, cont;
     node_p x, y;
     double sim2;

     if (simM1 == NULL)
	  fatal("Matrix 1 is NULL");
     if (simM2 == NULL)
	  fatal("Matrix 2 is NULL");
     n = VEC_SIZE(*vn);
     cont = 0;
     for (i = 0; i < n-1; i++) {
	  x = VEC_GET(*vn, i);
	  a = x->pos1;
	  b = x->pos2;
	  for (j = i+1; j < n; j++) {
	       y = VEC_GET(*vn, j);
	       c = y->pos1;
	       d = y->pos2;
	       assert(a <= c);
	       if (b < d)
		    sim2 = simM2[b][d];
	       else
		    sim2 = simM2[d][b];
	       if ( (simM1[a][c] <= threshold1) || (sim2 <= threshold2)) {
		    add_arc_to_graph(gc, cont, i, j, COST);
		    add_arc_to_graph(gc, cont, j, i, COST);
		    cont++;
	       }
	  }
     }
}

/*************************************
 *************************************
 **
 **  Coloration Solver
 **
 ************************************
 ************************************/

static inline int compare_saturation(const graph_t *g, const saturation_t *a, const saturation_t *b)
{
     int r;

     if (a->n_adj_colors > b->n_adj_colors) { 
	  r = 1;
     } else if (a->n_adj_colors < b->n_adj_colors) {
	  r = -1;
     } else {
	  if (g->degree[a->node].din > g->degree[b->node].din) {
	       r = 1;
	  } else if (g->degree[a->node].din < g->degree[b->node].din) { 
	       r = -1;
	  } else {
	       /* Alphabetical order */
	       if (a->node < b->node)
		    r = 1;
	       else
		    r = -1;
	  }
     }
     return r;
}

static void free_saturation_node(saturation_t *node)
{
     if (node) {
	  if (node->color_used == NULL)
	       fatal("Error in delete table of used color");
	  free(node->color_used);
	  free(node);
     }
}

static inline void pq_init(pqueue_t *pq)
{
     pq->size = 0;
     pq->heap = NULL;
}

static inline void pq_delete(pqueue_t *pq)
{
     long i;
     
     for(i = 0; i < pq->size; i++)
	  free_saturation_node(pq->heap[i]);
     free(pq->heap);
}

static inline void pq_insert(const struct graph *g, pqueue_t *pq, saturation_t *node)
{
     long i, p;
     saturation_t **tmp;

     tmp = xrealloc(pq->heap, (pq->size+1)*sizeof(saturation_t *));
     pq->heap = tmp;
     pq->heap[pq->size] = node;
     i = pq->size;
     p = PARENT(i);
     while((i > 0) &&  (compare_saturation(g, pq->heap[p], pq->heap[i]) < 0)){
	  SWAP(pq->heap[p], pq->heap[i]);
	  i = p;
	  p = PARENT(i);
     }
     pq->size++;
}

static long extract_max(const struct graph *g, pqueue_t *pq, saturation_t **node)
{
     long i, j, l, r;
     saturation_t *aux;
     saturation_t **tmp;

     if(pq->size == 0)
	  return -1;

     *node = pq->heap[0];
     aux =  pq->heap[pq->size-1];
     if((pq->size - 1) > 0){
	  tmp = xrealloc(pq->heap, (pq->size-1)*sizeof(saturation_t *));
	  pq->heap = tmp;
	  pq->size--;
     } else {
	  free(pq->heap);
	  pq->heap = NULL;
	  pq->size = 0;
	  return 0;
     }
     pq->heap[0] = aux;
     i = 0;
     while (true) {
	  l = LEFT(i);
	  r = RIGHT(i);
	  if((l < pq->size) && (compare_saturation(g, pq->heap[l], pq->heap[i]) > 0))
	       j = l;
	  else
	       j = i;

	  if((r < pq->size) && (compare_saturation(g, pq->heap[r], pq->heap[j]) > 0))
	       j = r;

	  if( j == i ) {
	       break;
	  } else {
	       SWAP(pq->heap[j], pq->heap[i]);
	       i = j;
	  }
     }
     return 0;
}

static long increase_key(const struct graph *g, pqueue_t *pq, long node, long color)
{
     long i, p, pos;

     if (pq->size == 0)
	  return -1;
     pos = -1;
     for (i = 0; i < pq->size; i++)
	  if (pq->heap[i]->node == node)
	       pos = i;
     if (pos == -1)
	  return -2;

     if (!pq->heap[pos]->color_used[color]) {
	  pq->heap[pos]->n_adj_colors++;
	  pq->heap[pos]->color_used[color] = true;
     } else {
	  return 0;
     }

     i = pos;
     p = PARENT(i);
     while((i > 0) && (compare_saturation(g, pq->heap[p], pq->heap[i]) < 0)){
	  SWAP(pq->heap[p], pq->heap[i]);
	  i = p;
	  p = PARENT(i);
     }
     return 0;
}

static void init_saturation_pq(const struct graph *g, pqueue_t *pq)
{
     long i, n;
     saturation_t *ns;
     size_t alloc;

     n = g->n_nodes;
     for (i = 0; i < n; i++) {
	  ns = (saturation_t *)xmalloc(sizeof(saturation_t));
	  ns->node = i;
	  ns->n_adj_colors = 0;
	  alloc = n*sizeof(bool)+1;
	  ns->color_used = xmalloc(alloc);
	  memset(ns->color_used, false, alloc);
	  pq_insert(g, pq, ns);
     }
}

static inline long get_color(const VEC(node_p) *node_color, long pos)
{
     long color;
     node_p node_ptr;

     if (pos < 0)
	  fatal("Invalid position");
     if (node_color == NULL)
	  fatal("Invalid pointer to color");
     node_ptr = VEC_GET(*node_color, pos);
     if (node_ptr->cp == NULL) {
	  color = NOCOLOR;
     } else {
	  color = node_ptr->cp->id;
     }
     return color;
}

static bool *get_ady_used_color(const struct graph *g, const VEC(node_p) *node_color, long node)
{
     struct edge_list *tmp;
     bool *color_used;
     long color;
     size_t alloc;

     alloc = g->n_nodes*sizeof(bool)+1;
     color_used = xmalloc(alloc);
     memset(color_used, false, alloc);
     adj_for_each(tmp, g->adj_list[node]) {
	  color = get_color(node_color, tmp->item.to);
	  assert(color < g->n_nodes);
	  if (color != NOCOLOR) {
	       color_used[color] = true;
	  }
     }
     return color_used;
}

static void update_saturation_degree(const struct graph *g, pqueue_t *pq, 
				     long node, const VEC(node_p) *node_color)
{
     struct edge_list *tmp;
     int r;
     long color;
     
     color = get_color(node_color, node);
     assert(color != NOCOLOR);
     adj_for_each(tmp, g->adj_list[node]) {
	  r = increase_key(g, pq, tmp->item.to, color);
	  if (r == -1)
	       fatal("Error in update saturation degree\n");
     }
}

static long greatest_saturation_node(const graph_t *g, pqueue_t *pq, const VEC(node_p) *node_color)
{
     long r, color, node;
     saturation_t *ns;

     node = NONODE;
     ns = NULL;
     color = INFTY;
     r = extract_max(g, pq, &ns);
     if (r == 1)
	  fatal("No node without color");
     if (ns) {
	  color = get_color(node_color, ns->node);
	  node = ns->node;
     } else {
	  fatal("Error in get the greatest saturation node");
     }
     if (color != NOCOLOR)
	  fatal("Error in node to coloring");
     
#ifdef PRGDEBUG
     printf("Node %ld degree %ld din %ld\n",node, ns->n_adj_colors, g->degree[node].din);
#endif
     free_saturation_node(ns); 
     return node;
}

static void get_free_colors(const graph_t *g, const VEC(node_p) *solution, long node,
			    VEC(long) *available_colors, VEC(color_p) *partitions)
{
     long i, cn, n;
     bool *color_used;
     long color_no_used;
     color_p ctmp;

     n = VEC_SIZE(*partitions);
     color_no_used = n+1;
     color_used = NULL;
     assert(g->n_nodes > node);
     assert(g->n_nodes >= color_no_used);
     color_used = get_ady_used_color(g, solution, node);
     cn = get_color(solution, node);
     if ((cn != NOCOLOR) && (!color_used[cn])) /* any adjacent vertex has the same color */
	  fatal("a adjacent node are using the same color");
     for (i = 0; i < n; i++) {
	  ctmp = VEC_GET(*partitions, i);
	  assert(ctmp->id == i);
	  if (!color_used[i]) {
	       VEC_PUSH_FAST(*available_colors, i);
	  }
     }
     if (VEC_SIZE(*available_colors) == 0)
	  VEC_PUSH_FAST(*available_colors, n);
     free(color_used);
}

static inline color_p new_color(long id, double sim_e1, double sim_e2, double sim_bt,
				long node, long e1, long e2)
{
     color_p new;

     new = (color_p)xcalloc(1, sizeof(color_t));
     new->id = id;
     new->sim_entity_1 = sim_e1;
     new->sim_entity_2 = sim_e2;
     new->sim_between = sim_bt;
     VEC_INIT(long, new->id_nodes);
     VEC_PUSH(long, new->id_nodes, node);
     VEC_INIT(long, new->entities1);
     VEC_PUSH(long, new->entities1, e1);
     VEC_INIT(long, new->entities2);
     VEC_PUSH(long, new->entities2, e2);

     return new;
}

static void free_color(color_p c)
{
     VEC_DESTROY(c->id_nodes);
     VEC_DESTROY(c->entities1);
     VEC_DESTROY(c->entities2);
     free(c);
}

static double get_partition_simililarity(color_p c)
{
     long n, m;
     double bpe, annt1, annt2;

     n = VEC_SIZE(c->entities1);
     m = VEC_SIZE(c->entities2);
     bpe = (double)(c->sim_between/VEC_SIZE(c->id_nodes));
     if (n > 1) {
	  annt1 = (double)(c->sim_entity_1/((n*n-n)/2));
     } else {
	  annt1 = c->sim_entity_1;
	  assert(annt1 == 0.0);
     }
     if (m > 1) {
	  annt2 = (double)(c->sim_entity_2/((m*m-m)/2));
     } else {
	  annt2 = c->sim_entity_2;
	  assert(annt2 == 0.0);
     }
     return (bpe + annt1 + annt2);
}

static double density_total_add(void *object, VEC(color_p) *partitions, node_p new_node,
				const VEC(long) *v1, const VEC(long) *v2, long color, bool matrix)
{
     long i, n, m, e1, e2, ep, n_colors;
     color_p cptr;
     double annt1_nc, annt2_nc, bpe_nc, dt, sim, bpe, annt1, annt2;
     graph_t *g;
     double **M;

     if (matrix)
	  M = (double **)object;
     else
	  g = (graph_t *)object;

     dt = 0.0;
     cptr = VEC_GET(*partitions, color);
     assert(cptr->id == color);
     n = VEC_SIZE(cptr->entities1);
     e1 = VEC_GET(*v1, new_node->pos1);
     annt1_nc = 0.0;
     for (i = 0; i < n; i++) {
	  ep = VEC_GET(cptr->entities1, i);
	  if (matrix)
	       annt1_nc += M[ep][e1];
	  else
	       annt1_nc += get_similarity(g, ep, e1);
     }

     n = VEC_SIZE(cptr->entities2);
     e2 = VEC_GET(*v2, new_node->pos2);
     annt2_nc = 0.0;
     for (i = 0; i < n; i++) {
	  ep = VEC_GET(cptr->entities2, i);
	  if (matrix)
	       annt2_nc += M[ep][e2];
	  else
	       annt2_nc += get_similarity(g, ep, e2);
     }
     bpe_nc = cptr->sim_between;
     n_colors = VEC_SIZE(*partitions);
     for (i = 0; i < n_colors; i++) {
	  cptr = VEC_GET(*partitions, color);
	  assert(cptr->id == color);
	  if (i == color) {
	       n = VEC_SIZE(cptr->entities1) + 1;
	       m = VEC_SIZE(cptr->entities2) + 1;
	       bpe = (double)((cptr->sim_between+bpe_nc)/(VEC_SIZE(cptr->id_nodes)+1));
	       if (n > 1) {
		    annt1 = (double)((cptr->sim_entity_1 + annt1_nc)/((n*n-n)/2));
	       } else {
		    annt1 = annt1_nc;
		    assert(cptr->sim_entity_1 == 0.0);
	       }
	       if (m > 1) {
		    annt2 = (double)((cptr->sim_entity_2 + annt2_nc)/((m*m-m)/2));
	       } else {
		    annt2 = annt2_nc;
		    assert(cptr->sim_entity_2 == 0.0);
	       }
	       sim = bpe + annt1 + annt2;
	  } else {
	       sim = get_partition_simililarity(cptr);
	  }
	  dt += 1.0 - (sim/3.0);
     }
     assert(dt/n_colors <= 1.0);
     return dt;
}

static double density_total(VEC(color_p) *partitions)
{
     long i, n_colors;
     double dt;
     color_p c;

     dt = 0.0;
     n_colors = VEC_SIZE(*partitions);
     for (i = 0; i < n_colors; i++) {
	  c = VEC_GET(*partitions, i);
	  dt += 1.0 - (get_partition_simililarity(c)/3.0);
     }
     assert(dt/n_colors <= 1.0);
     return dt;
}

static long get_color_of_minimum_density(void *object, const VEC(node_p) *vn,
					 VEC(color_p) *partitions,
					 long new_node, const VEC(long) *free_colors,
					 const VEC(long) *v1, const VEC(long) *v2, bool matrix)
{
     long i, n, best_color, n_colors, new_color;
     double nc_best, nc, current_density;
     node_p nptr;

     n = VEC_SIZE(*free_colors);
     best_color = -1;
     nc_best = INFTY;
     n_colors = VEC_SIZE(*partitions);
     current_density = density_total(partitions);
     nptr = VEC_GET(*vn, new_node);
     for (i = 0; i < n; i++) {
	  new_color = VEC_GET(*free_colors, i);
	  assert((new_color >= 0) && (new_color <= n_colors));
	  DEBUG("Color to evaluate %ld ", new_color);
	  if (n_colors == new_color) {
	       /* mew color */
	       nc = (double)((current_density + (1.0 - nptr->sim/3.0))/(n_colors+1));
	  } else {
	       nc = (double)(density_total_add(object, partitions, nptr, v1, v2, new_color, matrix)/n_colors);
	  }
	  if (nc_best > nc) {
	       nc_best = nc;
	       best_color = new_color;
	  }
     }
     if (n_colors < best_color) {
	  n = n_colors+1;
     } else {
	  n = n_colors;
     }
     if ((double)(nc_best/n) > 1.0) {
	  fatal("computed a density not valid");
     }
     return best_color;
}

static double compute_pairwise_sim(void *object, const VEC(long) *v, bool matrix)
{
     long i, j, n, m;
     double sim;
     double **M;
     graph_t *g;

     if (matrix)
	  M = (double **)object;
     else
	  g = (graph_t *)object;
     sim = 0.0;
     m = VEC_SIZE(*v);
     n = m - 1;
     for (i = 0; i < n; i++)  {
	  for (j = i+1; j < m; j++) {
	       if (matrix)
		    sim += M[VEC_GET(*v, i)][VEC_GET(*v, j)];
	       else
		    sim += get_similarity(g, VEC_GET(*v, i), VEC_GET(*v, j));
	  }
     }
     return sim;
}

/*
 * Return true if the entity was added
 */
static bool add_if_not_contains(VEC(long) *v, long entity)
{
     unsigned i;

     for (i = 0; i < VEC_SIZE(*v); i++) {
	  if (VEC_GET(*v, i) == entity)
	       return false;
     }
     VEC_PUSH(long, *v, entity);
     return true;
}

static void set_colors(void *object, VEC(color_p) *partitions, long color, node_p nptr,
		       const VEC(long) *v1, const VEC(long) *v2, bool matrix)
{
     long n_colors, e1, e2;
     color_p cptr;
     bool added;

     if (color < 0)
	  fatal("Invalid color");
     if (partitions == NULL)
	  fatal("Invalid pointer to partition");
     n_colors = VEC_SIZE(*partitions);
     if (n_colors == color) {
	  e1 = VEC_GET(*v1, nptr->pos1);
	  e2 = VEC_GET(*v2, nptr->pos2);
	  cptr = new_color(color, 0.0, 0.0, nptr->sim, nptr->id, e1, e2);
	  VEC_PUSH(color_p, *partitions, cptr);
	  nptr->cp = cptr;
     } else if (n_colors > color) {
	  cptr = VEC_GET(*partitions, color);
	  VEC_PUSH(long, cptr->id_nodes, nptr->id);

	  e1 = VEC_GET(*v1, nptr->pos1);
	  added = add_if_not_contains(&cptr->entities1, e1);
	  if (added)
	       cptr->sim_entity_1 = compute_pairwise_sim(object, &cptr->entities1, matrix);

	  e2 = VEC_GET(*v2, nptr->pos2);
	  added = add_if_not_contains(&cptr->entities2, e2);
	  if (added)
	       cptr->sim_entity_2 = compute_pairwise_sim(object, &cptr->entities2, matrix);

	  cptr->sim_between += nptr->sim;
	  nptr->cp = cptr;
     } else {
	  fatal("Bad color");
     }
}

static void coloring(void *object, const graph_t *g, const VEC(node_p) *nodes,
		     VEC(color_p) *partitions, const VEC(long) *v1, const VEC(long) *v2, bool matrix)
{
     color_p cptr;
     node_p nptr;
     long colored_nodes, new_node, e1, e2, color;
     VEC(long) free_colors;
     pqueue_t pq_saturation;

     colored_nodes = 0;
     VEC_INIT_N(long, free_colors, g->n_nodes);
     pq_init(&pq_saturation);
     init_saturation_pq(g, &pq_saturation);

     if (VEC_SIZE(*partitions) == 0)  {
	  new_node = greatest_saturation_node(g, &pq_saturation, nodes);
	  if (new_node == NONODE)
	       fatal("Error getting the greatest saturation node");
	  nptr = VEC_GET(*nodes, new_node);
	  e1 = VEC_GET(*v1, nptr->pos1);
	  e2 = VEC_GET(*v2, nptr->pos2);
	  assert(new_node == nptr->id);
	  cptr = new_color(0, 0.0, 0.0, nptr->sim, nptr->id, e1, e2);
	  VEC_PUSH(color_p, *partitions, cptr);
	  nptr->cp = cptr;
	  colored_nodes++;
	  if (pq_saturation.size != 0)
	       update_saturation_degree(g, &pq_saturation, new_node, nodes);
     }
     
     while (colored_nodes < g->n_nodes) {
	  new_node = greatest_saturation_node(g, &pq_saturation, nodes);
	  if (new_node == NONODE)
	       fatal("Error getting the greatest saturation node");
	  nptr = VEC_GET(*nodes, new_node);
	  assert(new_node == nptr->id);
	  VEC_CLEAR(free_colors);
	  get_free_colors(g, nodes, new_node, &free_colors, partitions);
#ifdef PRGDEBUG
	  print_vec_long(&free_colors);
#endif
	  color = get_color_of_minimum_density(object, nodes, partitions, new_node, &free_colors, v1, v2, matrix);
	  set_colors(object, partitions, color, nptr, v1, v2, matrix);
	  colored_nodes++;
	  if (pq_saturation.size != 0)
	       update_saturation_degree(g, &pq_saturation, new_node, nodes);
     }
     assert(pq_saturation.size == 0);
     pq_delete(&pq_saturation);
     VEC_DESTROY(free_colors);
}

#ifdef PRGDEBUG
static void print_coloring(VEC(color_p) *partitions)
{
     color_p c;
     for (unsigned i = 0; i < VEC_SIZE(*partitions); i++) {
	  c = VEC_GET(*partitions, i);
	  assert(c->id == i);
	  printf("%ld %.4f %.4f %.4f %zu %zu %zu\n", c->id, c->sim_entity_1, c->sim_entity_2, c->sim_between, VEC_SIZE(c->id_nodes), VEC_SIZE(c->entities1), VEC_SIZE(c->entities2));
     }
}
#endif

static double get_density_average(VEC(color_p) *partitions)
{
     return density_total(partitions) / VEC_SIZE(*partitions);
}

/*************************************
 *************************************
 **
 **  Predited Links
 **
 ************************************
 ************************************/

static void init_prediction_array(prediction_array_t *a, unsigned n)
{
     unsigned i;

     a->nr = n;
     a->pred = (VEC(prediction_t) *)xcalloc(n, sizeof(VEC(prediction_t)));
     for (i = 0; i < n; i++)
	  VEC_INIT(prediction_t, a->pred[i]);
}

static void free_prediction_array(prediction_array_t *a)
{
     unsigned i;

     for (i = 0; i < a->nr; i++)
	  VEC_DESTROY(a->pred[i]);
     free(a->pred);
}

static void free_hash_map_item(struct hash_map *hmap)
{
     item_entry_t *item;
     struct hash_entry *hentry;
     struct hlist_node *n;

     hmap_for_each_safe(hentry, n, hmap) {
	  item = hash_entry(hentry, item_entry_t, entry);
	  hmap_delete(hmap, hentry);
	  free(item);
     }
     hmap_destroy(hmap);
}

static inline double get_cluster_probability(long nodes1, long nodes2, long n_edges)
{
     double n_nodes = nodes1*nodes2;
     if (n_nodes == 0)
	  fatal("Zero Division");
     return n_edges / n_nodes;
}

static void get_predicted_links(prediction_array_t *cluster_pred,
				VEC(node_p) *color_nodes, VEC(color_p) *partitions,
				const VEC(long) *v1, const VEC(long) *v2)
{
     unsigned i, j, k, n_links, m, n, n_clusters;
     struct hash_map edges_obs[VEC_SIZE(*partitions)];
     struct hash_entry *hentry;
     char *buf;
     item_entry_t *item;
     int l;
     node_p edge;
     color_p cluster;
     long node1, node2;
     prediction_t ptemp;
     long n_edges[VEC_SIZE(*partitions)];
     double prob;

     n_clusters = VEC_SIZE(*partitions);
     n_links = VEC_SIZE(*color_nodes);
     for (i = 0; i < n_clusters; i++) {
	  hmap_create(&edges_obs[i], n_links);
	  n_edges[i] = 0;
     }

     for (i = 0; i < n_links; i++) {
	  edge = VEC_GET(*color_nodes, i);
	  l = asprintf(&buf, "%ld-%ld", VEC_GET(*v1, edge->pos1), VEC_GET(*v2, edge->pos2));
	  if (l == -1)
	       fatal("Error in edge key creation");
	  assert((unsigned)l == strlen(buf));
	  item = (item_entry_t *)xmalloc(sizeof(item_entry_t));
	  assert(edge->cp->id <= n_clusters); 
	  n_edges[edge->cp->id]++;
	  if (hmap_add_if_not_member(&edges_obs[edge->cp->id], &item->entry, buf, l) != NULL)
	       fatal("Error, repeat edge in bipartite graph\n");
	  free(buf);
     }

     for (i = 0; i < n_clusters; i++) {
	  cluster = VEC_GET(*partitions, i);
	  n = VEC_SIZE(cluster->entities1);
	  m = VEC_SIZE(cluster->entities2);
	  prob = get_cluster_probability(n, m, n_edges[i]);
	  assert((prob >= 0.0) && (prob <= 1.0));
	  for (j = 0; j < n; j++) {
	       node1 = VEC_GET(cluster->entities1, j);
	       for (k = 0; k < m; k++) {
		    node2 = VEC_GET(cluster->entities2, k);
		    l = asprintf(&buf, "%ld-%ld",node1, node2);
		    if (l == -1)
			 fatal("Error in make a new link");
		    assert((unsigned)l == strlen(buf));
		    hentry = hmap_find_member(&edges_obs[i], buf, l);
		    if (hentry == NULL) {
			 ptemp.entity1 = node1;
			 ptemp.entity2 = node2;
			 ptemp.prob = prob;
			 VEC_PUSH(prediction_t, cluster_pred->pred[i], ptemp);
		    }
		    free(buf);
	       }
	  }
     }
     for (i = 0; i < n_clusters; i++)
	  free_hash_map_item(&edges_obs[i]);
}

static void print_predicted_links(prediction_array_t *cluster_pred,
				  double threshold_E1, double threshold_E2,
				  const char name1[], const char name2[], char **desc)
{
     unsigned i, j, n;
     prediction_t ptemp;
     FILE *f;
     char *output;
     long cont = 0;

     if (asprintf(&output, "%s-%s-%.4f-%.4f-Pred.txt",
		  name1, name2, threshold_E1, threshold_E2) == -1)
	  fatal("Error in prediction file");
     f = fopen(output, "w");
     free(output);
     if (!f)
	  fatal("No descriptor file specified, abort\n");

     for (i = 0; i < cluster_pred->nr; i++) {
	  n = VEC_SIZE(cluster_pred->pred[i]);
	  if (n > 0) {
	       fprintf(f, "Cluster\t%u\n", i);
	       cont += n;
	  }
	  for (j = 0; j < n; j++) {
	       ptemp = VEC_GET(cluster_pred->pred[i], j);
	       fprintf(f, "%s\t%s\t%.4f\n", desc[ptemp.entity1], desc[ptemp.entity2], ptemp.prob);
	  }
     }
     printf("Number of links predicted: %ld\n", cont);
     fclose(f);
}

/*************************************
 *************************************
 **
 **  Output files
 **
 ************************************
 ************************************/

static unsigned print_singleton(FILE *f, bool *in_v, const VEC(long) *v, char **desc)
{
     unsigned i, n;
     
     n = 0;
     for (i = 0; i < VEC_SIZE(*v); i++) {
	  if (!in_v[i]) {
	       fprintf(f ,"%s\n", desc[VEC_GET(*v, i)]);
	       n++;
	  }
     }
     return n;
}

static char *print_output_files(VEC(node_p) *color_nodes, VEC(color_p) *partitions,
				const VEC(long) *v1, const VEC(long) *v2,
				double threshold_E1, double threshold_E2,
				const char name1[], const char name2[], char **desc)
{
     FILE *f;
     unsigned i, j, n, m;
     char *output1, *output2, *message;
     struct stat st;
     node_p edge;
     color_p cluster;
     long id_node;
     bool *in_v1, *in_v2; 

     if (asprintf(&output1, "%s-%s-%.4f-%.4f-Subgr", name1, name2, threshold_E1, threshold_E2) == -1)
	  fatal("Error in output directory");
     if (stat(output1, &st) == -1)
	  mkdir(output1, 0700);

     if (asprintf(&output2, "%s-%s-%.4f-%.4f.txt", name1, name2, threshold_E1, threshold_E2) == -1)
	  fatal("Error in output file");

     if (asprintf(&message, "Output files:\n\t%s\n\t%s\n", output1, output2) == -1)
	  fatal("Error in output message");
     
     f = fopen(output2, "w");
     free(output2);

     in_v1 = (bool *)xcalloc(VEC_SIZE(*v1), sizeof(bool));
     memset(in_v1, false, VEC_SIZE(*v1)*sizeof(bool));
     in_v2 = (bool *)xcalloc(VEC_SIZE(*v2), sizeof(bool));
     memset(in_v2, false, VEC_SIZE(*v2)*sizeof(bool));

     if (!f)
	  fatal("No descriptor file specified, abort\n");
     for (i = 0; i < VEC_SIZE(*color_nodes); i++) {
	  edge = VEC_GET(*color_nodes, i);
	  /*assert(edge->id == i);*/
	  fprintf(f ,"%s\t%s\t%.4f\n", desc[VEC_GET(*v1, edge->pos1)],
		  desc[VEC_GET(*v2, edge->pos2)], edge->sim);
	  in_v1[edge->pos1] = true;
	  in_v2[edge->pos2] = true;
     }
     n = print_singleton(f, in_v1, v1, desc);
     m = print_singleton(f, in_v2, v2, desc);
     fclose(f);

     printf("Number of partitions: %zu\n", VEC_SIZE(*partitions));
     for (i = 0; i< VEC_SIZE(*partitions); i++) {
	  cluster = VEC_GET(*partitions, i);
	  /* assert(cluster->id == i); */
	  if (asprintf(&output2, "%s/%s-%s-%u-%.4f-%.4f.txt",
		       output1, name1, name2, i, threshold_E1, threshold_E2) == -1)
	       fatal("Error in cluster file");
          f = fopen(output2, "w");
	  if (!f)
	       fatal("No descriptor file specified, abort\n");
	  for (j = 0; j < VEC_SIZE(cluster->id_nodes); j++) {
	       id_node = VEC_GET(cluster->id_nodes, j);
	       edge = VEC_GET(*color_nodes, id_node);
	       assert(cluster->id == edge->cp->id);
	       fprintf(f ,"%s\t%s\t%.4f\n", desc[VEC_GET(*v1, edge->pos1)],
		       desc[VEC_GET(*v2, edge->pos2)], edge->sim);
	  }
	  fclose(f);
	  free(output2);
     }

     /* Print a file with the singleton nodes */
     if ((n+m) > 0) {  
	  if (asprintf(&output2, "%s/%s-%s-%.4f-%.4f-singleton.txt", 
		       output1, name1, name2, threshold_E1, threshold_E2) == -1)
	       fatal("Error in singleton file");
	  f = fopen(output2, "w");
	  if (!f)
	       fatal("No descriptor file specified, abort\n");
	  print_singleton(f, in_v1, v1, desc);
	  print_singleton(f, in_v2, v2, desc);
	  fclose(f);
	  free(output2);
     }
     free(in_v1);
     free(in_v2);
     free(output1);

     return message;
}

/*************************************
 *************************************
 **
 **  Partition detection solver
 **
 ************************************
 ************************************/

double annotation_partition(void *object, long n_nodes,
			    const VEC(long) *v1, const VEC(long) *v2,
			    double threshold_E1, double threshold_E2,
			    double threshold_bt, const char name1[], const char name2[],
			    char **desc, bool prediction, bool matrix, enum measure d)
{
     graph_t gc;
     bool *in_term1, *in_term2;
     VEC(node_p) color_nodes;
     VEC(color_p) partitions;
     unsigned i, n;
     prediction_array_t cluster_pred;
     double **sim_mat, **simM1, **simM2;
     graph_t *g;
     clock_t ti, tf;
     struct lpairs nd1, nd2;
     long max_depth;
     char *message;
     double sim;
     
     ti = clock();
     sim_mat = NULL;
     simM1 = NULL;
     simM2 = NULL;
     sim = 0.0;

     if (matrix) {
	  sim_mat = (double **)object;
     } else {
	  n = 2*VEC_SIZE(*v1)*VEC_SIZE(*v2);
	  hmap_create(&sim_map, n);
	  g = (graph_t *)object;
	  init_metric_data(g);
	  if (d == DTAX) {
	       metricPtr = &sim_dtax;
	  } else if (d == DPS) {
	       metricPtr = &sim_dps;
	  } else {
	       nd1 = get_max_group_depth(v1);
	       nd2 = get_max_group_depth(v2);
	       max_depth = MAX(nd1.y, nd2.y);
	       if (max_depth == nd1.y)
		    printf("Deepest node in the annotations is %s with depth %ld\n", desc[nd1.x], nd1.y);
	       else 
		    printf("Deepest node in the annotations is %s with depth %ld\n", desc[nd2.x], nd2.y);
	       set_max_depth(max_depth);
	       metricPtr = &sim_str;
	  }
	  simM1 = similarity_between_all(g, v1);
	  simM2 = similarity_between_all(g, v2);
     }
     VEC_INIT(color_p, partitions);
     in_term1 = nodes_in_the_term(v1, n_nodes);
     in_term2 = nodes_in_the_term(v2, n_nodes);
     VEC_INIT(node_p, color_nodes);
     compute_color_nodes(object, v1, v2, &color_nodes, threshold_bt, in_term1, in_term2, matrix);
     free(in_term1);
     free(in_term2);
     init_graph(&gc, VEC_SIZE(color_nodes));
     if (matrix) {
	  build_graph_to_coloring_matrix(&gc, &color_nodes, sim_mat, v1, v2, threshold_E1, threshold_E2);
     } else {
	  build_graph_to_coloring_graph(&gc, &color_nodes, simM1, simM2, threshold_E1, threshold_E2);
	  free_double_matrix(simM1, 0, 0);
	  free_double_matrix(simM2, 0, 0);
     }
     tf = clock();
     printf("Bipartite Graph data - Nodes A: %zu; Nodes B: %zu; Edges: %zu\n",
	    VEC_SIZE(*v1), VEC_SIZE(*v2), VEC_SIZE(color_nodes));
     printf("Time to build the graph to coloring: %.3f secs\n", (double)(tf-ti)/CLOCKS_PER_SEC);
#ifdef PRGDEBUG
     print_graph(&gc);
#endif
     printf("Graph to Colorinng - Nodes: %ld; Edges: %ld\n", gc.n_nodes, gc.n_edges/2);
     ti = clock();
     if (gc.n_nodes != 0) {
	  coloring(object, &gc, &color_nodes, &partitions, v1, v2, matrix);
	  sim = get_density_average(&partitions);
#ifdef PRGDEBUG
	  print_coloring(&partitions);
#endif
     } else {
	  fatal("Graph to coloring has no nodes");
     }
     tf = clock();
     printf("Coloring solver time %.3f secs\n", (double)(tf-ti)/CLOCKS_PER_SEC);
     message = print_output_files(&color_nodes, &partitions, v1, v2, threshold_E1, threshold_E2, name1, name2, desc);
     printf("%s", message);
     if (prediction) {
	  init_prediction_array(&cluster_pred, VEC_SIZE(partitions));
	  get_predicted_links(&cluster_pred, &color_nodes, &partitions, v1, v2);
	  print_predicted_links(&cluster_pred, threshold_E1, threshold_E2, name1, name2, desc);
	  free_prediction_array(&cluster_pred);
     }

     /*
      * free all memory allocated
      */
     free(message);
     if (!matrix) {
	  free_metric();
	  free_hash_map_item(&sim_map);
     }
     for (i = 0; i < VEC_SIZE(partitions); i++)
	  free_color(VEC_GET(partitions, i));
     VEC_DESTROY(partitions);
     for (i = 0; i < VEC_SIZE(color_nodes); i++)
	  free(VEC_GET(color_nodes, i));
     VEC_DESTROY(color_nodes);
     free_graph(&gc);
     
     return sim;
}
