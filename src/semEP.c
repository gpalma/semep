/**
 * Copyright (C) 2013, 2014, 2015 Universidad Simón Bolívar
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

struct color {
     long id;
     double sim_entity_1;
     double sim_entity_2;
     double sim_between;
     struct long_array id_nodes;
     struct long_array entities1;
     struct long_array entities2;
};

struct color_ptr_array{
     unsigned nr;
     unsigned alloc;
     struct color **data;
};

struct node {
     long id;
     long pos1;
     long pos2;
     double sim;
     struct color *cp;
};

struct node_ptr_array{
     unsigned nr;
     unsigned alloc;
     struct node **data;
};

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

typedef struct prediction {
     long entity1;
     long entity2;
     double prob;
} prediction_t;

typedef struct prediction_array {
     unsigned nr;
     unsigned alloc;
     prediction_t *data;
} prediction_array_t;

typedef struct prediction_partition {
     unsigned nr;
     prediction_array_t *pred; 
} prediction_partition_t;

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
     return (*metricPtr)(g, x, y);
}

static struct lpairs get_max_group_depth(const struct long_array *v1)
{
     long max_depth = 0; /* ROOT depth */
     long max_node = ROOT;
     long i, n, node;
     const long *depth;
     struct lpairs p;

     depth = get_nodes_depth();
     n = v1->nr;
     for (i = 0; i < n; i++) {
	  node = v1->data[i];
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

static bool *nodes_in_the_term(const struct long_array *terms, long max)
{
     bool *in_term;

     in_term = (bool *)xcalloc(max, sizeof(bool));
     memset(in_term, false, max*sizeof(bool));

     for (size_t i = 0; i < terms->nr; i++)
	  in_term[terms->data[i]] = true;

     return in_term;
}

static void compute_color_nodes(void *object, const struct long_array *v1,
				const struct long_array *v2,
				struct node_ptr_array *color_nodes, double threshold,
				const bool* in_terms1, const bool* in_terms2, bool matrix)
{
     long i, j, n, m, x, y, cont;
     double sim;
     struct node *new;
     double **M = NULL;
     struct graph *g = NULL;

     if (matrix) {
	  M = (double **)object;
     } else {
	  g = (struct graph *)object;
     }

     cont = 0;
     n = v1->nr;
     m = v2->nr;
     for (i = 0; i < n; i++) {
	  x = v1->data[i];
	  for (j = 0; j < m; j++) {
	       y = v2->data[j];
	       if (matrix)
		    sim = M[x][y];
	       else
		    sim = get_similarity(g,x,y);
	       if  (sim > threshold) {
		    if (!((in_terms2[x] && in_terms1[y]) && (x != y))) {
			 new = xcalloc(1, sizeof(struct node));
			 new->id = cont;
			 new->pos1 = i;
			 new->pos2 = j;
			 new->sim = sim;
			 new->cp = NULL;
			 ARRAY_PUSH(*color_nodes, new);
			 cont++;
		    }
	       }
	  }
     }
}

static void build_graph_to_coloring_matrix(struct graph *gc, 
					   struct node_ptr_array *vn, double **M, 
					   const struct long_array *v1, 
					   const struct long_array *v2,
					   double threshold1, double threshold2)
{
     long i, j, n, a, b, c, d, cont;
     struct node *x, *y;

     n = vn->nr;
     cont = 0;
     for (i = 0; i < n-1; i++) {
	  x = vn->data[i];
	  a = x->pos1;
	  b = x->pos2;
	  for (j = i+1; j < n; j++) {
	       y = vn->data[j];
	       c = y->pos1;
	       d = y->pos2;
	       assert(a <= c);
	       if ( (M[v1->data[a]][v1->data[c]] <= threshold1) ||
		    (M[v2->data[b]][v2->data[d]] <= threshold2)) {
		    add_arc_to_graph(gc, cont, i, j, COST);
		    add_arc_to_graph(gc, cont, j, i, COST);
		    cont++;
	       }
	  }
     }
}

static double **similarity_between_all(const struct graph *g, const struct long_array *v)
{
     long i, j, n;
     double **sim_matrix;

     n = v->nr;
     sim_matrix = double_matrix(0, n, 0, n);
     for (i = 0; i < n; i++) {
	  for (j = i; j < n; j++) {
	       sim_matrix[i][j] = get_similarity(g, v->data[i], v->data[j]);
	       /* I suppose that the similarity measure is symmetric */
	       sim_matrix[j][i] = sim_matrix[i][j];
	  }
     }
     return sim_matrix;
}

static void build_graph_to_coloring_graph(struct graph *gc, struct node_ptr_array *vn,
					  double **simM1, double **simM2,
					  double threshold1, double threshold2)
{
     long i, j, n, a, b, c, d, cont;
     struct node *x, *y;

     if (simM1 == NULL)
	  fatal("Matrix 1 is NULL");
     if (simM2 == NULL)
	  fatal("Matrix 2 is NULL");
     n = vn->nr;
     cont = 0;
     for (i = 0; i < n-1; i++) {
	  x = vn->data[i];
	  a = x->pos1;
	  b = x->pos2;
	  for (j = i+1; j < n; j++) {
	       y = vn->data[j];
	       c = y->pos1;
	       d = y->pos2;
	       assert(a <= c);
	       if ( (simM1[a][c] <= threshold1) || (simM2[b][d] <= threshold2)) {
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

static inline int compare_saturation(const struct graph *g, const saturation_t *a, const saturation_t *b)
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

static inline long get_color(const struct node_ptr_array *node_color, long pos)
{
     long color;
     struct node *node_ptr;

     if (pos < 0)
	  fatal("Invalid position");
     if (node_color == NULL)
	  fatal("Invalid pointer to color");
     node_ptr = node_color->data[pos];
     if (node_ptr->cp == NULL) {
	  color = NOCOLOR;
     } else {
	  color = node_ptr->cp->id;
     }
     return color;
}

static bool *get_ady_used_color(const struct graph *g, 
				const struct node_ptr_array *node_color, long node)
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
				     long node, const struct node_ptr_array *node_color)
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

static long greatest_saturation_node(const struct graph *g, pqueue_t *pq, 
				     const struct node_ptr_array *node_color)
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

static void get_free_colors(const struct graph *g, const struct node_ptr_array *solution,
			    long node, struct long_array *available_colors,
			    struct color_ptr_array *partitions)
{
     long i, cn, n;
     bool *color_used;
     long color_no_used;
     struct color *ctmp;

     n = partitions->nr;
     color_no_used = n+1;
     color_used = NULL;
     assert(g->n_nodes > node);
     assert(g->n_nodes >= color_no_used);
     color_used = get_ady_used_color(g, solution, node);
     cn = get_color(solution, node);
     if ((cn != NOCOLOR) && (!color_used[cn])) /* any adjacent vertex has the same color */
	  fatal("a adjacent node are using the same color");
     for (i = 0; i < n; i++) {
	  ctmp = partitions->data[i];
	  assert(ctmp->id == i);
	  if (!color_used[i]) {
	      available_colors->data[available_colors->nr++] = i;
	  }
     }
     if (available_colors->nr == 0)
	  available_colors->data[available_colors->nr++] = i;
     free(color_used);
}

static inline struct color *new_color(long id, double sim_e1, double sim_e2, 
				      double sim_bt, long node, long e1, long e2)
{
     struct color *new;
     struct long_array aux = {0, 0, NULL};

     new = xcalloc(1, sizeof(struct color));
     new->id = id;
     new->sim_entity_1 = sim_e1;
     new->sim_entity_2 = sim_e2;
     new->sim_between = sim_bt;
     new->id_nodes = aux;
     ARRAY_PUSH(new->id_nodes, node);
     new->entities1 = aux;
     ARRAY_PUSH(new->entities1, e1);
     new->entities2 = aux;
     ARRAY_PUSH(new->entities2, e2);

     return new;
}

static void free_color(struct color *c)
{
     free(c->id_nodes.data);
     free(c->entities1.data);
     free(c->entities2.data);
     free(c);
}

static double get_partition_similarity(struct color *c)
{
     long n, m;
     double bpe, annt1, annt2;

     n = c->entities1.nr;
     m = c->entities2.nr;
     bpe = (double)c->sim_between/c->id_nodes.nr;
     if (n > 1) {
	  annt1 = (double)c->sim_entity_1/((n*n-n)/2.0);
     } else {
	  annt1 = c->sim_entity_1;
	  assert(annt1 == 0.0);
     }
     if (m > 1) {
	  annt2 = (double)c->sim_entity_2/((m*m-m)/2.0);
     } else {
	  annt2 = c->sim_entity_2;
	  assert(annt2 == 0.0);
     }
     return (bpe + annt1 + annt2);
}


static double density_total_add(void *object, struct color_ptr_array *partitions,
				struct node *new_node, 
				const struct long_array *v1,
				const struct long_array *v2, 
				long color, bool matrix)
{
     long i, n, m, e1, e2, ep, n_colors;
     struct color *cptr;
     double annt1_nc, annt2_nc, bpe_nc, dt, sim, bpe, annt1, annt2;
     struct graph *g = NULL;
     double **M = NULL;

     if (matrix)
	  M = (double **)object;
     else
	  g = (struct graph *)object;

     dt = 0.0;
     cptr = partitions->data[color];
     assert(cptr->id == color);
     n = cptr->entities1.nr;
     e1 = v1->data[new_node->pos1];
     annt1_nc = 0.0;
     for (i = 0; i < n; i++) {
	  ep = cptr->entities1.data[i];
	  if (matrix)
	       annt1_nc += M[ep][e1];
	  else
	       annt1_nc += get_similarity(g, ep, e1);
     }

     n = cptr->entities2.nr;
     e2 = v2->data[new_node->pos2];
     annt2_nc = 0.0;
     for (i = 0; i < n; i++) {
	  ep = cptr->entities2.data[i];
	  if (matrix)
	       annt2_nc += M[ep][e2];
	  else
	       annt2_nc += get_similarity(g, ep, e2);
     }
     bpe_nc = cptr->sim_between;
     n_colors = partitions->nr;
     for (i = 0; i < n_colors; i++) {
	  cptr = partitions->data[color];
	  assert(cptr->id == color);
	  if (i == color) {
	       n = cptr->entities1.nr + 1;
	       m = cptr->entities2.nr + 1;
	       bpe = (double)(cptr->sim_between + bpe_nc)/(cptr->id_nodes.nr + 1.0);
	       if (n > 1) {
		    annt1 = (double)(cptr->sim_entity_1 + annt1_nc)/((n*n-n)/2.0);
	       } else {
		    annt1 = annt1_nc;
		    assert(cptr->sim_entity_1 == 0.0);
	       }
	       if (m > 1) {
		    annt2 = (double)(cptr->sim_entity_2 + annt2_nc)/((m*m-m)/2.0);
	       } else {
		    annt2 = annt2_nc;
		    assert(cptr->sim_entity_2 == 0.0);
	       }
	       sim = bpe + annt1 + annt2;
	  } else {
	       sim = get_partition_similarity(cptr);
	  }
	  dt += 1.0 - (sim/3.0);
     }
     assert(dt/n_colors <= 1.0);
     return dt;
}

static double density_total(struct color_ptr_array *partitions)
{
     long i, n_colors;
     double dt;
     struct color *c;

     dt = 0.0;
     n_colors = partitions->nr;
     for (i = 0; i < n_colors; i++) {
	  c = partitions->data[i];
	  dt += 1.0 - (get_partition_similarity(c)/3.0);
     }
     assert(dt/n_colors <= 1.0);
     return dt;
}

static long get_color_of_minimum_density(void *object, const struct node_ptr_array *vn,
					 struct color_ptr_array *partitions,
					 long new_node, struct long_array *free_colors,
					 const struct long_array *v1, 
					 const struct long_array *v2, bool matrix)
{
     long i, n, best_color, n_colors, new_color;
     double nc_best, nc, current_density;
     struct node *nptr;

     n = free_colors->nr;
     best_color = -1;
     nc_best = INFTY;
     n_colors = partitions->nr;
     current_density = density_total(partitions);
     nptr = vn->data[new_node];
     for (i = 0; i < n; i++) {
	  new_color = free_colors->data[i];
	  assert((new_color >= 0) && (new_color <= n_colors));
	  DEBUG("Color to evaluate %ld ", new_color);
	  if (n_colors == new_color) {
	       /* mew color */
	       nc = (double)(current_density + (1.0 - nptr->sim/3.0)) / (n_colors+1.0);
	  } else {
	       nc = (double)density_total_add(object, partitions, nptr, v1, v2, new_color, matrix)/n_colors;
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
     if ( (double)nc_best/n > 1.0) {
	  fatal("computed a density not valid");
     }
     return best_color;
}

static double compute_pairwise_sim(void *object, const struct long_array *v, bool matrix)
{
     long i, j, n, m;
     double sim;
     double **M = NULL;
     struct graph *g = NULL;

     if (matrix)
	  M = (double **)object;
     else
	  g = (struct graph *)object;
     sim = 0.0;
     m = v->nr;
     n = m - 1;
     for (i = 0; i < n; i++)  {
	  for (j = i+1; j < m; j++) {
	       if (matrix) {
		    sim += M[v->data[i]][v->data[j]];
	       } else {
		    sim += get_similarity(g, v->data[i], v->data[j]);
	       }
	  }
     }
     return sim;
}

/*
 * Return true if the entity was added
 */
static bool add_if_not_contains(struct long_array *v, long entity)
{
     for (unsigned i = 0; i < v->nr; i++) {
	  if (v->data[i] == entity)
	       return false;
     }
     ARRAY_PUSH(*v, entity);
     return true;
}

static void set_colors(void *object, struct color_ptr_array *partitions, 
		       long color, struct node *nptr,
		       const struct long_array *v1, const struct long_array *v2, bool matrix)
{
     long n_colors, e1, e2;
     struct color *cptr;
     bool added;

     if (color < 0)
	  fatal("Invalid color");
     if (partitions == NULL)
	  fatal("Invalid pointer to partition");
     n_colors = partitions->nr;
     if (n_colors == color) {
	  e1 = v1->data[nptr->pos1];
	  e2 = v2->data[nptr->pos2];
	  cptr = new_color(color, 0.0, 0.0, nptr->sim, nptr->id, e1, e2);
	  ARRAY_PUSH(*partitions, cptr);
	  nptr->cp = cptr;
     } else if (n_colors > color) {
	  cptr = partitions->data[color];
	  ARRAY_PUSH(cptr->id_nodes, nptr->id);

	  e1 = v1->data[nptr->pos1];
	  added = add_if_not_contains(&cptr->entities1, e1);
	  if (added)
	       cptr->sim_entity_1 = compute_pairwise_sim(object, &cptr->entities1, matrix);

	  e2 = v2->data[nptr->pos2];
	  added = add_if_not_contains(&cptr->entities2, e2);
	  if (added)
	       cptr->sim_entity_2 = compute_pairwise_sim(object, &cptr->entities2, matrix);

	  cptr->sim_between += nptr->sim;
	  nptr->cp = cptr;
     } else {
	  fatal("Bad color");
     }
}

static void coloring(void *object, const struct graph *g, 
		     const struct node_ptr_array *nodes,
		     struct color_ptr_array *partitions, 
		     const struct long_array *v1, 
		     const struct long_array *v2, bool matrix)
{
     struct color *cptr;
     struct node *nptr;
     long colored_nodes, new_node, e1, e2, color;
     pqueue_t pq_saturation;
     struct long_array free_colors = {0, 0, NULL};
     
     colored_nodes = 0;
     ALLOC_GROW(free_colors.data, (unsigned)g->n_nodes, free_colors.alloc);
     pq_init(&pq_saturation);
     init_saturation_pq(g, &pq_saturation);

     if (partitions->nr == 0)  {
	  new_node = greatest_saturation_node(g, &pq_saturation, nodes);
	  if (new_node == NONODE)
	       fatal("Error getting the greatest saturation node");
	  nptr = nodes->data[new_node];
	  e1 = v1->data[nptr->pos1];
	  e2 = v2->data[nptr->pos2];
	  assert(new_node == nptr->id);
	  cptr = new_color(0, 0.0, 0.0, nptr->sim, nptr->id, e1, e2);
	  ARRAY_PUSH(*partitions, cptr);
	  nptr->cp = cptr;
	  colored_nodes++;
	  if (pq_saturation.size != 0)
	       update_saturation_degree(g, &pq_saturation, new_node, nodes);
     }
     
     while (colored_nodes < g->n_nodes) {
	  new_node = greatest_saturation_node(g, &pq_saturation, nodes);
	  if (new_node == NONODE)
	       fatal("Error getting the greatest saturation node");
	  nptr = nodes->data[new_node];
	  assert(new_node == nptr->id);
	  free_colors.nr = 0;
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
     free(free_colors.data);
}

#ifdef PRGDEBUG
static void print_coloring(struct color_ptr_array *partitions)
{
     struct color *c;
     for (unsigned i = 0; i < partitions->nr; i++) {
	  c = partitions->data[i];
	  assert(c->id == i);
	  printf("%ld %.4f %.4f %.4f %zu %zu %zu\n", 
		 c->id, c->sim_entity_1, c->sim_entity_2, c->sim_between, 
		 c->id_nodes.nr, c->entities1.nr, c->entities2.nr);
     }
}
#endif

static double get_density_average(struct color_ptr_array *partitions)
{
     long i, n_colors;
     double dt;
     struct color *c;

     dt = 0.0;
     n_colors = partitions->nr;
     
     for (i = 0; i < n_colors; i++) {
	  c = partitions->data[i];
	  dt += (get_partition_similarity(c)/3.0);;
     }
     assert(dt/n_colors <= 1.0);
     return dt / partitions->nr;
 }

/*************************************
 *************************************
 **
 **  Get Predited Links
 **
 ************************************
 ************************************/

static void init_prediction_array(prediction_partition_t *a, unsigned n)
{
     a->pred = xmalloc(n * sizeof(*(a->pred)));
     for (unsigned i = 0; i < n; i++) {
	  a->pred[i].alloc = 0;
	  a->pred[i].nr = 0;
	  a->pred[i].data = NULL; 
     }
     a->nr = n;
}

static void free_prediction_array(prediction_partition_t *a)
{
     for (unsigned i = 0; i < a->nr; i++) {
	  free(a->pred[i].data);
	  a->pred[i].nr = 0;
	  a->pred[i].alloc = 0;
     }
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

static double get_cluster_probability(long nodes1, long nodes2, long n_edges)
{
     long n_nodes = nodes1 * nodes2;
     if (n_nodes == 0)
	  fatal("Zero Division");
     return (double)n_edges / n_nodes;
}

static void get_predicted_links(prediction_partition_t *cluster_pred,
				struct node_ptr_array *color_nodes,
				struct color_ptr_array *partitions,
				const struct long_array *v1, const struct long_array *v2)
{
     unsigned i, j, k, n_links, m, n, n_clusters;
     struct hash_map edges_obs[partitions->nr];
     struct hash_entry *hentry;
     char *buf;
     item_entry_t *item;
     int l;
     struct node *edge;
     struct color *cluster;
     long node1, node2;
     prediction_t ptemp;
     long n_edges[partitions->nr];
     double prob;

     n_clusters = partitions->nr;
     n_links = color_nodes->nr;
     for (i = 0; i < n_clusters; i++) {
	  hmap_create(&edges_obs[i], n_links);
	  n_edges[i] = 0;
     }
     for (i = 0; i < n_links; i++) {
	  edge = color_nodes->data[i];
	  l = asprintf(&buf, "%ld-%ld", v1->data[edge->pos1], v2->data[edge->pos2]);
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
	  cluster = partitions->data[i];
	  n = cluster->entities1.nr;
	  m = cluster->entities2.nr;
	  prob = get_cluster_probability(n, m, n_edges[i]);
	  assert((prob >= 0.0) && (prob <= 1.0));
	  for (j = 0; j < n; j++) {
	       node1 = cluster->entities1.data[j];
	       for (k = 0; k < m; k++) {
		    node2 = cluster->entities2.data[k];
		    l = asprintf(&buf, "%ld-%ld",node1, node2);
		    if (l == -1)
			 fatal("Error in make a new link");
		    assert((unsigned)l == strlen(buf));
		    hentry = hmap_find_member(&edges_obs[i], buf, l);
		    if (hentry == NULL) {
			 ptemp.entity1 = node1;
			 ptemp.entity2 = node2;
			 ptemp.prob = prob;
			 ARRAY_PUSH(cluster_pred->pred[i], ptemp);
		    }
		    free(buf);
	       }
	  }
     }
     for (i = 0; i < n_clusters; i++)
	  free_hash_map_item(&edges_obs[i]);
}

static void print_predicted_links(prediction_partition_t *cluster_pred,
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
	  n = cluster_pred->pred[i].nr;
	  if (n > 0) {
	       fprintf(f, "Cluster\t%u\n", i);
	       cont += n;
	  }
	  for (j = 0; j < n; j++) {
	       ptemp = cluster_pred->pred[i].data[j];
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

static unsigned print_singleton(FILE *f, bool *in_v, const struct long_array *v, char **desc)
{
     unsigned i, n;
     
     n = 0;
     for (i = 0; i < v->nr; i++) {
	  if (!in_v[i]) {
	       fprintf(f ,"%s\n", desc[v->data[i]]);
	       n++;
	  }
     }
     return n;
}

static char *print_output_files(struct node_ptr_array *color_nodes,
				struct color_ptr_array *partitions,
				const struct long_array *v1, const struct long_array *v2,
				double threshold_E1, double threshold_E2,
				const char name1[], const char name2[], char **desc)
{
     FILE *f;
     unsigned i, j, n, m;
     char *output1, *output2, *message;
     struct stat st;
     struct node *edge;
     struct color *cluster;
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

     in_v1 = (bool *)xcalloc(v1->nr, sizeof(bool));
     memset(in_v1, false, v1->nr*sizeof(bool));
     in_v2 = (bool *)xcalloc(v2->nr, sizeof(bool));
     memset(in_v2, false, v2->nr*sizeof(bool));

     if (!f)
	  fatal("No descriptor file specified, abort\n");
     for (i = 0; i < color_nodes->nr; i++) {
	  edge = color_nodes->data[i];
	  fprintf(f ,"%s\t%s\t%.4f\n", 
		  desc[v1->data[edge->pos1]],
		  desc[v2->data[edge->pos2]], edge->sim);
	  in_v1[edge->pos1] = true;
	  in_v2[edge->pos2] = true;
     }
     n = print_singleton(f, in_v1, v1, desc);
     m = print_singleton(f, in_v2, v2, desc);
     fclose(f);

     printf("Number of partitions: %u\n", partitions->nr);
     for (i = 0; i < partitions->nr; i++) {
	  cluster = partitions->data[i];
	  if (asprintf(&output2, "%s/%s-%s-%u-%.4f-%.4f.txt",
		       output1, name1, name2, i, threshold_E1, threshold_E2) == -1)
	       fatal("Error in cluster file");
          f = fopen(output2, "w");
	  if (!f)
	       fatal("No descriptor file specified, abort\n");
	  for (j = 0; j < cluster->id_nodes.nr; j++) {
	       id_node = cluster->id_nodes.data[j];
	       edge = color_nodes->data[id_node];
	       assert(cluster->id == edge->cp->id);
	       fprintf(f ,"%s\t%s\t%.4f\n", 
		       desc[v1->data[edge->pos1]],
		       desc[v2->data[edge->pos2]], edge->sim);
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
			    const struct long_array *v1, const struct long_array *v2,
			    double threshold_E1, double threshold_E2,
			    double threshold_bt, const char name1[], const char name2[],
			    char **desc, bool prediction, bool matrix, enum measure d)
{
     struct graph gc;
     bool *in_term1, *in_term2;
     unsigned i;
     double **sim_mat, **simM1, **simM2;
     struct graph *g;
     clock_t ti, tf;
     struct lpairs nd1, nd2;
     long max_depth;
     char *message;
     double sim;
     prediction_partition_t cluster_pred;
     struct color_ptr_array partitions = {0, 0, NULL};
     struct node_ptr_array color_nodes = {0, 0, NULL};
 
     ti = clock();
     sim_mat = NULL;
     simM1 = NULL;
     simM2 = NULL;
     sim = 0.0;

     if (matrix) {
	  sim_mat = (double **)object;
     } else {
	  g = (struct graph *)object;
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
     
     in_term1 = nodes_in_the_term(v1, n_nodes);
     in_term2 = nodes_in_the_term(v2, n_nodes);
     compute_color_nodes(object, v1, v2, &color_nodes, threshold_bt, in_term1, in_term2, matrix);
     free(in_term1);
     free(in_term2);
     init_graph(&gc, color_nodes.nr);
     if (matrix) {
	  build_graph_to_coloring_matrix(&gc, &color_nodes, sim_mat, v1, v2, threshold_E1, threshold_E2);
     } else {
	  build_graph_to_coloring_graph(&gc, &color_nodes, simM1, simM2, threshold_E1, threshold_E2);
	  free_double_matrix(simM1, 0, 0);
	  free_double_matrix(simM2, 0, 0);
     }
     tf = clock();
     printf("Bipartite Graph data - Nodes A: %u; Nodes B: %u; Edges: %u\n",
	    v1->nr, v2->nr, color_nodes.nr);
     printf("Time to build the graph to coloring: %.3f secs\n", (double)(tf-ti)/CLOCKS_PER_SEC);
#ifdef PRGDEBUG
     print_graph(&gc);
#endif
     printf("Graph to Coloring - Nodes: %ld; Edges: %ld\n", gc.n_nodes, gc.n_edges/2);
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
     message = print_output_files(&color_nodes, &partitions, v1, v2, 
				  threshold_E1, threshold_E2, name1, name2, desc);
     printf("%s", message);
     if (prediction) {
	  init_prediction_array(&cluster_pred, partitions.nr);
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
     }
     for (i = 0; i < partitions.nr; i++)
	  free_color(partitions.data[i]);
     free(partitions.data);
     for (i = 0; i < color_nodes.nr; i++)
	  free(color_nodes.data[i]);
     free(color_nodes.data);
     free_graph(&gc);

     return sim;
}
