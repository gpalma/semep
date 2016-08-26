/**
 * Copyright (C) 2013-2015 Universidad Simón Bolívar
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
#include "graph_adj.h"
#include "memory.h"
#include "hash_map.h"
#include "util.h"
#include "rbtree.h"
#include "semEP.h"

#define COST          0
#define NOCOLOR      -1
#define NONODE       -1
#define NOCLUSTER    -1
#define INFTY         INT_MAX

/*
 * Macros used for binaries heaps implementation
 */
#define PARENT(i)     ((int)(((i) - 1) / 2))
#define LEFT(i)       (((i) * 2) + 1)
#define RIGHT(i)      (((i) * 2) + 2)

struct color_ptr_array{
     unsigned nr;
     unsigned alloc;
     struct color **data;
};

typedef struct item_entry {
     double value;
     struct hash_entry entry;
} item_entry_t;

typedef struct set_type {
     int key_int;
     struct rb_node node;
} set_type_t;

typedef struct saturation {
     int node;
     int n_adj_colors;
     struct rb_root color_used;
} saturation_t;

typedef struct pqueue {
     int size;
     int *node_pos;
     saturation_t **heap;
} pqueue_t; 

typedef struct prediction {
     int cluster;
     int entity1;
     int entity2;
     double prob;
} prediction_t;

typedef struct prediction_array {
     unsigned nr;
     unsigned alloc;
     prediction_t *data;
} prediction_array_t;

/*
 * Variables to get the simimilarity
 */

static int sl, sr, el, er;
static double **ML;
static double **MR;

/***************************************************
 **  Function to get the similirity between terms
 ***************************************************/

static inline double similarity(int a, int b)
{
     double sim = 0.0;
     
     if (((a >= sl) && (a < el)) && ((b >= sl) && (b < el))) {
	  sim = ML[a][b];
     } else if ((a >= sr) && (a < er) && (b >= sr) && (b < er)) {
	  sim = MR[a][b];
     } else {
	  fatal("Error in computing the similarity between %d and %d\n", a, b);
     }
     return sim;
}

/*************************************
 *************************************
 **
 **  Build the graph to coloring
 **
 ************************************
 ************************************/

static void build_graph_to_coloring_matrix(struct graph_adj *gc, 
					   struct node_ptr_array *vn, 
					   const struct int_array *v1, 
					   const struct int_array *v2,
					   double threshold1, double threshold2)
{
     int i, j, n, a, b, c, d, cont;
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
	       if ( (similarity(v1->data[a], v1->data[c]) <= threshold1) ||
		    (similarity(v2->data[b], v2->data[d]) <= threshold2) ) {
		    add_arc(gc, cont, i, j);
		    add_arc(gc, cont, j, i);
		    cont++;
	       }
	  }
     }
}

/*************************************
 *************************************
 **
 **  Operations of a set implemented  
 ** as red black tree
 **
 ************************************
 ************************************/

static bool search_and_insert(struct rb_root *root, int key)
{
     struct rb_node **new;
     struct rb_node *parent;
     int result;
     set_type_t *tmp;
     set_type_t *data;
     
     parent = NULL;
     new = &(root->rb_node); 
     while (*new) {
	  tmp = container_of(*new, set_type_t, node);
	  result = tmp->key_int;
	  parent = *new; 
	  if (key < result)
	       new = &((*new)->rb_left);
	  else if (key > result)
	       new = &((*new)->rb_right);
	  else
	       return false;
     }
     data = (set_type_t *)xmalloc(sizeof(set_type_t)); 
     data->key_int = key;
     rb_link_node(&data->node, parent, new);
     rb_insert_color(&data->node, root);
     return true;
}

static void destroy_set(struct rb_root *mytree)
{
     struct rb_node *node, *ant;

     set_type_t *tmp;
     for (node = rb_first(mytree); node; ) {
	  tmp = rb_entry(node, set_type_t, node);
	  ant = node;
	  node = rb_next(node);
	  rb_erase(ant, mytree);
	  free(tmp);
     }
}

/*************************************
 *************************************
 **
 **  Coloration Solver
 **
 ************************************
 ************************************/

#ifdef PRGDEBUG
static void print_int_array(struct int_array *v)
{
     printf("\n");
     for (size_t i = 0; i < v->nr; i++) {
	  printf("%d ", v->data[i]);
     }
     printf("\n");
}
#endif

static inline int compare_saturation(const struct graph_adj *g, const saturation_t *a,
				     const saturation_t *b)
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
	  destroy_set(&node->color_used);
	  free(node);
     }
}

static inline void pq_init(pqueue_t *pq)
{
     pq->size = 0;
     pq->heap = NULL;
     pq->node_pos = NULL;
}

static inline void pq_delete(pqueue_t *pq)
{
     int i;
     
     for(i = 0; i < pq->size; i++)
	  free_saturation_node(pq->heap[i]);
     free(pq->node_pos);
     free(pq->heap);
}

static inline void pq_insert(const struct graph_adj *g, pqueue_t *pq, saturation_t *node)
{
     int i, p;
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

static int extract_max(const struct graph_adj *g, pqueue_t *pq, saturation_t **node)
{
     int i, j, l, r;
     saturation_t *aux;
     saturation_t **tmp;
          
     if(pq->size == 0)
	  return -1;

     *node = pq->heap[0];
     aux =  pq->heap[pq->size-1];
     SWAP(pq->node_pos[pq->heap[0]->node],
	  pq->node_pos[pq->heap[pq->size-1]->node]); /* SWAP the positions*/
     if((pq->size - 1) > 0){
	  tmp = (saturation_t **)xrealloc(pq->heap, (pq->size-1)*sizeof(saturation_t *));
	  pq->heap = tmp;
	  pq->size--;
     } else {
	  free(pq->heap);
	  pq->heap = NULL;
	  free(pq->node_pos);
	  pq->node_pos = NULL;
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
	       SWAP(pq->node_pos[pq->heap[j]->node],
		    pq->node_pos[pq->heap[i]->node]); /* SWAP the positions*/
	       SWAP(pq->heap[j], pq->heap[i]);
	       i = j;
	  }
     }
     return 0;
}

static int increase_key(const struct graph_adj *g, pqueue_t *pq, int node, int color)
{
     int i, p, pos;
	  
     if (pq->size == 0)
	  return -1;

     pos = pq->node_pos[node];
     if (pos >= pq->size)
	  pos = -1;
   
     if (pos == -1)
	  return -2;

     if (search_and_insert(&(pq->heap[pos]->color_used), color))
	  pq->heap[pos]->n_adj_colors++;
     else
	  return 0;
     
     i = pos;
     p = PARENT(i);
     while((i > 0) && (compare_saturation(g, pq->heap[p], pq->heap[i]) < 0)){
	  SWAP(pq->node_pos[pq->heap[p]->node],
	       pq->node_pos[pq->heap[i]->node]); /* SWAP the positions*/
	  SWAP(pq->heap[p], pq->heap[i]);
	  i = p;
	  p = PARENT(i);
     }
     return 0;
}

static void init_saturation_pq(const struct graph_adj *g, pqueue_t *pq)
{
     int i, n;
     saturation_t *ns;
     
     n = g->n_nodes;
     for (i = 0; i < n; i++) {
	  ns = (saturation_t *)xmalloc(sizeof(saturation_t));
	  ns->node = i;
	  ns->n_adj_colors = 0;
	  ns->color_used = RB_ROOT;
	  pq_insert(g, pq, ns);	  
     }
     assert(pq->size == n);
     pq->node_pos = (int *)xmalloc(n*sizeof(int));
     for (i = 0; i < n; i++) {
	  pq->node_pos[pq->heap[i]->node] = i;
     }
}

static inline int get_color(const struct node_ptr_array *node_color, int pos)
{
     int color;
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
 
static bool *get_ady_used_color(const struct graph_adj *g, 
				const struct node_ptr_array *node_color, int node)
{
     struct arc *current;
     bool *color_used;
     int color;
     size_t alloc;

     alloc = g->n_nodes*sizeof(bool)+1;
     color_used = xmalloc(alloc);
     memset(color_used, false, alloc);
     current = g->adj_list[node];
     while (current != NULL) {
	  color = get_color(node_color, current->to);
	  assert(color < g->n_nodes);
	  if (color != NOCOLOR) {
	       color_used[color] = true;
	  }
	  current = current->next;
     }
     return color_used;
}

static void update_saturation_degree(const struct graph_adj *g, pqueue_t *pq, 
				     int node, const struct node_ptr_array *node_color)
{
     int r;
     int color;
     struct arc *current;
     
     color = get_color(node_color, node);
     assert(color != NOCOLOR);
     current = g->adj_list[node];
     while(current != NULL) {
	  r = increase_key(g, pq, current->to, color);
	  if (r == -1)
	       fatal("Error in update saturation degree\n");
	  current = current->next;
     }
}

static int greatest_saturation_node(const struct graph_adj *g, pqueue_t *pq, 
				    const struct node_ptr_array *node_color)
{
     int r, color, node;
     saturation_t *ns;

     node = NONODE;
     ns = NULL;
     color = INFTY;
     r = extract_max(g, pq, &ns);
     if (r == -1)
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
     printf("Node %d degree %d din %d\n",node, ns->n_adj_colors, g->degree[node].din);
#endif
     free_saturation_node(ns); 
     return node;
}

static void get_free_colors(const struct graph_adj *g, const struct node_ptr_array *solution,
			    int node, struct int_array *available_colors,
			    struct color_ptr_array *partitions)
{
     int i, cn, n;
     bool *color_used;
     int color_no_used;
     struct color *ctmp;

     n = partitions->nr;
     color_no_used = n+1;
     color_used = NULL;
     assert(g->n_nodes > node);
     assert(g->n_nodes >= color_no_used);
     color_used = get_ady_used_color(g, solution, node);
     cn = get_color(solution, node);
     if (cn != NOCOLOR)
	  if (!color_used[cn]) /* any adjacent vertex has the same color */
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

static inline struct color *new_color(int id, double sim_e1, double sim_e2, 
				      double sim_bt, int node, int e1, int e2)
{
     struct color *new;
     struct int_array aux = {0, 0, NULL};

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
     int n, m;
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


static double density_total_add(struct color_ptr_array *partitions,
				struct node *new_node, 
				const struct int_array *v1,
				const struct int_array *v2, 
				int color)
{
     int i, n, m, e1, e2, ep, n_colors;
     struct color *cptr;
     double annt1_nc, annt2_nc, bpe_nc, dt, sim, bpe, annt1, annt2;
     
     dt = 0.0;
     cptr = partitions->data[color];
     assert(cptr->id == color);
     n = cptr->entities1.nr;
     e1 = v1->data[new_node->pos1];
     annt1_nc = 0.0;
     for (i = 0; i < n; i++) {
	  ep = cptr->entities1.data[i];
	  annt1_nc += similarity(ep, e1);
     }

     n = cptr->entities2.nr;
     e2 = v2->data[new_node->pos2];
     annt2_nc = 0.0;
     for (i = 0; i < n; i++) {
	  ep = cptr->entities2.data[i];
	  annt2_nc += similarity(ep, e2);
     }
     bpe_nc = new_node->sim;
     n_colors = partitions->nr;
     for (i = 0; i < n_colors; i++) {
	  cptr = partitions->data[i];
	  assert(cptr->id == i);
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
     int i, n_colors;
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

static int get_color_of_minimum_density(const struct node_ptr_array *vn,
					struct color_ptr_array *partitions,
					int new_node, struct int_array *free_colors,
					const struct int_array *v1, 
					const struct int_array *v2)
{
     int i, n, best_color, n_colors, new_color;
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
	  DEBUG("Color to evaluate %d ", new_color);
	  if (n_colors == new_color) {
	       /* mew color */
	       nc = (double)(current_density + (1.0 - nptr->sim/3.0)) / (n_colors+1.0);
	  } else {
	       nc = (double)density_total_add(partitions, nptr, v1, v2, new_color)/n_colors;
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

static double compute_pairwise_sim(const struct int_array *v)
{
     int i, j, n, m;
     double sim;

     sim = 0.0;
     m = v->nr;
     n = m - 1;
     for (i = 0; i < n; i++)  {
	  for (j = i+1; j < m; j++) {
	       sim += similarity(v->data[i], v->data[j]);
	  }
     }
     return sim;
}

/*
 * Return true if the entity was added
 */
static bool add_if_not_contains(struct int_array *v, int entity)
{
     for (unsigned i = 0; i < v->nr; i++) {
	  if (v->data[i] == entity)
	       return false;
     }
     ARRAY_PUSH(*v, entity);
     return true;
}

static void set_colors(struct color_ptr_array *partitions, 
		       int color, struct node *nptr,
		       const struct int_array *v1, const struct int_array *v2)
{
     int n_colors, e1, e2;
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
	       cptr->sim_entity_1 = compute_pairwise_sim(&cptr->entities1);

	  e2 = v2->data[nptr->pos2];
	  added = add_if_not_contains(&cptr->entities2, e2);
	  if (added)
	       cptr->sim_entity_2 = compute_pairwise_sim(&cptr->entities2);

	  cptr->sim_between += nptr->sim;
	  nptr->cp = cptr;
     } else {
	  fatal("Bad color");
     }
}

static void coloring(const struct graph_adj *g, 
		     const struct node_ptr_array *nodes,
		     struct color_ptr_array *partitions, 
		     const struct int_array *v1, 
		     const struct int_array *v2)
{
     struct color *cptr;
     struct node *nptr;
     int colored_nodes, new_node, e1, e2, color;
     pqueue_t pq_saturation;
     struct int_array free_colors = {0, 0, NULL};
	  
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
	  printf("Free colors:\n");
	  print_int_array(&free_colors);
#endif
	  color = get_color_of_minimum_density(nodes, partitions, new_node, &free_colors, v1, v2);
	  set_colors(partitions, color, nptr, v1, v2);
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
	  assert((unsigned)c->id == i);
	  printf("%d %.4f %.4f %.4f %u %u %u\n", 
		 c->id, c->sim_entity_1, c->sim_entity_2, c->sim_between, 
		 c->id_nodes.nr, c->entities1.nr, c->entities2.nr);
     }
}
#endif

static double get_density_average(struct color_ptr_array *partitions)
{
     int i, n_colors;
     double dt;
     struct color *c;

     dt = 0.0;
     n_colors = partitions->nr;
     
     for (i = 0; i < n_colors; i++) {
	  c = partitions->data[i];
	  dt += (get_partition_similarity(c)/3.0);;
     }
     return dt / partitions->nr;
}

/*************************************
 *************************************
 **
 **  Get Predited Links
 **
 ************************************
 ************************************/

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

static double get_cluster_probability(int nodes1, int nodes2, int n_edges)
{
     int n_nodes = nodes1 * nodes2;

     if (n_nodes == 0)
	  fatal("Zero Division");
     return (double)n_edges / n_nodes;
}

static void get_predicted_links(prediction_array_t *cluster_pred,
				struct node_ptr_array *color_nodes,
				struct color_ptr_array *partitions,
				const struct int_array *v1, const struct int_array *v2)
{
     unsigned int i, j, k, n_links, m, n, n_clusters;
     struct hash_map edges_obs;
     struct hash_entry *hentry;
     char *buf;
     item_entry_t *item;
     struct node *edge;
     struct color *cluster;
     int node1, node2, l;
     prediction_t ptemp;
     int n_edges[partitions->nr];
     double prob;

     n_clusters = partitions->nr;
     n_links = color_nodes->nr;
     hmap_create(&edges_obs, 2*n_links+1);
     memset(n_edges, 0, partitions->nr*sizeof(int));     
     for (i = 0; i < n_links; i++) {
	  edge = color_nodes->data[i];
	  l = asprintf(&buf, "%d-%d-%d", edge->cp->id, v1->data[edge->pos1], v2->data[edge->pos2]);
	  if (l == -1)
	       fatal("Error in edge key creation");
	  item = (item_entry_t *)xmalloc(sizeof(item_entry_t));
	  assert((unsigned)edge->cp->id <= n_clusters); 
	  n_edges[edge->cp->id]++;
	  if (hmap_add_if_not_member(&edges_obs, &item->entry, buf, l) != NULL)
	       fatal("Error, repeat edge in bipartite graph\n");
	  free(buf);
     }
     for (i = 0; i < n_clusters; i++) {
	  cluster = partitions->data[i];
	  assert((unsigned)cluster->id == i);
	  n = cluster->entities1.nr;
	  m = cluster->entities2.nr;
	  prob = get_cluster_probability(n, m, n_edges[i]);
	  assert((prob >= 0.0) && (prob <= 1.0));
	  for (j = 0; j < n; j++) {
	       node1 = cluster->entities1.data[j];
	       for (k = 0; k < m; k++) {
		    node2 = cluster->entities2.data[k];
		    l = asprintf(&buf, "%d-%d-%d", i, node1, node2);
		    if (l == -1)
			 fatal("Error in made a new link");
		    hentry = hmap_find_member(&edges_obs, buf, l);
		    if (hentry == NULL) {
			 ptemp.cluster = i;
			 ptemp.entity1 = node1;
			 ptemp.entity2 = node2;
			 ptemp.prob = prob;
			 ARRAY_PUSH(*cluster_pred, ptemp);
		    }
		    free(buf);
	       }
	  }
     }
     free_hash_map_item(&edges_obs);
}

static void print_predicted_links(prediction_array_t *cluster_pred,
				  double threshold_lf, double threshold_rg,
				  const char *name, char **desc)
{
     unsigned i, n;
     prediction_t ptemp;
     FILE *f;
     char *output;
     int current;
     
     if (asprintf(&output, "%s-%.4f-%.4f-Predictions.txt", name, threshold_lf, threshold_rg) == -1)
	  fatal("Error in prediction file");
     f = fopen(output, "w");
     free(output);
     if (!f)
	  fatal("No descriptor file specified, abort\n");
     n = cluster_pred->nr;
     current = NOCLUSTER;
     for (i = 0; i < n; i++) {
	  ptemp = cluster_pred->data[i];
	  if (current != ptemp.cluster) {
	       current = ptemp.cluster;
	       fprintf(f, "Cluster\t%u\n", current);
	  }
	  fprintf(f, "%s\t%s\t%.4f\n", desc[ptemp.entity1], desc[ptemp.entity2], ptemp.prob);
     }
     printf("Number of links predicted: %d\n", n);
     fclose(f);
}

/*************************************
 *************************************
 **
 **  Output files
 **
 ************************************
 ************************************/

static char *print_output_files(struct node_ptr_array *color_nodes,
				struct color_ptr_array *partitions,
				const struct int_array *v1, const struct int_array *v2,
				double threshold_lf, double threshold_rg,
				const char *name, char **desc)
{
     FILE *f;
     unsigned i, j, n, m;
     char *output1, *output2, *message;
     struct stat st;
     struct node *edge;
     struct color *cluster;
     int id_node;
    
     if (asprintf(&output1, "%s-%.4f-%.4f-Clusters", name, threshold_lf, threshold_rg) == -1)
	  fatal("Error in output directory");
     if (stat(output1, &st) == -1)
	  mkdir(output1, 0700);
 
     if (asprintf(&message, "Cluster directory: %s\n", output1) == -1)
	  fatal("Error in output message");

     printf("Number of partitions: %u\n", partitions->nr);
     n = partitions->nr;
     for (i = 0; i < n; i++) {
	  cluster = partitions->data[i];
	  if (asprintf(&output2, "%s/%s-%u-%.4f-%.4f.txt",
		       output1, name, i, threshold_lf, threshold_rg) == -1)
	       fatal("Error in cluster file");
          f = fopen(output2, "w");
	  if (!f)
	       fatal("No descriptor file specified, abort\n");
	  m = cluster->id_nodes.nr;
	  for (j = 0; j < m; j++) {
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
     free(output1);
     return message;
}

/*************************************
 *************************************
 **
 **  Solver
 **
 ************************************
 ************************************/

double semEP_solver(const struct matrix *lmatrix, const struct matrix *rmatrix,
		    const struct int_array *lterms, const struct int_array *rterms,
		    const struct string_array *desc, struct node_ptr_array *color_nodes,
		    double lthreshold, double rthreshold, const char *bpgraph_name)

{
     struct graph_adj gc;
     unsigned i;
     clock_t ti, tf;
     char *message;
     double density;
     prediction_array_t cluster_pred = {0, 0, NULL};
     struct color_ptr_array partitions = {0, 0, NULL};
     
     ti = clock();
     density = 0.0;
     ML =  lmatrix->data;
     sl = lmatrix->start;
     el = lmatrix->end;
     MR =  rmatrix->data;
     sr = rmatrix->start;
     er = rmatrix->end;
     init_graph_adj(&gc, color_nodes->nr);
     printf("Bipartite Graph data - Nodes A: %u; Nodes B: %u; Edges: %u\n",
	    	    lterms->nr, rterms->nr, color_nodes->nr);
     build_graph_to_coloring_matrix(&gc, color_nodes, lterms, rterms, lthreshold, rthreshold);
     tf = clock();
     printf("Time to build the graph to coloring: %.3f secs\n", (double)(tf-ti)/CLOCKS_PER_SEC);
#ifdef PRGDEBUG
     print_graph_adj(&gc);
#endif
     printf("Graph to Coloring - Nodes: %d; Edges: %d\n", gc.n_nodes, gc.n_arcs/2);
     ti = clock();
     if (gc.n_nodes != 0) {
	  coloring(&gc, color_nodes, &partitions, lterms, rterms);
	  density = get_density_average(&partitions);
#ifdef PRGDEBUG
	  print_coloring(&partitions);
#endif
     } else {
	  fatal("Graph to coloring has no nodes");
     }
     tf = clock();
     printf("Coloring solver time %.3f secs\n", (double)(tf-ti)/CLOCKS_PER_SEC);
     message = print_output_files(color_nodes, &partitions, lterms, rterms, 
				  lthreshold, rthreshold, bpgraph_name, desc->data);
     printf("%s", message);
     get_predicted_links(&cluster_pred, color_nodes, &partitions, lterms, rterms);
     print_predicted_links(&cluster_pred, lthreshold, rthreshold, bpgraph_name, desc->data);
     free_array(cluster_pred);
     free(message);
     for (i = 0; i < partitions.nr; i++)
	  free_color(partitions.data[i]);
     free(partitions.data);
     free_graph_adj(&gc);

     return density;
}
