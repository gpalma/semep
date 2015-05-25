/**
 * Copyright (C) 2012, Universidad Simón Bolívar
 *
 * @brief Implementation of Abstract Data Type Graph
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
#include "graph.h"
#include "util.h"
#include "memory.h"
#include "graph.h"

#define COST        1
#define ROOT        0
#define INFTY       INT_MAX
#define NS          -1
#define ZERO        0
#define BUFSZ       256

typedef enum {WHITE, GRAY, BLACK} color_e;

/******************************************************
*******************************************************
**
** Graph Basic Operations
**
*******************************************************
*******************************************************/

void set_edge(struct edge *e, long id, long from, long to, long c)
{
     e->id = id;
     e->from = from;
     e->to = to;
     e->cost = c;
}

void init_graph(struct graph *g, long nodes)
{
     long i, n;

     n = nodes + 1;
     g->n_buckets = n;
     g->n_nodes = nodes;
     g->n_edges = 0;
     g->adj_list = (struct edge_list *)xcalloc(n, sizeof(struct edge_list));
     g->degree = (struct valency *)xcalloc(n, sizeof(struct valency));
     for (i = 0; i < n; i++) {
	  set_edge(&(g->adj_list[i].item), -1, 0, 0, 0); /* dumb edge */
	  INIT_LIST_HEAD(&(g->adj_list[i].list));
	  g->degree[i].din = 0;
	  g->degree[i].dout = 0;
     }
}

struct edge_list *get_adjacent_list(struct graph *g, long from)
{

     if ((g->n_nodes <= from) || (from < 0))
	  fatal("Error, node %ld not found in the graph\n", g->n_nodes);

     return &(g->adj_list[from]);
}

struct edge_list *new_edge_list(long id, long from, long to, long c)
{
     struct edge_list *new;

     new = (struct edge_list *)xmalloc(sizeof(struct edge_list));
     set_edge(&(new->item), id, from, to, c);
     INIT_LIST_HEAD(&new->list);

     return new;
}

extern void add_arc_to_graph(struct graph *g, long id, long from, long to, long c)
{
     struct edge_list *adj_list, *new;

     adj_list = get_adjacent_list(g, from);
     new = new_edge_list(id, from, to, c);
     list_add_tail(&new->list, &adj_list->list);
     g->degree[from].dout++;
     g->degree[to].din++;
     g->n_edges++;
}

void add_reprensentative_ancestor(struct graph *g)
{
     unsigned long i, n;
     struct long_array zdegree = {0, 0, NULL};

     n = g->n_nodes;
     if (g->degree[0].din != 0)
	  fatal("The root no has in-degree equal to zero\n");
  
     for (i = 1; i < n; i++) 
	  if (g->degree[i].din == 0) 
	       ARRAY_PUSH(zdegree, i);
     
     if (zdegree.nr > 0) {
#ifdef PRGDEBUG
	  printf("\n ** Graph with %ld nodos with in-degree zero ** \n", zdegree.nr);
#endif
	  for (i = 0; i < zdegree.nr; i++) 
	       add_arc_to_graph(g, g->n_edges, ROOT, zdegree.data[i], COST);
     }
#ifdef PRGDEBUG
     else {
	  printf("\n ** DAG with one root ** \n");
     }
#endif
     free(zdegree.data);
}

void print_graph(const struct graph *g)
{
     long i, cont, n;
     struct edge_list *tmp;

     n = g->n_nodes;
     printf("\n***** Graph *****\n");
     printf("Num. of Nodes %ld --- Num. of Edges %ld\n",g->n_nodes, g->n_edges);
     for (i = 0; i < n; i++) {
	  printf("Node %ld - d in-out (%ld, %ld): ", i,  g->degree[i].din, g->degree[i].dout);
	  cont = 0;
	  adj_for_each(tmp, g->adj_list[i]) {
	       if (cont % 4 == 0 && cont != 0)
		    printf("\n");
	       printf("(id %ld f %ld t %ld) ", tmp->item.id, tmp->item.from, tmp->item.to);
	       cont++;
	  }
	  printf("\n");
     }
}

void set_edges_cost(const struct graph *g)
{
     long i, n, c1, c2;
     struct edge_list *tmp;

     n = g->n_buckets;
     for (i = 0; i < n; i++) {
	  adj_for_each(tmp, g->adj_list[i]) {
	       if (g->degree[tmp->item.from].dout == 0)
		    c1 = 1;
	       else
		    c1 = g->degree[tmp->item.from].dout;

	       if (g->degree[tmp->item.to].dout == 0)
		    c2 = 1;
	       else
		    c2 = g->degree[tmp->item.to].dout;
	       tmp->item.cost = (long) (log10(c1) + log10(c2)) + 1;
	       assert(tmp->item.cost > 0);
	  }
     }
}

void free_edge_item(struct edge_list *el)
{
     list_del(&el->list);
}

void free_graph(struct graph *g)
{
     struct edge_list *tmp;
     long i, n;

     graph_for_each_safe(tmp, g) {
	  free_edge_item(tmp);
	  free(tmp);
     }
     n = g->n_buckets;
     for (i = 0; i < n; i++) {
	  list_del(&(g->adj_list[i].list));
     }
     free(g->degree);
     free(g->adj_list);
     g->n_nodes = 0;
     g->n_edges = 0;
     g->n_buckets = 0;
}

bool find_edge(struct graph *g, long v1, long v2)
{
     bool is_edge;
     struct edge_list *tmp;

     is_edge = false;
     adj_for_each(tmp, g->adj_list[v1]) {
	  if (tmp->item.to == v2) {
	       is_edge = true;
	       break;
	  }
     }
     return is_edge;
}

/******************************************************
*******************************************************
**
** Graph Search Algorithms
**
*******************************************************
*******************************************************/

static void dfs_visit(const struct graph *g, long u, long *d,
		      long *f, long *pred, color_e *color, long *ctr)
{
     long v;
     struct edge_list *etmp;

     color[u] = GRAY;
     d[u] = ++(*ctr);
     adj_for_each(etmp, g->adj_list[u]) {
	  v = etmp->item.to;
	  if (color[v] == WHITE) {
	       pred[v] = u;
	       dfs_visit(g, v, d, f, pred, color, ctr);
	  }
     }
     color[u] = BLACK;
     f[u] = ++(*ctr);
}

void dfs_search(const struct graph *g, long s, long *d, long *f, long *pred)
{
     long i, n, ctr;
     color_e *color;

     n = g->n_nodes;
     ctr = 0;
     color = (color_e *)xmalloc(n*sizeof(color_e));
     for (i = 0; i < n; i++) {
	  color[i] = WHITE;
	  pred[i] = NS;
	  d[i] = NS;
	  f[i] = NS;
     }

     dfs_visit(g, s, d, f, pred, color, &ctr);
     /*for (i = 0; i < n; i++) {
       if (color[i] == WHITE) {
       dfs_visit(g, i, d, f, pred, color, &ctr);
       }
       }*/
     free(color);
}

static void dfs_min(const struct graph *g, long x, long t, long *cmin, bool *visit)
{
     long y;
     struct edge_list *etmp;

     visit[x] = true;
     if (x == t) {
	  cmin[x] = 0;
     } else {
	  adj_for_each(etmp, g->adj_list[x]) {
	       y = etmp->item.to;
	       if (!visit[y]) {
		    dfs_min(g, y, t, cmin, visit);
	       }
	  }
	  adj_for_each(etmp, g->adj_list[x]) {
	       y = etmp->item.to;
	       if ((cmin[y] != INFTY) && (cmin[x] > cmin[y] + etmp->item.cost)) {
		    cmin[x] = cmin[y] + etmp->item.cost;
	       }
	  }
     }
}

long min_distance(const struct graph *g, long s, long t)
{
     long i, n, min;
     long *cmin;
     bool *visit;

     n = g->n_nodes;
     assert(n > 0);
     visit = (bool *)xmalloc(n * sizeof(bool));
     cmin = (long *)xmalloc(n * sizeof(long));
     memset(visit, false, n*sizeof(bool));
     for (i = 0; i < n; i++) {
	  cmin[i] = INFTY;
     }
     dfs_min(g, s, t, cmin, visit);
     min = cmin[s];
     free(visit);
     free(cmin);

     return min;
}

void all_pairs_shortest(const struct graph *g, long **dist)
{
     long i, j, k;
     long long new_dist;
     struct edge_list *etmp;

     for (i = 0; i < g->n_nodes; i++) {
	  for (j = 0; j < g->n_nodes; j++) {
	       dist[i][j] = INFTY;
	  }
	  dist[i][i] = 0;
	  adj_for_each(etmp, g->adj_list[i]){
	       k = etmp->item.to;
	       dist[i][k] = etmp->item.cost;
	  }
     }
     for (k = 0; k < g->n_nodes; k++) {
	  for (i = 0; i < g->n_nodes; i++) {
	       if (dist[i][k] == INFTY)
		    continue;
	       for (j = 0; j < g->n_nodes; j++) {
		    new_dist = dist[i][k]; /* avoid overflow of infinity distance */
		    new_dist += dist[k][j];
		    if (new_dist < dist[i][j]) {
			 dist[i][j] = new_dist;
		    }
	       }
	  }
     }
}

void graph_inverse(const struct graph *orig, struct graph *inv)
{
     struct edge_list *tmp;

     init_graph(inv, orig->n_nodes);
     graph_for_each(tmp, orig) {
	  add_arc_to_graph(inv, tmp->item.id, tmp->item.to, tmp->item.from, tmp->item.cost);
     }
     inv->n_nodes = orig->n_nodes;
}

void graph_undirect(const struct graph *orig, struct graph *ud)
{
     struct edge_list *tmp;

     init_graph(ud, orig->n_nodes);
     graph_for_each(tmp, orig) {
	  add_arc_to_graph(ud, tmp->item.id, tmp->item.to, tmp->item.from, tmp->item.cost);
	  add_arc_to_graph(ud, tmp->item.id, tmp->item.from, tmp->item.to, tmp->item.cost);
     }
     ud->n_nodes = orig->n_nodes;
     ud->n_edges = ud->n_edges/2;
}

static void dfs_max(const struct graph *g, long x, long t, long *cmax, bool *visit)
{
     long y;
     struct edge_list *etmp;

     visit[x] = true;
     if (x == t) {
	  cmax[x] = 0;
     } else {
	  adj_for_each(etmp, g->adj_list[x]) {
	       y = etmp->item.to;
	       if (!visit[y]) {
		    dfs_max(g, y, t, cmax, visit);
	       }
	  }
	  adj_for_each(etmp, g->adj_list[x]) {
	       y = etmp->item.to;
	       //		printf("Entrando x %ld - y %ld\n", x, y);
	       if ((cmax[y] != NS) && (cmax[x] <= cmax[y] + etmp->item.cost)) {
		    //printf("x %ld - y %ld - cy %ld nc %ld\n", x, y, cmax[y], (cmax[x] + COST));
		    cmax[x] = cmax[y] + etmp->item.cost;
	       }

	  }

     }
}

long max_distance(const struct graph *g, long s, long t)
{
     long i, n, max;
     long *cmax;
     bool *visit;

     n = g->n_nodes;
     visit = (bool *)xmalloc(n * sizeof(bool));
     cmax = (long *)xmalloc(n * sizeof(long));
     memset(visit, false, n*sizeof(bool));
     for (i = 0; i < n; i++) {
	  cmax[i] = NS;
     }
     dfs_max(g, s, t, cmax, visit);
     max = cmax[s];
     free(visit);
     free(cmax);

     return max;
}

/******************************************************
*******************************************************
**
** Dijkstra algorithm
**
*******************************************************
*******************************************************/

/*
 * Macros used for binaries heaps implementations
 */
#define PARENT(i)     ((long)(((i) - 1) / 2))
#define LEFT(i)       (((i) * 2) + 1)
#define RIGHT(i)      (((i) * 2) + 2)

/*********************************
 **  Constants
 *********************************/

#define ARC_ID    -1

/*********************************
 ** Structures
 *********************************/

enum node_set {
     OPEN,
     CLOSED,
     NONE
};

struct node_record {
     long node;
     long from_node;
     struct edge *edge;
     long cost_so_far;
     enum node_set category;
};

struct pqueue {
     long size;
     struct node_record **heap;
};

/*********************************
 ** Priority Queue (Binary Heap)
 *********************************/


static inline void pq_init(struct pqueue *pq)
{
     pq->size = 0;
     pq->heap = NULL;
}

/*
  static inline void pq_delete( struct pqueue *pq )
  {
  int i;
  for( i = 0; i < pq->size; i++ )
  free( pq->heap[i] );
  free( pq->heap );
  memset( pq, 0, sizeof( struct pqueue ) );
  }
*/

static inline void pq_insert(struct pqueue *pq, struct node_record *node)
{
     long i, p;
     struct node_record **tmp;

     tmp = xrealloc(pq->heap, (pq->size+1)*sizeof(struct node_record *));
     pq->heap = tmp;
     pq->heap[pq->size] = node;
     i = pq->size;
     p = PARENT(i);
     while((i > 0) && ((pq->heap[i]->cost_so_far
			- pq->heap[p]->cost_so_far) < 0)){
	  SWAP(pq->heap[p], pq->heap[i]);
	  i = p;
	  p = PARENT(i);
     }
     pq->size++;
}

static long extract_min(struct pqueue *pq, struct node_record **node)
{
     long i, j, l, r;
     struct node_record *aux;
     struct node_record **tmp;

     if( pq->size == 0 )
	  return -1;
     *node = pq->heap[0];
     aux =  pq->heap[pq->size-1];
     if( (pq->size - 1) > 0 ){
	  tmp = xrealloc(pq->heap, (pq->size-1)*sizeof(struct node_record *));
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
	  if((l < pq->size) && ((pq->heap[i]->cost_so_far
				 - pq->heap[l]->cost_so_far) > 0))
	       j = l;
	  else
	       j = i;

	  if((r < pq->size) && ((pq->heap[j]->cost_so_far
				 - pq->heap[r]->cost_so_far) > 0))
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

static long decrease_key(struct pqueue *pq, long node, long new_cost)
{
     long i, p, pos;

     if( pq->size == 0 )
	  return -1;
     pos = -1;
     for( i = 0; i < pq->size; i++ )
	  if( pq->heap[i]->node == node )
	       pos = i;
     if( pos == -1 )
	  return -1;

     if( new_cost > pq->heap[pos]->cost_so_far )
	  return -1;
     pq->heap[pos]->cost_so_far = new_cost;
     i = pos;
     p = PARENT(i);
     while((i > 0) && ((pq->heap[i]->cost_so_far
			- pq->heap[p]->cost_so_far) < 0)){
	  SWAP(pq->heap[p], pq->heap[i]);
	  i = p;
	  p = PARENT(i);
     }
     return 0;
}

struct edge *new_arc(int id, long from, long to, long c)
{
     struct edge *e;

     e = xcalloc(1, sizeof(struct edge));
     e->id = id;
     e->from = from;
     e->to = to;
     e->cost = c;

     return e;
}
/*
  static void copy_arc(const struct edge *from, struct edge *to)
  {
  to->id = from->id;
  to->from = from->from;
  to->to = from->to;
  to->cost = from->cost;
  }
*/
/*********************************
 ** Dijkstra Algorithm
 *********************************/

/**
 * Dijkstra edge path with duplicate edges
 */
long min_path(const struct graph *g, long start, long goal)
{
     long min_dist;
     struct node_record *start_r, *current, *end_node_r;
     struct node_record **statusT;
     struct pqueue open;
     enum node_set contains;
     long result, end_node, end_node_cost;
     struct edge_list *connections;
     struct edge *etmp;
     long i;

     etmp = new_arc(ARC_ID, start, start, 0);

     start_r = xmalloc(sizeof(struct node_record));
     start_r->node  = start;
     start_r->from_node = 0;
     start_r->edge = etmp;
     start_r->cost_so_far = 0;
     start_r->category = OPEN;

     statusT = xcalloc(g->n_nodes, sizeof(struct node_record *));
     statusT[start_r->node] =  start_r;
     pq_init(&open);
     pq_insert(&open, start_r);
     current = NULL;
     end_node_r = NULL;

     while ( open.size > 0 ) {
	  result = extract_min(&open, &current);
	  assert(result != -1);
	  if (current->node == goal)
	       break;
	  adj_for_each(connections, g->adj_list[current->node]) {
	       end_node = connections->item.to;
	       end_node_cost = current->cost_so_far + connections->item.cost;
	       if (statusT[end_node] == NULL) {
		    contains = NONE;
	       } else {
		    end_node_r = statusT[end_node];
		    contains = end_node_r->category;
	       }

	       if (contains == CLOSED) {
		    continue;
	       } else if (contains == OPEN) {
		    if (end_node_r->cost_so_far <= end_node_cost)
			 continue;
		    end_node_r->node = end_node;
		    end_node_r->from_node = current->node;
		    end_node_r->edge = &connections->item;
		    end_node_r->cost_so_far = end_node_cost;
		    result = decrease_key(&open, end_node_r->node, end_node_r->cost_so_far);
		    assert(result != -1);
	       } else {
		    end_node_r = xmalloc(sizeof(struct node_record));
		    end_node_r->node = end_node;
		    end_node_r->from_node = current->node;
		    end_node_r->edge = &connections->item;
		    end_node_r->cost_so_far = end_node_cost;
		    end_node_r->category = OPEN;
		    pq_insert(&open, end_node_r);
		    statusT[end_node_r->node] = end_node_r;
	       }
	  }
	  current->category = CLOSED;
     }
     if (current->node != goal) {
	  min_dist = error("Error, goal no reached. Last Node: %u\n", current->node);
     } else {
	  min_dist = current->cost_so_far;
     }
     free(etmp);

     for (i = 0; i < g->n_nodes; i++)
	  if (statusT[i] != NULL)
	       free(statusT[i]);
     free(statusT);
     free(open.heap);

     return  min_dist;
}

/******************************************************
*******************************************************
**
** Topological sorting
**
*******************************************************
*******************************************************/

static void dfs_topological_sort(const struct graph *g, long u,
				 color_e *color, struct long_list *tpl_sort)
{
     long v;
     struct edge_list *etmp;
     struct long_list *ntmp;

     color[u] = GRAY;
     adj_for_each(etmp, g->adj_list[u]) {
	  v = etmp->item.to;
	  if (color[v] == WHITE) {
	       dfs_topological_sort(g, v, color, tpl_sort);
	  }
     }
     color[u] = BLACK;
     ntmp = (struct long_list *)xmalloc(sizeof(struct long_list));
     ntmp->item = u;
     list_add(&(ntmp->list), &(tpl_sort->list));
}

struct long_list *topological_sort(const struct graph *g, long s)
{
     long i, n;
     color_e *color;
     struct long_list *tpl_sort;

     n = g->n_nodes;
     color = (color_e *)xmalloc(n*sizeof(color_e));
     tpl_sort = (struct long_list *)xmalloc(sizeof(struct long_list));
     INIT_LIST_HEAD(&(tpl_sort->list));
     for (i = 0; i < n; i++) {
	  color[i] = WHITE;
     }
     dfs_topological_sort(g, s, color, tpl_sort);
     /*for (i = 0; i < n; i++) {
       if (color[i] == WHITE) {
       dfs_visit(g, i, d, f, pred, color, &ctr);
       }
       }*/

     free(color);

     return tpl_sort;
}

long *calculate_depth(const struct graph *g)
{
     long i, nodes, m, u, v;
     long *depth;
     struct long_list *tpl_sort;
     struct long_list *tmp;
     struct list_head *pos, *q;
     struct edge_list *etmp;

     m = 0;
     nodes = g->n_nodes;
     depth = (long *)xmalloc(nodes*sizeof(long));
     tpl_sort = topological_sort(g, ROOT);

     for (i = 0; i < nodes; i++) {
	  depth[i] = ZERO;
     }
     depth[ROOT] = ZERO;
     list_for_each_entry(tmp, &(tpl_sort->list), list) {
	  u = tmp->item;
	  adj_for_each(etmp, g->adj_list[u]) {
	       v = etmp->item.to;
	       if (depth[v] < (depth[u] + etmp->item.cost))
		    depth[v] = depth[u] + etmp->item.cost;
	  }
     }
     list_for_each_safe(pos, q, &(tpl_sort->list)){
	  tmp = list_entry(pos, struct long_list, list);
	  list_del(pos);
	  free(tmp);
	  m++;
     }
     assert(m = nodes);
     free(tpl_sort);
     return depth;
}

long *calculate_depth_bfs(const struct graph *g)
{
     long i, n, u, v;
     color_e *color;
     long *depth;
     struct long_list queue;
     struct long_list  *ntmp;
     struct edge_list *etmp;
     struct list_head *pos;

     n = g->n_nodes;
     color = (color_e *)xmalloc(n*sizeof(color_e));
     depth = (long *)xmalloc(n*sizeof(long));
     for (i = 0; i < n; i++) {
	  color[i] = WHITE;
	  depth[i] = NS;
     }
     INIT_LIST_HEAD(&(queue.list));
     ntmp = (struct long_list *)xmalloc(sizeof(struct long_list));
     ntmp->item = ROOT;
     list_add_tail(&(ntmp->list), &(queue.list));
     depth[ROOT] = 0;
     while ( !list_empty(&(queue.list)) ) {
	  pos = (&queue.list)->next;
	  ntmp = list_entry(pos, struct long_list, list);
	  u = ntmp->item;
	  list_del(pos);
	  free(ntmp);
	  adj_for_each(etmp, g->adj_list[u]) {
	       v = etmp->item.to;
	       if (color[v] == WHITE) {
		    depth[v] = depth[u] + 1;
		    color[v] = GRAY;
		    ntmp = (struct long_list *)xmalloc(sizeof(struct long_list));
		    ntmp->item = v;
		    list_add_tail(&(ntmp->list), &(queue.list));
	       } else {
		    if (depth[u] >= depth[v] ) {
			 depth[v] = depth[u] + 1;
			 ntmp = (struct long_list *)xmalloc(sizeof(struct long_list));
			 ntmp->item = v;
			 list_add_tail(&(ntmp->list), &(queue.list));
		    }
	       }
	  }
	  color[u] = BLACK;
     }
     return depth;
}

static bool visit(const struct graph *g, color_e *color, long v)
{
     long u;
     struct edge_list *etmp;

     color[v] = GRAY;
     adj_for_each(etmp, g->adj_list[v]) {
	  u = etmp->item.to;
	  if (color[u] == GRAY) {
	       return true;
	  } else {
	       if (color[u] == WHITE) {
		    if (visit(g, color, u)) {
			 return true;
		    }
	       }
	  }
     }
     color[v] = BLACK;
     return false;
}

bool detect_cycle(const struct graph *g)
{
     long i, n;
     color_e *color;

     n = g->n_nodes;
     color = (color_e *)xmalloc(n*sizeof(color_e));

     for (i = 0; i < n; i++) {
	  color[i] = WHITE;
     }
     for (i = ROOT; i < n; i++) {
	  if (color[i] == WHITE) {
	       if (visit(g, color, i)) {
		    free(color);
		    return true;
	       }
	  }
     }
     free(color);
     return false;
}

static void dfs_spanning_tree(const struct graph *g, color_e *color, long v, bool *st)
{
     long u;
     struct edge_list *etmp;

     color[v] = BLACK;
     adj_for_each(etmp, g->adj_list[v]) {
	  u = etmp->item.to;
	  if (color[u] == WHITE) {
	       assert(!st[etmp->item.id]);
	       st[etmp->item.id] = true;
	       dfs_spanning_tree(g, color, u, st);
	  }
     }
}

/**
 * Return a array with the id of the edges
 * that are part of the spanning_tree of the graph
 */
bool *get_spanning_tree(const struct graph *g)
{
     long i, n;
     color_e *color;
     bool *st;
     size_t alloc;

     alloc = g->n_edges*sizeof(bool);
     st = (bool *)xmalloc(alloc);
     memset(st, false, alloc);
     n = g->n_nodes;
     color = (color_e *)xmalloc(n*sizeof(color_e));

     for (i = 0; i < n; i++) {
	  color[i] = WHITE;
     }
     for (i = ROOT; i < n; i++) {
	  if (color[i] == WHITE) {
	       dfs_spanning_tree(g, color, i, st);
	  }
     }
     free(color);
     return st;
}

static void dfs_euler_tour(const struct graph *g, color_e *color,
			   long v, struct long_array *et)
{
     long u;
     struct edge_list *etmp;

     color[v] = BLACK;
     adj_for_each(etmp, g->adj_list[v]) {
	  u = etmp->item.to;
	  if (color[u] == WHITE) {
	       printf("depth %ld\n", u);
	       ARRAY_PUSH(*et, u);
	       dfs_euler_tour(g, color, u, et);
	       printf("bt %ld\n", v);
	       ARRAY_PUSH(*et, v);
	  }
     }
}

/**
 * Return a vector with the edges of the euler tour
 */

struct long_array *get_euler_tour(const struct graph *g)
{
     long i, n;
     color_e *color;
     struct long_array *etour;

     ALLOC_STRUCT(etour);
     n = g->n_nodes;
     color = (color_e *)xmalloc(n*sizeof(color_e));

     for (i = 0; i < n; i++) {
	  color[i] = WHITE;
     }
     for (i = ROOT; i < n; i++) {
	  if (color[i] == WHITE) {
	       printf("init add %ld\n", i);
	       ARRAY_PUSH(*etour, i);
	       dfs_euler_tour(g, color, i, etour);
	  }
     }
     free(color);
     return etour;
}

