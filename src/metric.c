/**
 * Copyright (C) 2012, 2013, 2014 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include "dlist.h"
#include "types.h"
#include "util.h"
#include "memory.h"
#include "graph.h"
#include "CA.h"
#include "metric.h"
#include "hash_map.h"

#define ROOT  0

typedef struct item_entry {
     long value;
     struct hash_entry entry;
} item_entry_t;

static long n;
static bool *visited;
static struct long_array **ancestors;
static long *depth;
static struct graph gi;
static bool init_metric = false;
static struct hash_map dis_map;
static long max_depth;

void init_metric_data(const struct graph *g)
{
     depth = calculate_depth(g);
     DEBUG("\n** Depth node calculation done ** \n");
     n = g->n_nodes;
     ancestors = xcalloc(n, sizeof(struct long_array *));
     graph_inverse(g, &gi);
     DEBUG("\n** Graph inverse done ** \n");
     visited = xcalloc(n, sizeof(bool));
     hmap_create(&dis_map, n);
     max_depth = INT_MAX;
     init_metric = true;
}

static long min_dag_distance(const struct graph *g, long x, long y)
{
     int l;
     char *key;
     item_entry_t *item;
     struct hash_entry *hentry;
     long d;
     
     if (x > y) {
	  l = asprintf(&key, "%ld-%ld", x, y);
     } else {
	  l = asprintf(&key, "%ld-%ld", y, x);
     }
     if (l == -1)
	  fatal("Error in output directory");
     
     hentry = hmap_find_member(&dis_map, key, l);
     if (hentry == NULL) {
	  d = min_distance(g, x, y);
	  item = (item_entry_t *)xmalloc(sizeof(item_entry_t));
	  item->value = d;
	  if (hmap_add(&dis_map, &item->entry, key, l) != 0)
	       fatal("Error in distance hash table");
     } else {
	  item = hash_entry(hentry, item_entry_t, entry);
	  d = item->value;
     }
     free(key);
     return d;
}

static inline double dtax(long dax, long day, long drx, long dry)
{
     double r = (double)(dax+day)/(drx+dry);

     return MIN(r, 1.0);
}

static struct long_array *get_list_ancestors(long node)
{
     struct long_array *la;

     if (!visited[node]) {
	  la = get_ancestors(&gi, node);
	  ancestors[node] = la;
	  visited[node] = true;
     } else {
	  la = ancestors[node];
     }

     return la;
}

double dist_tax(const struct graph *g, long x, long y)
{
     long lca, dax, day, drx, dry;
   
     if (!init_metric)
	  fatal("Error, uninitialized data for metric calcule");

     lca = LCA_CA(get_list_ancestors(x), get_list_ancestors(y), depth);
     dax = min_dag_distance(g, lca, x);
     day = min_dag_distance(g, lca, y);
     drx = min_dag_distance(g, ROOT, x);
     dry = min_dag_distance(g, ROOT, y);

     return dtax(dax, day, drx, dry);
}

double dist_tax_lca(const struct graph *g, long x, long y, long *lcap)
{
     long lca, dax, day, drx, dry;
      
     if (!init_metric)
	  fatal("Error, uninitialized data for metric calcule");

     lca = LCA_CA(get_list_ancestors(x), get_list_ancestors(y), depth);
     *lcap = lca;

     dax = min_dag_distance(g, lca, x);
     day = min_dag_distance(g, lca, y);
     drx = min_dag_distance(g, ROOT, x);
     dry = min_dag_distance(g, ROOT, y);

     return dtax(dax, day, drx, dry);
}

double sim_dtax(const struct graph *g, long x, long y)
{
     return (1.0 - dist_tax(g, x, y));
}

static inline double dps(long dax, long day, long dra)
{
     return (1.0 - ((double)dra/(dax + day + dra)));
}

double dist_ps(const struct graph *g, long x, long y)
{
     long lca, dax, day, dra;

     if (!init_metric)
	  fatal("Error, uninitialized data for metric calcule");

     lca = LCA_CA(get_list_ancestors(x), get_list_ancestors(y), depth);
     dax = min_dag_distance(g, lca, x);
     day = min_dag_distance(g, lca, y);
     dra = max_distance(g, ROOT, lca);

     return dps(dax, day, dra);
}

double dist_ps_lca(const struct graph *g, long x, long y, long *lcap)
{
     long lca, dax, day, dra;

     if (!init_metric)
	  fatal("Error, uninitialized data for metric calcule");

     lca = LCA_CA(get_list_ancestors(x), get_list_ancestors(y), depth);
     *lcap = lca;
     dax = min_dag_distance(g, lca, x);
     day = min_dag_distance(g, lca, y);
     dra = max_distance(g, ROOT, lca);

     return dps(dax, day, dra);
}

double sim_dps(const struct graph *g, long x, long y)
{
     return (1.0 - dist_ps(g, x, y));
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

void free_metric(void)
{
     for (long i = 0; i < n; i++) {
	  if (visited[i]) {
	       free(ancestors[i]->data);
	       free(ancestors[i]);
	  }
     }
     free(ancestors);
     free(depth);
     free(visited);
     free_graph(&gi);
     free_hash_map_item(&dis_map);
}

const long *get_nodes_depth(void)
{
     return depth;
}

static inline double decresing_factor(long node_depth, long max_depth)
{
     return (double)(max_depth - node_depth) / max_depth;
}

void set_max_depth(long new_depth)
{
     max_depth = new_depth;
}

double sim_str(const struct graph *g, long x, long y)
{
     double dfx, dfy;
     double sim;
  
     if (x == y) { 
	  sim  = sim_dtax(g, x, y);
     } else {
	  dfx = decresing_factor(depth[x], max_depth);
	  dfy = decresing_factor(depth[y], max_depth);    
	  sim  = sim_dtax(g, x, y) * (1.0 - MAX(dfx, dfy));
     }
     return sim;
}

double dist_str(const struct graph *g, long x, long y)
{
     return  (1.0 - sim_str(g, x, y));
}

struct long_array *lca_vector(long x, long y)
{
     if (!init_metric)
	  fatal("Error, uninitialized data for metric calcule");

     return LCA_CA_SET(get_list_ancestors(x), get_list_ancestors(y), depth);
}
