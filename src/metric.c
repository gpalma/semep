/**
 * Copyright (C) 2012, 2013, 2014 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#include "dlist.h"
#include "types.h"
#include "util.h"
#include "memory.h"
#include "graph.h"
#include "CA.h"
#include "metric.h"

#define ROOT  0

static long n;
static bool *visited;
static VEC(long) **ancestors;
static long *depth;
static struct graph gi;
static bool init_metric = false;
static long max_depth;

void init_metric_data(const struct graph *g)
{
  depth = calculate_depth(g);
  DEBUG("\n** Depth node calculation done ** \n");
  n = g->n_nodes;
  ancestors = xmalloc(n*sizeof(VEC(long)));
  graph_inverse(g, &gi);
  DEBUG("\n** Graph inverse done ** \n");
  visited = xcalloc(n, sizeof(bool));
  max_depth = INT_MAX;
  init_metric = true;
}

static inline double dtax(long dax, long day, long drx, long dry)
{
  /*     fprintf(stderr, "\n dtax %ld %ld %ld %ld\n", dax, day, drx, dry);*/
  double r = (double)(dax+day)/(drx+dry);

  return MIN(r, 1.0);
}

static VEC(long) *get_list_ancestors(long node)
{
  VEC(long) *la;

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
  VEC(long) *lx, *ly;

  if (!init_metric)
    fatal("Error, uninitialized data for metric calcule");

  lx = get_list_ancestors(x);
  ly = get_list_ancestors(y);
  lca = LCA_CA(lx, ly, depth);
  dax = min_distance(g, lca, x);
  day = min_distance(g, lca, y);
  drx = min_distance(g, ROOT, x);
  dry = min_distance(g, ROOT, y);

  return dtax(dax, day, drx, dry);
}

double dist_tax_lca(const struct graph *g, long x, long y, long *lcap)
{
  long lca, dax, day, drx, dry;
  VEC(long) *lx, *ly;

  if (!init_metric)
    fatal("Error, uninitialized data for metric calcule");

  lx = get_list_ancestors(x);
  ly = get_list_ancestors(y);
  lca = LCA_CA(lx, ly, depth);
  *lcap = lca;
  dax = min_distance(g, lca, x);
  day = min_distance(g, lca, y);
  drx = min_distance(g, ROOT, x);
  dry = min_distance(g, ROOT, y);

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
  VEC(long) *lx, *ly;

  if (!init_metric)
    fatal("Error, uninitialized data for metric calcule");

  lx = get_list_ancestors(x);
  ly = get_list_ancestors(y);
  lca = LCA_CA(lx, ly, depth);
  dax = min_distance(g, lca, x);
  day = min_distance(g, lca, y);
  dra = max_distance(g, ROOT, lca);

  return dps(dax, day, dra);
}

double dist_ps_lca(const struct graph *g, long x, long y, long *lcap)
{
  long lca, dax, day, dra;
  VEC(long) *lx, *ly;

  if (!init_metric)
    fatal("Error, uninitialized data for metric calcule");

  lx = get_list_ancestors(x);
  ly = get_list_ancestors(y);
  lca = LCA_CA(lx, ly, depth);
  *lcap = lca;
  dax = min_distance(g, lca, x);
  day = min_distance(g, lca, y);
  dra = max_distance(g, ROOT, lca);

  return dps(dax, day, dra);
}

double sim_dps(const struct graph *g, long x, long y)
{
  return (1.0 - dist_ps(g, x, y));
}

void free_metric(void)
{
  long i;

  for (i = 0; i < n; i++) {
    if (visited[i]) {
      VEC_DESTROY(*ancestors[i]);
      free(ancestors[i]);
    }
  }
  free(ancestors);
  free(depth);
  free(visited);
  free_graph(&gi);
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

const long *get_nodes_depth(void)
{
  return depth;
}

VEC(long) *lca_vector(long x, long y)
{
  VEC(long) *lx, *ly;

  if (!init_metric)
    fatal("Error, uninitialized data for metric calcule");

  lx = get_list_ancestors(x);
  ly = get_list_ancestors(y);

  return LCA_CA_SET(lx, ly, depth);
}
