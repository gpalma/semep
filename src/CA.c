/**
 * Copyright (C) 2012, 2014, 2015 Universidad Simón Bolívar
 *
 * @brief Common ancestors
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <limits.h>

#include "dlist.h"
#include "types.h"
#include "util.h"
#include "memory.h"
#include "graph.h"
#include "CA.h"

#define INVALID -1

struct long_array *get_ancestors(const struct graph *gi, long node)
{
     long i, n;
     struct long_array *ancs;
     long *pred;
     long *discovered;
     long *finished;

     ALLOC_STRUCT(ancs);
     n = gi->n_nodes;
     pred = (long *)xcalloc(n, sizeof(long));
     discovered = (long *)xcalloc(n, sizeof(long));
     finished = (long *)xcalloc(n, sizeof(long));
     
     ARRAY_PUSH(*ancs, node);
     dfs_search(gi, node, discovered, finished, pred);
     for (i = 0; i < n; i++) {
	  if (pred[i] != -1) {
	       ARRAY_PUSH(*ancs, i);
	  }
     }
     free(pred);
     free(discovered);
     free(finished);

     return ancs;
}

struct long_array **get_all_ancestors(const struct graph *g)
{
     long i, j, n;
     struct long_array **all_a;
     long *pred;
     long *discovered;
     long *finished;
     struct long_array aux = {0, 0, NULL};

     n = g->n_nodes;
     all_a = xcalloc(n, sizeof(struct long_array *));
     
     pred = (long *)xcalloc(n, sizeof(long));
     discovered = (long *)xcalloc(n, sizeof(long));
     finished = (long *)xcalloc(n, sizeof(long));
     for (i = 0; i < n; i++) {
	  all_a[i] = xmalloc(sizeof(struct long_array));
	  *all_a[i] = aux;
	  ARRAY_PUSH(*all_a[i], i);
     }
     /* start search ancestors */
     for (i = 0; i < n; i++) {
	  fprintf(stderr, "** Search in %ld  ** \n", i);
	  dfs_search(g, i, discovered, finished, pred);
	  for (j = 0; j < n; j++) {
	       if (pred[j] != -1) {
		    ARRAY_PUSH(*all_a[i], i);
	       }
	  }
     }
     free(pred);
     free(discovered);
     free(finished);

     return all_a;
}

long LCA_CA(struct long_array *lx, struct long_array *ly, long *depth)
{
     long i, j, nlx, nly, vx;
     long max, lca;

     max = INVALID;
     lca = INVALID;
     nlx = lx->nr;
     nly = ly->nr;
     for (i = 0; i < nlx; i++) {
	  vx = lx->data[i];
	  for (j = 0; j < nly; j++) {
	       if ((vx == ly->data[j]) && (depth[vx] > max)) {
		    max = depth[vx];
		    lca = vx;
	       }
	  }
     }
     if (lca == -1)
	  fatal("Error with the lowest common ancestor");
     return lca;
}

struct long_array *LCA_CA_SET(struct long_array *lx, 
			      struct long_array *ly, long *depth)
{
     long i, j, nlx, nly, vx, max, lcam;
     struct long_array *lca;

     ALLOC_STRUCT(lca);
     max = INVALID;
     lcam = INVALID;
     nlx = lx->nr;
     nly = ly->nr;
     for (i = 0; i < nlx; i++) {
	  vx = lx->data[i];
	  for (j = 0; j < nly; j++) {
	       if ((vx == ly->data[j]) && (depth[vx] > max)) {
		    max = depth[vx];
		    lcam = vx;
	       }
	  }
     }

     if (lcam == -1)
	  fatal("Error with the lowest common ancestor");

     for (i = 0; i < nlx;  i++) {
	  vx = lx->data[i];
	  for (j = 0; j < nly; j++) {
	       if ((vx == ly->data[j]) && (depth[vx] == max)) {
		    ARRAY_PUSH(*lca, vx);
	       }
	  }
     }
     return lca;
}
