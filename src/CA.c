/**
 * Copyright (C) 2012, 2014 Universidad Simón Bolívar
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

VEC(long) *get_ancestors(const struct graph *gi, long node)
{
     long i, n;
     VEC(long) *ancs;
     long *pred;
     long *discovered;
     long *finished;

     ancs = (VEC(long) *)xmalloc(sizeof(VEC(long)));
     n = gi->n_nodes;
     pred = (long *)xcalloc(n, sizeof(long));
     discovered = (long *)xcalloc(n, sizeof(long));
     finished = (long *)xcalloc(n, sizeof(long));
     VEC_INIT(long, *ancs);
     VEC_PUSH(long, *ancs, node);

     dfs_search(gi, node, discovered, finished, pred);
     for (i = 0; i < n; i++) {
	  if (pred[i] != -1) {
	       VEC_PUSH(long, *ancs, i);
	  }
     }
     free(pred);
     free(discovered);
     free(finished);

     return ancs;
}

VEC(long) **get_all_ancestors(const struct graph *g)
{
     long i, j, n;
     VEC(long) **all_a;
     long *pred;
     long *discovered;
     long *finished;

     n = g->n_nodes;
     all_a = xmalloc(n*sizeof(VEC(long)));
     pred = (long *)xcalloc(n, sizeof(long));
     discovered = (long *)xcalloc(n, sizeof(long));
     finished = (long *)xcalloc(n, sizeof(long));
     for (i = 0; i < n; i++) {
	  all_a[i] = (VEC(long) *)xmalloc(sizeof(VEC(long)));
	  VEC_INIT(long, *all_a[i]);
	  VEC_PUSH(long, *all_a[i], i);
     }
     /* start search ancestors */
     for (i = 0; i < n; i++) {
	  fprintf(stderr, "** Search in %ld  ** \n", i);
	  dfs_search(g, i, discovered, finished, pred);
	  for (j = 0; j < n; j++) {
	       if (pred[j] != -1) {
		    VEC_PUSH(long, *all_a[j], i);
	       }
	  }
     }
     free(pred);
     free(discovered);
     free(finished);

     return all_a;
}

long LCA_CA(VEC(long) *lx, VEC(long) *ly, long *depth)
{
     long i, j, nlx, nly, vx;
     long max, lca;

     max = -1;
     lca = -1;
     nlx = VEC_SIZE(*lx);
     nly = VEC_SIZE(*ly);
     for (i = 0; i < nlx; i++) {
	  vx = VEC_GET(*lx, i);
	  for (j = 0; j < nly; j++) {
	       if ((vx == VEC_GET(*ly, j)) && (depth[vx] > max)) {
		    max = depth[vx];
		    lca = vx;
	       }
	  }
     }
     if (lca == -1)
	  fatal("Error with the lowest common ancestor");
     return lca;
}

VEC(long) *LCA_CA_SET(VEC(long) *lx, VEC(long) *ly, long *depth)
{
     long i, j, nlx, nly, vx, max, lcam;
     VEC(long) *lca;

     lca = (VEC(long) *)xmalloc(sizeof(VEC(long)));
     VEC_INIT(long, *lca);

     max = -1;
     lcam = -1;
     nlx = VEC_SIZE(*lx);
     nly = VEC_SIZE(*ly);
     for (i = 0; i < nlx; i++) {
	  vx = VEC_GET(*lx, i);
	  for (j = 0; j < nly; j++) {
	       if ((vx == VEC_GET(*ly, j)) && (depth[vx] > max)) {
		    max = depth[vx];
		    lcam = vx;
	       }
	  }
     }

     if (lcam == -1)
	  fatal("Error with the lowest common ancestor");

     for (i = 0; i < nlx;  i++) {
	  vx = VEC_GET(*lx,  i);
	  for (j = 0; j < nly; j++) {
	       if ((vx == VEC_GET(*ly, j)) && (depth[vx] == max)) {
		    VEC_PUSH(long, *lca, vx);
	       }
	  }
     }
     return lca;
}
