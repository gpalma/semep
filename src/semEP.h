/**
 * Copyright (C) 2013-2015 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#ifndef ___SEMEP_H
#define ___SEMEP_H

#include "types.h"

struct color {
     int id;
     double sim_entity_1;
     double sim_entity_2;
     double sim_between;
     struct int_array id_nodes;
     struct int_array entities1;
     struct int_array entities2;
};

struct node {
     int id;
     int pos1;
     int pos2;
     double sim;
     struct color *cp;
};

struct node_ptr_array{
     unsigned nr;
     unsigned alloc;
     struct node **data;
};

double semEP_solver(const struct matrix *lmatrix, const struct matrix *rmatrix,
		    const struct int_array *lterms, const struct int_array *rterms,
		    const struct string_array *desc, struct node_ptr_array *color_nodes,
		    double lthreshold, double rthreshold, const char *bpgraph_name);

#endif /* ___SEMEP_H */
