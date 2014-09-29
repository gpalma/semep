/**
 * Copyright (C) 2012, 2013 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#ifndef ___TYPES_H
#define ___TYPES_H

#include "vec.h"
#include "dlist.h"

/**
 * Structures that define a vector
 */
DEFINE_VEC(long);

typedef VEC(long) *VEC_LONG_PTR;
DEFINE_VEC(VEC_LONG_PTR);

typedef char *string;
DEFINE_VEC(string);

/**
 * Lists of basic types
 */
struct long_list {
     long item;
     struct list_head list;
};

/**
 * Long Pairs
 */
struct lpairs {
     long x;
     long y;
};

/**
 * Solver output
 */
enum output {
     MFILE,
     DRAW
};

/**
 * Measure
 */
enum measure {
     DTAX,
     DSTR,
     DPS
};

void print_long_list(struct long_list *l) ;

void destroy_long_list(struct long_list *l);

void copy_long_list(struct long_list *dst, struct long_list *src);

void print_vec_long(VEC(long) *v);

#endif /* ___TYPES_H */
