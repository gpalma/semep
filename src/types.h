/**
 * Copyright (C) 2012, 2013, 2015 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#ifndef ___TYPES_H
#define ___TYPES_H

#include "dlist.h"

struct long_array {
     unsigned nr;
     unsigned alloc;
     long *data;
};

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

void print_vec_long(struct long_array *v);

#endif /* ___TYPES_H */
