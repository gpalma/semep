/**
 * Copyright (C) 2012, 2015 Universidad Simón Bolívar
 *
 * @brief Common ancestors
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#ifndef ___CA_H
#define ___CA_H

struct long_array *get_ancestors(const struct graph *gi, long node);

struct long_array **get_all_ancestors(const struct graph *g);

long LCA_CA(struct long_array *lx, struct long_array *ly, long *depth);

struct long_array *LCA_CA_SET(struct long_array *lx, struct long_array *ly, long *depth);

#endif /* ___CA_H */
