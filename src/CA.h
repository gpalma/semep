/**
 * Copyright (C) 2012, Universidad Simón Bolívar
 *
 * @brief Common ancestors
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#ifndef ___CA_H
#define ___CA_H

VEC(long) *get_ancestors(const struct graph *gi, long node);

VEC(long) **get_all_ancestors(const struct graph *g);

long LCA_CA(VEC(long) *lx, VEC(long) *ly, long *depth);

VEC(long) *LCA_CA_SET(VEC(long) *lx, VEC(long) *ly, long *depth);

#endif /* ___CA_H */
