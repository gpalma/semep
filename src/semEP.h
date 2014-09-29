/**
 * Copyright (C) 2014, Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#ifndef ___SEMEP_H
#define ___SEMEP_H

double annotation_partition(void *object, long n_nodes,
			    const VEC(long) *v1, const VEC(long) *v2,
			    double threshold1, double threshold2, double threshold3,
			    const char name1[], const char name2[], char **desc, 
			    bool prediction, bool in_matrix, enum measure d);

#endif /* ___SEMEP_H */
