/**
 * Copyright (C) 2014, 2015 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#ifndef ___SEMEP_H
#define ___SEMEP_H

double annotation_partition(void *object, long n_nodes,
			    const struct long_array *v1, const struct long_array *v2,
			    double threshold_E1, double threshold_E2,
			    double threshold_bt, const char name1[], const char name2[],
			    char **desc, bool prediction, bool matrix, enum measure d);

#endif /* ___SEMEP_H */
