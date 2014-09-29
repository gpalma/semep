/**
 * Copyright (C) 2013, 2014 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#ifndef ___INPUT_H
#define ___INPUT_H

struct input_data {
     void *object;
     long n;
     VEC(long) anntt1, anntt2;
     char **descriptions;
};

struct input_data get_input_data(const char *graph_filename, const char *desc_filename,
				 const char *annt1_filename, const char *annt2_filename,
				 bool matrix, bool description);

void free_input_data(struct input_data *in, bool matrix);

#endif /* ___INPUT_H */
