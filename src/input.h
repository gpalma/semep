/**
 * Copyright (C) 2013, 2014, 2015 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#ifndef ___INPUT_H
#define ___INPUT_H

struct input_data {
     long n;
     void *object;
     struct long_array anntt1;
     struct long_array anntt2;
     char **descriptions;
};

struct input_data get_input_data(const char *graph_filename, 
				 const char *desc_filename,
				 const char *annt1_filename,
				 const char *annt2_filename,
				 bool matrix, bool description);

void free_input_data(struct input_data *in, bool matrix);

#endif /* ___INPUT_H */
