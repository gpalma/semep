/**
 * Copyright (C) 2012, 2013 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#ifndef ___UTIL_H
#define ___UTIL_H

#ifdef PRGDEBUG
#define DEBUG(...) printf(__VA_ARGS__)
#else
#define DEBUG(...)
#endif

#define MIN(x, y) ({                            \
      typeof(x) _min1 = (x);			\
      typeof(y) _min2 = (y);			\
      (void) (&_min1 == &_min2);		\
      _min1 < _min2 ? _min1 : _min2; })

#define MAX(x, y) ({                            \
      typeof(x) _max1 = (x);			\
      typeof(y) _max2 = (y);			\
      (void) (&_max1 == &_max2);		\
      _max1 > _max2 ? _max1 : _max2; })

#define SWAP(a, b)				\
  do { typeof(a) __tmp = (a);                   \
    (a) = (b);                                  \
    (b) = __tmp;				\
  } while (0)

#define ABS(val)      ((val) < 0 ? -(val):(val))

#define ISEVEN(x)     (!((x)&0x01))

void fatal(const char *msg, ...);

int error(const char *msg, ...);

int **int_matrix(long nrl, long nrh, long ncl, long nch);

void free_int_matrix(int **m, long nrl, long ncl);

long **long_matrix(long nrl, long nrh, long ncl, long nch);

void free_long_matrix(long **m, long nrl, long ncl);

char *ltostr(char *buf, long number, size_t digits);

long number_digits(long n);

char **new_cmatrix(long n, long m);

double **double_matrix(long nrl, long nrh, long ncl, long nch);

void free_double_matrix(double **m, long nrl, long ncl);

#endif /* ___UTIL_H */
