/**
 * Copyright (C) 2012 Universidad Simon Bolivar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>

#include "memory.h"
#include "util.h"

#define BUF         1024
#define END         1
#define BUFSIZE1    256

static void write_out(const char *prefix, const char *msg, va_list params)
{
     char buf[BUF];
     vsnprintf(buf, sizeof(buf), msg, params);
     fprintf(stderr,"%s%s\n", prefix, buf);
}

void fatal(const char *msg, ...)
{
     va_list params;
     
     va_start(params, msg);
     write_out("Fatal: ", msg, params);
     va_end(params);
     
     exit(1);
}

int error(const char *msg, ...)
{
     va_list params;

     va_start(params, msg);
     write_out("Error: ", msg, params);
     va_end(params);

     return -1;
}

/*
  Allocate a int matrix with
  subscript range m[nrl..nrh][ncl..nch]
*/
int **int_matrix(long nrl, long nrh, long ncl, long nch)

{
     long i, nrow, ncol;
     int **m;

     nrow = nrh-nrl+1;
     ncol = nch-ncl+1;

     /* allocate pointers to rows */
     m = (int **) xmalloc((size_t)((nrow+END)*sizeof(int*)));
     
     m += END;
     m -= nrl;

     /* allocate rows and set pointers to them */
     m[nrl] = (int *) xmalloc((nrow*ncol+END)*sizeof(int));
     m[nrl] += END;
     m[nrl] -= ncl;
     
     for( i = nrl+1; i <= nrh; i++)
	  m[i] = m[i-1] + ncol;

     /* return pointer to array of pointers to rows */
     return m;
}


void free_int_matrix(int **m, long nrl, long ncl)
{
     free((m[nrl]+ncl-END));
     free((m+nrl-END));
}

/*
  Allocate a long matrix with
  subscript range m[nrl..nrh][ncl..nch]
*/
long **long_matrix(long nrl, long nrh, long ncl, long nch)
     
{
     long i, nrow, ncol;
     long **m;
     
     nrow = nrh-nrl+1;
     ncol = nch-ncl+1;
     
     /* allocate pointers to rows */
     m = (long **) xmalloc((size_t)((nrow+END)*sizeof(long*)));

     m += END;
     m -= nrl;

     /* allocate rows and set pointers to them */
     m[nrl] = (long *) xmalloc((size_t)((nrow*ncol+END)*sizeof(long)));
     m[nrl] += END;
     m[nrl] -= ncl;

     for( i = nrl+1; i <= nrh; i++)
	  m[i] = m[i-1] + ncol;

     /* return pointer to array of pointers to rows */
     return m;
}

void free_long_matrix(long **m, long nrl, long ncl)
{
     free((m[nrl]+ncl-END));
     free((m+nrl-END));
}

char *ltostr(char *buf, long number, size_t digits)
{
     long i, nz;
     char bufn[BUFSIZE1];

     sprintf(bufn,"%ld", number);
     nz = digits - strlen(bufn);
     for (i = 0; i < nz; i++)
	  buf[i] = 48;
     buf[i] = '\0';

     return strcat(buf, bufn);
}

long number_digits(long n)
{
     long cont;

     cont = 0;
     if (n == 0) {
	  cont = 1;
     } else {
	  while (n > 0) {
	       cont++;
	       n /= 10;
	  }
     }
     return cont;
}

char **new_cmatrix(long n, long m)
{
     long i;
     char **cmatrix;

     cmatrix = (char **)xmalloc(sizeof(char) * n * m + sizeof(char *) * n);
     for (i = 0; i < n; i++)
	  cmatrix[i] = (char *)(cmatrix + n) + i * m;

     return cmatrix;
}


/*
  Allocate a double matrix with
  subscript range m[nrl..nrh][ncl..nch]
*/
double **double_matrix(long nrl, long nrh, long ncl, long nch)

{
     long i, nrow, ncol;
     double **m;

     nrow = nrh-nrl+1;
     ncol = nch-ncl+1;

     /* allocate pointers to rows */
     m = (double **) xmalloc((size_t)((nrow+END)*sizeof(double*)));

     m += END;
     m -= nrl;

     /* allocate rows and set pointers to them */
     m[nrl] = (double *) xmalloc((size_t)((nrow*ncol+END)*sizeof(double)));
     m[nrl] += END;
     m[nrl] -= ncl;

     for( i = nrl+1; i <= nrh; i++)
	  m[i] = m[i-1] + ncol;

     /* return pointer to array of pointers to rows */
     return m;
}


void free_double_matrix(double **m, long nrl, long ncl)
{
     free((m[nrl]+ncl-END));
     free((m+nrl-END));
}
