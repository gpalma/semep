/**
 * Copyright (C) 2011, 2014 Universidad Simón Bolívar
 *
 * @brief Dynamic array library
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 * @version 0.6
 */

#ifndef __VEC_H
#define __VEC_H

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define __MIN_ALLOC     10
#define __MIN_SIZE      55

#define VEC(type) VEC_##type

#define DEFINE_VEC(type)			\
     typedef struct VEC(type) {			\
	  size_t nr;				\
	  size_t alloc;				\
	  type *data;				\
     } VEC(type)

#define VEC_INIT(type, vec)						\
     do {								\
	  (vec).data = (type *)malloc(__MIN_ALLOC * sizeof(type));	\
	  (vec).alloc = __MIN_ALLOC;					\
	  (vec).nr = 0;							\
     } while(0)

#define VEC_INIT_N(type, vec, n)				\
     do {							\
	  (vec).data = (type *)malloc((n) * sizeof(type));	\
	  (vec).alloc = (n);					\
	  (vec).nr = 0;						\
     } while(0)

#define VEC_INIT_EMPTY_ZERO(type, vec)					\
     do {								\
     size_t _alloc = sizeof(__MIN_ALLOC * sizeof(type);			\
			    (vec).data = (type *)malloc(_alloc);	\
			    memset((vec).data, 0, _alloc);		\
			    (vec).alloc = __MIN_ALLOC;			\
			    (vec).nr = 0;				\
			    } while (0)

#define VEC_INIT_EMPTY_ZERO_N(type, vec, n)				\
     do {								\
     size_t _alloc = sizeof((n) * sizeof(type);				\
			    (vec).data = (type *)malloc(_alloc);	\
			    memset((vec).data, 0, _alloc);		\
			    (vec).alloc = __MIN_ALLOC;			\
			    (vec).nr = 0;				\
			    } while (0)

#define VEC_EMPTY(vec) ((vec).nr == 0 ? 1: 0)

#define VEC_GET(vec, index) ((vec).data[(index)])

#define VEC_GET_SAFE(vec, index) ((vec).nr <= (index) ?			\
				  (fprintf(stderr,"Fatal array bound exceeded in file %s line %d\n", __FILE__, __LINE__), \
				   exit(1), 1)				\
				  : (vec).data[(index)])

#define VEC_FRONT(vec) ((vec).data[0])

#define VEC_LAST(vec) ((vec).data[(vec).nr-1])

#define VEC_DESTROY(vec)			\
     do {					\
	  free((vec).data);			\
	  (vec).alloc = 0;			\
	  (vec).nr = 0;				\
     } while (0)

#define VEC_POP(vec) ((vec).data[--((vec).nr)])

#define VEC_SIZE(vec) ((vec).nr)

#define VEC_LENGTH(vec) ((vec).alloc)

#define VEC_CLEAR(vec) ((vec).nr = 0)

#define _ALLOC_NR(x) (((x)+16)*3/2)

#define VEC_PUSH(type, vec, object)					\
     do {								\
	  size_t _nr = (vec).nr;					\
	  size_t _alloc = (vec).alloc;					\
	  if (_nr == _alloc) {						\
	       _alloc = _ALLOC_NR(_alloc);				\
	       (vec).data = (type *)realloc((vec).data, _alloc * sizeof(type)); \
	       (vec).alloc = _alloc;					\
	  }								\
	  (vec).data[_nr] = (object);					\
	  (vec).nr++;							\
     } while (0)

#define VEC_PUSH_FAST(vec, object) (vec).data[(vec).nr++] = (object)

#define VEC_INSERT(type, vec, object, index)				\
     do {								\
	  size_t _alloc = (vec).alloc;					\
	  size_t _nr = (vec).nr;					\
	  if (_nr < (index)) {						\
	       fprintf(stderr,"Error array bound exceeded in file %s line %d\n", __FILE__, __LINE__); \
	       exit(1);							\
	  }								\
	  if (_nr == _alloc) {						\
	       _alloc = _ALLOC_NR(_alloc);				\
	       (vec).data = (type *)realloc((vec).data, _alloc * sizeof(type)); \
	       (vec).alloc = _alloc;					\
	  }								\
	  if (_nr > (index))						\
	       memmove(((vec).data+(index)+1), ((vec).data+(index)), (_nr-(index))*sizeof(*((vec).data))); \
	  (vec).data[(index)] = (object);				\
	  (vec).nr++;							\
     } while (0)

#define VEC_SET(vec, index, object) (vec).data[(index)] = (object)

#define VEC_SET_SAFE(vec, index, object)				\
     do {								\
	  if ((vec).nr <= (index)) {					\
	       fprintf(stderr,"Error fatal array bound exceeded in file %s line %d -- i %u -- v %d \n", __FILE__, __LINE__, (index), (object)); \
	       exit(1);							\
	  }								\
	  (vec).data[(index)] = (object);				\
     } while (0)

#define VEC_COPY(type, from, to)					\
     do {								\
	  size_t _alloc = VEC_LENGTH((to));				\
	  size_t _size = VEC_SIZE((from));				\
	  VEC_CLEAR((to));						\
	  if (_alloc < _size) {						\
	       (to).data = (type *) realloc((to).data, _size * sizeof(type)); \
	       (to).alloc = _size;					\
	  }								\
	  memcpy((to).data, (from).data, _size*sizeof(type));		\
	  (to).nr = _size;						\
     } while (0)

#define VEC_ERASE(vec, index)						\
     do {								\
	  size_t _nr = (vec).nr;					\
	  if (_nr <= (index)) {						\
	       fprintf(stderr,"Error fatal array bound exceeded in file %s line %d\n", __FILE__, __LINE__); \
	       exit(1);							\
	  }								\
	  memmove(((vec).data+(index)), ((vec).data+(index)+1), (_nr-1-(index))*sizeof(*((vec).data))); \
	  (vec).nr--;							\
     } while (0)

#define VEC_RESIZE(type, vec, n)					\
     do {								\
	  if ((vec).nr > (n))						\
	       (vec).nr = (n);						\
	  (vec).alloc = (n);						\
	  (vec).data = (type *) realloc((vec).data, (n) * sizeof(type)); \
     } while(0)

#define VEC_EXCH(type, a, b)			\
     do {					\
	  type _aux = (a);			\
	  (a) = (b);				\
	  (b) = _aux;				\
     } while(0)

#define __VEC_SORT2(vec, type, _cmp)			\
     do {						\
	  type _tmps2;					\
	  if (_cmp((vec)[0], (vec)[1]) <= 0) break;	\
	  _tmps2 = (vec)[0];				\
	  (vec)[0] = (vec)[1];				\
	  (vec)[1] = _tmps2;				\
     } while(0)

#define __VEC_SORT3(vec, type, _cmp)				\
     do {							\
	  type _tmps3;						\
	  if (_cmp((vec)[0], (vec)[1]) <= 0) {			\
	       if (_cmp((vec)[1], (vec)[2]) <= 0) break;	\
	       if (_cmp((vec)[2], (vec)[0]) <= 0) {		\
		    _tmps3 = (vec)[0];				\
		    (vec)[0] = (vec)[2];			\
		    (vec)[2] = (vec)[1];			\
		    (vec)[1] = _tmps3;				\
		    break;					\
	       }						\
	       _tmps3 = (vec)[1];				\
	  } else {						\
	       _tmps3 = (vec)[0];				\
	       if (_cmp((vec)[0], (vec)[2]) <= 0) {		\
		    (vec)[0] = (vec)[1];			\
		    (vec)[1] = _tmps3;				\
		    break;					\
	       }						\
	       if (_cmp((vec)[2], (vec)[1]) <= 0) {		\
		    (vec)[0] = (vec)[2];			\
		    (vec)[2] = _tmps3;				\
		    break;					\
	       }						\
	       (vec)[0] = (vec)[1];				\
	  }							\
	  (vec)[1] = (vec)[2];					\
	  (vec)[2] = _tmps3;					\
     } while(0)

#define __VEC_SORT4(vec, type, _cmp)					\
     do {								\
	  type _tmps4;							\
	  if (_cmp((vec)[0], (vec)[1]) < 0) {				\
	       if (_cmp((vec)[1], (vec)[2]) < 0) {			\
		    if (_cmp((vec)[1], (vec)[3]) < 0) {			\
			 if (_cmp((vec)[2], (vec)[3]) >= 0) {		\
			      _tmps4 = (vec)[2];			\
			      (vec)[2] = (vec)[3];			\
			      (vec)[3] = _tmps4;			\
			 }						\
		    } else {						\
			 _tmps4 = (vec)[1];				\
			 if (_cmp((vec)[0], (vec)[3]) < 0) {		\
			      (vec)[1] = (vec)[3];			\
			 } else {					\
			      (vec)[1] = (vec)[0];			\
			      (vec)[0] = (vec)[3];			\
			 }						\
			 (vec)[3] = (vec)[2];				\
			 (vec)[2] = _tmps4;				\
		    }							\
	       } else {							\
		    if (_cmp((vec)[0], (vec)[2]) < 0) {			\
			 if (_cmp((vec)[2], (vec)[3]) < 0) {		\
			      if (_cmp((vec)[1], (vec)[3]) < 0) {	\
				   _tmps4 = (vec)[1];			\
			      } else {					\
				   _tmps4 = (vec)[3];			\
				   (vec)[3] = (vec)[1];			\
			      }						\
			      (vec)[1] = (vec)[2];			\
			      (vec)[2] = _tmps4;			\
			 } else {					\
			      if (_cmp((vec)[0], (vec)[3]) < 0) {	\
				   _tmps4 = (vec)[3];			\
			      } else {					\
				   _tmps4 = (vec)[0];			\
				   (vec)[0] = (vec)[3];			\
			      }						\
			      (vec)[3] = (vec)[1];			\
			      (vec)[1] = _tmps4;			\
			 }						\
		    } else {						\
			 if (_cmp((vec)[0], (vec)[3]) < 0) {		\
			      _tmps4 = (vec)[0];			\
			      (vec)[0] = (vec)[2];			\
			      if (_cmp((vec)[1], (vec)[3]) < 0) {	\
				   (vec)[2] = (vec)[1];			\
			      } else {					\
				   (vec)[2] = (vec)[3];			\
				   (vec)[3] = (vec)[1];			\
			      }						\
			      (vec)[1] = _tmps4;			\
			 } else {					\
                           if (_cmp((vec)[2], (vec)[3]) < 0) {		\
				   _tmps4 = (vec)[0];			\
				   (vec)[0] = (vec)[2];			\
				   (vec)[2] = _tmps4;			\
				   _tmps4 = (vec)[1];			\
				   (vec)[1] = (vec)[3];			\
			      } else {					\
				   _tmps4 = (vec)[1];			\
				   (vec)[1] = (vec)[2];			\
				   (vec)[2] = (vec)[0];			\
				   (vec)[0] = (vec)[3];			\
			      }						\
			      (vec)[3] = _tmps4;			\
			 }						\
		    }							\
	       }							\
	  } else {							\
	       _tmps4 = (vec)[0];					\
	       if (_cmp(_tmps4, (vec)[2]) < 0) {			\
		    if (_cmp(_tmps4, (vec)[3]) < 0) {			\
			 (vec)[0] = (vec)[1];				\
			 (vec)[1] = _tmps4;				\
			 if (_cmp((vec)[2], (vec)[3]) >= 0) {		\
			      _tmps4 = (vec)[2];			\
			      (vec)[2] = (vec)[3];			\
			      (vec)[3] = _tmps4;			\
			 }						\
		    } else {						\
			 if (_cmp((vec)[1], (vec)[3]) < 0) {		\
			      (vec)[0] = (vec)[1];			\
			      (vec)[1] = (vec)[3];			\
			 } else {					\
			      (vec)[0] = (vec)[3];			\
			 }						\
			 (vec)[3] = (vec)[2];				\
			 (vec)[2] = _tmps4;				\
		    }							\
	       } else {							\
		    if (_cmp((vec)[1], (vec)[2]) < 0) {			\
			 if (_cmp((vec)[2], (vec)[3]) < 0) {		\
			      (vec)[0] = (vec)[1];			\
			      (vec)[1] = (vec)[2];			\
			      if (_cmp(_tmps4, (vec)[3]) < 0) {		\
				   (vec)[2] = _tmps4;			\
			      } else {					\
				   (vec)[2] = (vec)[3];			\
				   (vec)[3] = _tmps4;			\
			      }						\
			 } else {					\
			      if (_cmp((vec)[1], (vec)[3]) < 0) {	\
				   (vec)[0] = (vec)[1];			\
				   (vec)[1] = (vec)[3];			\
			      } else {					\
				   (vec)[0] = (vec)[3];			\
			      }						\
			      (vec)[3] = _tmps4;			\
			 }						\
		    } else {						\
			 if (_cmp((vec)[1],(vec)[3]) < 0) {		\
			      (vec)[0] = (vec)[2];			\
			      if (_cmp(_tmps4, (vec)[3]) < 0) {		\
				   (vec)[2] = _tmps4;			\
			      } else {					\
				   (vec)[2] = (vec)[3];			\
				   (vec)[3] = _tmps4;			\
			      }						\
			 } else {					\
			      if (_cmp((vec)[2],(vec)[3]) < 0) {	\
				   (vec)[0] = (vec)[2];			\
				   (vec)[2] = (vec)[1];			\
				   (vec)[1] = (vec)[3];			\
				   (vec)[3] = _tmps4;			\
			      } else {					\
				   (vec)[0] = (vec)[3];			\
				   (vec)[3] = _tmps4;			\
				   _tmps4 = (vec)[1];			\
				   (vec)[1] = (vec)[2];			\
				   (vec)[2] = _tmps4;			\
			      }						\
			 }						\
		    }							\
	       }							\
	  }								\
     } while(0)

#define __VEC_INSERTION_SORT_AUX(vec, type, left, right, _cmp)		\
     do {								\
     size_t _loc;                                                       \
     int _i, __n = 0;                                                   \
     type __data;                                                       \
	  __n = right - left + 1;					\
	  if (__n > 4) {						\
	       for (_loc = (left)+1; _loc <= (right); _loc++) {		\
		    _i = _loc-1;					\
		    __data = (vec)[_loc];				\
		    while (_i >= 0 && _cmp((vec)[_i], __data) > 0) {	\
			 (vec)[_i+1] = (vec)[_i];			\
			 _i--;						\
		    }							\
		    (vec)[_i+1] = __data;				\
	       }							\
	  } else if (__n == 4) {					\
	       __VEC_SORT4(vec+left, type, _cmp);			\
	  } else if (__n == 3) {					\
	       __VEC_SORT3(vec+left, type, _cmp);			\
	  } else if (__n == 2) {					\
	       __VEC_SORT2(vec+left, type, _cmp);			\
	  } else {							\
	       break;							\
	  }								\
     } while(0)

#define VEC_INSERTION_SORT(type, vec, _cmp) __VEC_INSERTION_SORT_AUX((vec).data, type, 0, VEC_SIZE((vec))-1, _cmp)

#define VEC_INIT_QSORT(type, name, _cmp)				\
     void vec_sort_aux_##name(type vec[], size_t left, size_t right) {	\
	  size_t i, j, mid, pos, l1, l2;				\
	  type pivot;							\
	  while (left+1 < right) {					\
	       mid = (left+right)/2;					\
	       if (_cmp(vec[left], vec[mid]) < 0) {			\
		    if (_cmp(vec[mid], vec[right-1]) < 0) {		\
			 pos = mid;					\
		    } else {						\
			 if (_cmp(vec[left], vec[right-1]) < 0) {	\
			      pos = right-1;				\
			 } else {					\
			      pos = left;				\
			 }						\
		    }							\
	       } else  {						\
                 if (_cmp(vec[left], vec[right-1]) < 0) {              \
			 pos = left;					\
		    } else {						\
			 if (_cmp(vec[mid], vec[right-1]) < 0) {	\
			      pos = right-1;				\
			 } else {					\
			      pos = mid;				\
			 }						\
		    }							\
	       }							\
	       pivot = vec[pos];					\
	       i = left;						\
	       j = right;						\
	       while (i < j) {						\
		    while (_cmp(vec[i], pivot) < 0)			\
			 i++;						\
		    while (_cmp(pivot, vec[j-1]) < 0)			\
			 j--;						\
		    if (i != j) {					\
			 VEC_EXCH(type, vec[i], vec[j-1]);		\
			 i++;						\
			 j--;						\
		    }							\
	       }							\
	       l1 = j-left;						\
	       l2 = right-i;						\
	       if (l1 <= l2) {						\
		    if (l1 <= __MIN_SIZE) {				\
			 __VEC_INSERTION_SORT_AUX(vec, type, left, j-1, _cmp); \
		    } else {						\
			 vec_sort_aux_##name(vec, left, j);		\
		    }							\
		    left = i;						\
	       } else {							\
		    if (l2 <= __MIN_SIZE) {				\
			 __VEC_INSERTION_SORT_AUX(vec, type, i, right-1, _cmp); \
		    } else {						\
			 vec_sort_aux_##name(vec, i, right);		\
		    }							\
		    right = j;						\
	       }							\
	  }								\
     }

#define VEC_QSORT(vec, name)						\
     do {								\
	  if (VEC_SIZE((vec)) >= 2)					\
	       vec_sort_aux_##name((vec).data, 0, VEC_SIZE((vec)));	\
     } while(0)


#endif /* VEC_H */
