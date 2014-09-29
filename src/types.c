/**
 * Copyright (C) 2012, Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "dlist.h"
#include "memory.h"
#include "types.h"


void print_long_list(struct long_list *l)
{
  struct long_list *tmp;

  printf("Print list ");
  list_for_each_entry(tmp, &(l->list), list)
      printf("%ld ", tmp->item);
  printf("\n");
}

void destroy_long_list(struct long_list *l)
{
  struct long_list *tmp;
  struct list_head *pos, *q;

  list_for_each_safe(pos, q, &(l->list)){
    tmp = list_entry(pos, struct long_list, list);
    list_del(pos);
    free(tmp);
  }
}

void copy_long_list(struct long_list *dst, struct long_list *src)
{
  struct long_list *tmp, *new;

  list_for_each_entry(tmp, &(src->list), list) {
    new = (struct long_list *)xmalloc(sizeof(struct long_list));
    new->item = tmp->item;
    list_add_tail(&(new->list), &(dst->list));
  }
}

void print_vec_long(VEC(long) *v)
{
  unsigned long i;

  printf("\n");
  for (i = 0; i < VEC_SIZE(*v); i++) {
    printf("%ld ", VEC_GET(*v, i));
  }
  printf("\n");
}
