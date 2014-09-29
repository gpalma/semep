/**
 * Copyright (C) 2011, Universidad Simón Bolívar
 *
 * @brief Wrappers for the C libraries allocation functions
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

static inline void *xmalloc(size_t size)
{
  void *ptr = malloc(size);
  if (!ptr) {
    fprintf(stderr, "Out of memory, malloc failed tried to allocate %lu bytes\n",
            (unsigned long)size);
    exit(1);
  }
  return ptr;
}

static inline void *xrealloc(void *old_ptr, size_t size)
{
  void *ptr = realloc(old_ptr, size);
  if (!ptr) {
    fprintf(stderr, "Out of memory, realloc\n");
    exit(1);
  }
  return ptr;
}

static inline void *xcalloc(size_t nmemb, size_t size)
{
  void *ptr = calloc(nmemb, size);
  if (!ptr) {
    fprintf(stderr, "Out of memory, calloc failed\n");
    exit(1);
  }
  return ptr;
}
