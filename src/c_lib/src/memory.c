/*
Wrappers for malloc(), calloc(), free(), realloc().


AUTHORS:

- Jeroen Demeyer (2011-01-13): initial version (#10258)

*/

/*****************************************************************************
 *       Copyright (C) 2011 Jeroen Demeyer <jdemeyer@cage.ugent.be>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/

#include <mpir.h>
#include "memory.h"

/* mpir memory functions */
void* sage_mpir_malloc(size_t size)
{
    return sage_malloc(size);
}

void* sage_mpir_realloc(void *ptr, size_t old_size, size_t new_size)
{
    return sage_realloc(ptr, new_size);
}

void sage_mpir_free(void *ptr, size_t size)
{
    sage_free(ptr);
}

void init_memory_functions()
{
#if 0
    void* (*mpir_malloc)(size_t);
    void* (*mpir_realloc)(void*, size_t, size_t);
    void (*mpir_free)(void*, size_t);
    mp_get_memory_functions(&mpir_malloc, &mpir_realloc, &mpir_free);
    printf("MPIR memory functions BEFORE: %p %p %p\n", mpir_malloc, mpir_realloc, mpir_free);
#endif
    mp_set_memory_functions(sage_mpir_malloc, sage_mpir_realloc, sage_mpir_free);
#if 0
    mp_get_memory_functions(&mpir_malloc, &mpir_realloc, &mpir_free);
    printf("MPIR memory functions AFTER:  %p %p %p\n", mpir_malloc, mpir_realloc, mpir_free);
#endif
}
