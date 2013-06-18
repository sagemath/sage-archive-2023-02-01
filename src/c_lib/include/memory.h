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

#ifndef C_LIB_INCLUDE_MEMORY_H
#define C_LIB_INCLUDE_MEMORY_H

#include "interrupt.h"

#ifdef __cplusplus
extern "C" {
#endif

static inline void* sage_malloc(size_t n)
{
    sig_block();
    void* ret = malloc(n);
    sig_unblock();
    return ret;
}
static inline void* sage_calloc(size_t nmemb, size_t size)
{
    sig_block();
    void* ret = calloc(nmemb, size);
    sig_unblock();
    return ret;
}
static inline void sage_free(void* ptr)
{
    sig_block();
    free(ptr);
    sig_unblock();
}
static inline void* sage_realloc(void *ptr, size_t size)
{
    sig_block();
    void* ret = realloc(ptr, size);
    sig_unblock();
    return ret;
}

void init_memory_functions();


#ifdef __cplusplus
}  /* extern "C" */
#endif
#endif /* C_LIB_INCLUDE_MEMORY_H */
