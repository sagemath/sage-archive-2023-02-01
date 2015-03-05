"""
Sage low-level memory allocation functions
"""

#*****************************************************************************
#       Copyright (C) 2014 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


cdef extern from "memory.h":
    size_t mul_overflowcheck(size_t a, size_t b)
    void  sage_free(void* ptr) nogil
    void* sage_malloc(size_t) nogil
    void* sage_realloc(void* ptr, size_t n) nogil
    void* sage_calloc(size_t nmemb, size_t size) nogil

cdef inline void* check_allocarray(size_t nmemb, size_t size) except? NULL:
    """
    Allocate memory for ``nmemb`` elements of size ``size``.
    """
    if nmemb == 0:
        return NULL
    cdef size_t n = mul_overflowcheck(nmemb, size)
    cdef void* ret = sage_malloc(n)
    if ret == NULL:
        raise MemoryError("failed to allocate %s * %s bytes" % (nmemb, size))
    return ret

cdef inline void* check_reallocarray(void* ptr, size_t nmemb, size_t size) except? NULL:
    """
    Re-allocate memory at ``ptr`` to hold ``nmemb`` elements of size
    ``size``. If ``ptr`` equals ``NULL``, this behaves as
    ``check_allocarray``.

    When ``nmemb`` equals 0, then free the memory at ``ptr``.
    """
    if nmemb == 0:
        sage_free(ptr)
        return NULL
    cdef size_t n = mul_overflowcheck(nmemb, size)
    cdef void* ret = sage_realloc(ptr, n)
    if ret == NULL:
        raise MemoryError("failed to allocate %s * %s bytes" % (nmemb, size))
    return ret

cdef inline void* check_malloc(size_t n) except? NULL:
    """
    Allocate ``n`` bytes of memory.
    """
    if n == 0:
        return NULL
    cdef void* ret = sage_malloc(n)
    if ret == NULL:
        raise MemoryError("failed to allocate %s bytes" % n)
    return ret

cdef inline void* check_realloc(void* ptr, size_t n) except? NULL:
    """
    Re-allocate memory at ``ptr`` to hold ``n`` bytes.
    If ``ptr`` equals ``NULL``, this behaves as ``check_malloc``.
    """
    if n == 0:
        sage_free(ptr)
        return NULL
    cdef void* ret = sage_realloc(ptr, n)
    if ret == NULL:
        raise MemoryError("failed to allocate %s bytes" % n)
    return ret

cdef inline void* check_calloc(size_t nmemb, size_t size) except? NULL:
    """
    Allocate memory for ``nmemb`` elements of size ``size``. The
    resulting memory is zeroed.
    """
    if nmemb == 0:
        return NULL
    cdef void* ret = sage_calloc(nmemb, size)
    if ret == NULL:
        raise MemoryError("failed to allocate %s * %s bytes" % (nmemb, size))
    return ret
