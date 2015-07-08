include 'sage/ext/interrupt.pxi'

from sage.ext.memory cimport malloc, free, realloc

cdef class MemoryAllocator:
    cdef int n
    cdef int max_size
    cdef void ** pointers
    cdef void * malloc(self, size_t size) except NULL
    cdef enlarge_if_needed(self)
