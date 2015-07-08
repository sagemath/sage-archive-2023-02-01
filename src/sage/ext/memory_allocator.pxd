from sage.ext.memory cimport malloc, free, realloc, check_calloc, check_allocarray, check_realloc, check_reallocarray

cdef class MemoryAllocator:
    cdef size_t n
    cdef size_t max_size
    cdef void ** pointers
    cdef void * malloc(self, size_t size) except NULL
    cdef void * calloc(self, size_t nmemb, size_t size) except NULL
    cdef enlarge_if_needed(self) except -1
