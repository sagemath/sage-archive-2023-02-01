# cython: profile=False

include "../ext/cdefs.pxi"
include "../ext/gmp.pxi"
include "../ext/stdsage.pxi"

cdef void* pymem_realloc(void *ptr, size_t old_size, size_t new_size):
    return sage_realloc(ptr, new_size)

cdef void pymem_free(void *ptr, size_t size):
    sage_free(ptr)

cdef void* pymem_malloc(size_t size):
    return sage_malloc(size)

cdef extern from "gmp.h":
    void mp_set_memory_functions (void *(*alloc_func_ptr) (size_t),  \
                                  void *(*realloc_func_ptr) (void *, size_t, size_t),    \
                                  void (*free_func_ptr) (void *, size_t))


def pmem_malloc():
    """
    Use our own memory manager for for GMP memory management.


    WARNING: Call this before the integer initialization or evil
    things might happen.
    """
    mp_set_memory_functions(sage_malloc, pymem_realloc, pymem_free)
    #mp_set_memory_functions(PyMem_Malloc, pymem_realloc, pymem_free)
    #mp_set_memory_functions(pymem_malloc, pymem_realloc, pymem_free)

pmem_malloc()

