from cpython.object cimport *

cdef class lazy_list_abstract(object):
    cdef list cache                  # the cache
    cdef lazy_list_abstract master   # a reference if self is a slice
    cdef Py_ssize_t start, stop, step

    cpdef get(self, Py_ssize_t i)
    cpdef int _fit(self, Py_ssize_t n) except -1
    cdef int update_cache_up_to(self, Py_ssize_t i) except -1

cdef class lazy_list_from_iterator(lazy_list_abstract):
    cdef object iterator

cdef class lazy_list_from_function(lazy_list_abstract):
    cdef object callable

cdef class lazy_list_from_update_function(lazy_list_abstract):
    cdef object update_function
