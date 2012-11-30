include "../ext/python_object.pxi"

cdef class lazy_list(object):
    cdef object iterator
    cdef list cache
    cdef Py_ssize_t start, stop, step

    cdef int update_cache_up_to(self, Py_ssize_t i) except -1

cdef class lazy_list_iterator(object):
    cdef lazy_list l
    cdef Py_ssize_t pos, step

cdef class stopped_lazy_list_iterator(object):
    cdef lazy_list l
    cdef Py_ssize_t pos, step, stop


