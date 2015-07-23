cdef class _LazyString(object):
    cdef func
    cdef args
    cdef kwargs
    cpdef update_lazy_string(self, tuple args, dict kwds)
