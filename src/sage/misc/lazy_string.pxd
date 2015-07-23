cdef class _LazyString(object):
    cdef object _func
    cdef tuple _args
    cdef dict _kwargs
    cpdef update_lazy_string(self, tuple args, dict kwds)
