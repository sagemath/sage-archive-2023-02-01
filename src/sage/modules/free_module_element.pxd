from sage.structure.element cimport Vector

cdef class FreeModuleElement(Vector):
    cdef bint _is_mutable
    cpdef int hamming_weight(self)

cdef class FreeModuleElement_generic_dense(FreeModuleElement):
    # data
    cdef list _entries

    # cdef'd methods
    cdef _new_c(self, object v)


cdef class FreeModuleElement_generic_sparse(FreeModuleElement):
    # data
    cdef dict _entries

    # cdef'd methods
    cdef _new_c(self, object v)

