from sage.structure.element cimport Vector

cdef class FreeModuleElement(Vector):
    cdef bint _is_mutable
    cdef int _cmp_same_ambient_c(left, FreeModuleElement right) except -2
    cpdef int hamming_weight(self)

cdef class FreeModuleElement_generic_dense(FreeModuleElement):
    # data
    cdef object _entries

    # cdef'd methods
    cdef _new_c(self, object v)


cdef class FreeModuleElement_generic_sparse(FreeModuleElement):
    # data
    cdef object _entries

    # cdef'd methods
    cdef _new_c(self, object v)

