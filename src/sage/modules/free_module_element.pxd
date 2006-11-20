from sage.structure.element cimport Vector

cdef class FreeModuleElement(Vector):
    pass

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

