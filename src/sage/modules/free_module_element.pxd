from sage.structure.element cimport ModuleElement
from sage.matrix.matrix cimport Matrix

cdef class FreeModuleElement(ModuleElement):
    pass
    #cdef FreeModuleElement _matrix_multiply(self, Matrix A)

cdef class FreeModuleElement_generic_dense(FreeModuleElement):
    cdef _new_c(self, object v)
    cdef object _entries

cdef class FreeModuleElement_generic_sparse(FreeModuleElement):
    cdef _new_c(self, object v)
    cdef object _entries
