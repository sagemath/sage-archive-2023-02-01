cimport sage.structure.element
cimport sage.matrix.matrix

cdef class FreeModuleElement(sage.structure.element.ModuleElement):
    cdef FreeModuleElement _matrix_multiply(self, sage.matrix.matrix.Matrix A)
    cdef FreeModuleElement _scalar_multiply(self, s)

cdef class FreeModuleElement_generic_dense(FreeModuleElement):
    cdef object _entries

cdef class FreeModuleElement_generic_sparse(FreeModuleElement):
    cdef object _entries
