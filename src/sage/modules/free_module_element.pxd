cimport sage.structure.element
import  sage.structure.element

cdef class FreeModuleElement(sage.structure.element.ModuleElement):
    pass

cdef class FreeModuleElement_generic_dense(FreeModuleElement):
    cdef object _entries

cdef class FreeModuleElement_generic_sparse(FreeModuleElement):
    cdef object _entries
