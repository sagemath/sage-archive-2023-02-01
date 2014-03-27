from sage.ext.mod_int cimport *

cimport free_module_element
import  free_module_element

cdef class Vector_modn_dense(free_module_element.FreeModuleElement):
    cdef mod_int* _entries
    cdef mod_int _p
    cdef object _base_ring

    cdef _new_c(self)
    cdef _init(self, Py_ssize_t degree, parent, mod_int p)

