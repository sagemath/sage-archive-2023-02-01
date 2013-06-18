cimport free_module_element
import  free_module_element

from sage.libs.m4ri cimport mzd_t

cdef class Vector_mod2_dense(free_module_element.FreeModuleElement):
    cdef mzd_t *_entries
    cdef object _base_ring

    cdef _new_c(self)
    cdef _init(self, Py_ssize_t degree, parent)
