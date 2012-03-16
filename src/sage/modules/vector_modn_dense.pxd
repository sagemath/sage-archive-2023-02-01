cdef extern from "stdint.h":
    ctypedef unsigned long uint_fast64_t

cdef extern from "../ext/multi_modular.h":
    ctypedef uint_fast64_t mod_int

cimport free_module_element
import  free_module_element

cdef class Vector_modn_dense(free_module_element.FreeModuleElement):
    cdef mod_int* _entries
    cdef mod_int _p
    cdef object _base_ring

    cdef _new_c(self)
    cdef _init(self, Py_ssize_t degree, parent, mod_int p)

