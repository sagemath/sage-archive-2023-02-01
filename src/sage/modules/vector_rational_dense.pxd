cimport free_module_element
import  free_module_element

include 'sage/ext/cdefs.pxi'

cdef class Vector_rational_dense(free_module_element.FreeModuleElement):
        cdef mpq_t* _entries
        cdef _new_c(self)
        cdef _init(self, Py_ssize_t degree, parent)

