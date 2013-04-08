include '../ext/cdefs.pxi'

cimport free_module_element
cimport numpy as cnumpy

cdef class Vector_double_dense(free_module_element.FreeModuleElement):
    cdef object _numpy_dtype
    # cdef cnumpy.NPY_TYPES _numpy_dtypeint
    cdef int _numpy_dtypeint
    cdef object _python_dtype
    cdef object _sage_dtype
    cdef object _sage_vector_dtype
    cdef cnumpy.ndarray _vector_numpy
    cdef set_unsafe(self, Py_ssize_t i, object value)
    cdef get_unsafe(self, Py_ssize_t i)
    cdef Vector_double_dense _new(self, cnumpy.ndarray vector_numpy)
    cdef _replace_self_with_numpy(self, cnumpy.ndarray numpy_array)
