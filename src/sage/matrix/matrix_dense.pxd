cimport matrix

cdef class Matrix_dense(matrix.Matrix):
    cdef set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value)

