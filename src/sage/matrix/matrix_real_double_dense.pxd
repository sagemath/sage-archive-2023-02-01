cimport matrix_double_dense
cimport matrix_dense
include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'
include '../ext/python.pxi'
include '../ext/numpy.pxd'
import matrix
cimport matrix

import numpy

cdef class Matrix_real_double_dense(matrix_double_dense.Matrix_double_dense):
    cdef set_unsafe_double(self, Py_ssize_t i, Py_ssize_t j, double value)
    cdef double get_unsafe_double(self, Py_ssize_t i, Py_ssize_t j)


