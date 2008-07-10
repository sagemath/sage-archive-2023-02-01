cimport matrix_double_dense
cimport matrix_dense
include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'
include '../ext/python.pxi'
include '../ext/numpy.pxd'
import matrix
cimport matrix

import numpy

cdef class Matrix_complex_double_dense(matrix_double_dense.Matrix_double_dense):
    pass

