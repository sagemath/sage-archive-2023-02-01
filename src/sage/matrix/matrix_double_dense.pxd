cimport matrix_dense
include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'
include '../ext/python.pxi'
include '../ext/numpy.pxd'
import matrix
cimport matrix

from sage.rings.ring cimport Field

import numpy

cdef class Matrix_double_dense(matrix_dense.Matrix_dense):
    cdef object _numpy_dtype
    cdef NPY_TYPES _numpy_dtypeint
    cdef object _python_dtype
    cdef object _sage_dtype
    cdef object _sage_vector_dtype
    cdef ndarray _L, _U, _P
    cdef int _LU_valid
    cdef Matrix_double_dense _new(self, int nrows=*, int ncols=*)
    cdef _c_compute_LU(self)
    cdef ndarray _matrix_numpy
