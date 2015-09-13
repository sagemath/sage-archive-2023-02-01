#*****************************************************************************
#       Copyright (C) 2015 Simon King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef class FieldConverter_class:
    cdef object field  # that's a function converting an int to a field element
    cdef object int_to_field(self, int x)
    cdef int field_to_int(self, x)

from sage.matrix.matrix_dense cimport Matrix_dense
from sage.structure.element cimport Matrix
from sage.libs.meataxe cimport *

cdef class Matrix_gfpn_dense(Matrix_dense):
    cdef Matrix_t *Data
    cdef FieldConverter_class _converter
    #cpdef Matrix_gfpn_dense normalized(Matrix_gfpn_dense self)
    #cpdef Matrix_gfpn_dense semi_echelon(Matrix_gfpn_dense self)
    #cpdef int nullity(Matrix_gfpn_dense self)
    #cpdef tuple lead(self)
    cdef set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value)
    cdef inline int get_unsafe_int(self, Py_ssize_t i, Py_ssize_t j)
    cdef Matrix _matrix_times_matrix_(self, Matrix right)
    cpdef Matrix_gfpn_dense _multiply_classical(Matrix_gfpn_dense self, Matrix_gfpn_dense right)
    cpdef Matrix_gfpn_dense _multiply_strassen(Matrix_gfpn_dense self, Matrix_gfpn_dense right, cutoff=*)
