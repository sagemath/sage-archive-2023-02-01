from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense

include "sage/libs/singular/singular-cdefs.pxi"

cdef class Matrix_mpolynomial_dense(Matrix_generic_dense):
    cdef ideal *_to_libsingular(self, bint module)
    cdef void _from_libsingular(self, ideal *m)
