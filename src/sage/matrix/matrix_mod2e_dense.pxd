# choose: dense or sparse

from sage.rings.finite_rings.element_givaro cimport GivaroGfq, Cache_givaro

from sage.libs.m4rie cimport mzed_t

cimport matrix_dense

cdef class Matrix_mod2e_dense(matrix_dense.Matrix_dense):
    cdef mzed_t *_entries
    cdef Cache_givaro cc
    cdef object _one
    cdef object _zero

    cpdef Matrix_mod2e_dense _multiply_newton_john(Matrix_mod2e_dense self, Matrix_mod2e_dense right)
    cpdef Matrix_mod2e_dense _multiply_karatsuba(Matrix_mod2e_dense self, Matrix_mod2e_dense right)
    cpdef Matrix_mod2e_dense _multiply_strassen(Matrix_mod2e_dense self, Matrix_mod2e_dense right, cutoff=*)
