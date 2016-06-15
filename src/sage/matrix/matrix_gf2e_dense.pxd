from sage.libs.m4rie cimport mzed_t

cimport matrix_dense

cdef class Matrix_gf2e_dense(matrix_dense.Matrix_dense):
    cdef mzed_t *_entries
    cdef object _one
    cdef object _zero

    cpdef Matrix_gf2e_dense _multiply_newton_john(Matrix_gf2e_dense self, Matrix_gf2e_dense right)
    cpdef Matrix_gf2e_dense _multiply_karatsuba(Matrix_gf2e_dense self, Matrix_gf2e_dense right)
    cpdef Matrix_gf2e_dense _multiply_strassen(Matrix_gf2e_dense self, Matrix_gf2e_dense right, cutoff=*)
