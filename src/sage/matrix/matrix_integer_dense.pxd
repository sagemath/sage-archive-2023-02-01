include "../ext/cdefs.pxi"

cimport matrix_integer

cdef class Matrix_integer_dense(matrix_integer.Matrix_integer):
    cdef mpz_t *_entries
    cdef mpz_t **_matrix
    cdef mpz_t tmp
    cdef int _nrows, _ncols
    cdef object __pivots
    cdef int initialized
    cdef int mpz_height(self, mpz_t height) except -1





