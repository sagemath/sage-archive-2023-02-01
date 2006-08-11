include "../ext/cdefs.pxi"

cimport matrix_pid

cdef class Matrix_integer_dense(matrix_pid.Matrix_pid):
    cdef mpz_t **matrix
    cdef mpz_t tmp
    cdef int _nrows, _ncols
    cdef object __pivots
    cdef int initialized
    cdef int mpz_height(self, mpz_t height) except -1
    cdef int allocate(self) except -1
    cdef set_matrix(Matrix_integer_dense self, mpz_t **m)





