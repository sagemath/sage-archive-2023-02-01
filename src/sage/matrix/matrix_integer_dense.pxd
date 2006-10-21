include "../ext/cdefs.pxi"

cimport matrix_pid

cdef class Matrix_integer_dense(matrix_pid.Matrix_pid):
    cdef mpz_t *_entries
    cdef mpz_t **_matrix
    cdef mpz_t tmp
    cdef object __pivots
    cdef int initialized
    cdef int mpz_height(self, mpz_t height) except -1





