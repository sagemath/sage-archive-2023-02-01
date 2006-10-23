include "../ext/cdefs.pxi"
include "../ext/interrupt.pxi"


cimport matrix_pid

cdef class Matrix_integer_dense(matrix_pid.Matrix_pid):
    cdef mpz_t *_entries
    cdef mpz_t **_matrix
    cdef mpz_t tmp
    cdef object _pivots
    cdef int mpz_height(self, mpz_t height) except -1

    cdef void _zero_out_matrix(self)





