include "../ext/cdefs.pxi"

cimport matrix_dense

cdef class Matrix_integer_dense(matrix_dense.Matrix_dense):
    cdef char _initialized
    cdef mpz_t *_entries
    cdef mpz_t **_matrix
    cdef object _pivots
    cdef int mpz_height(self, mpz_t height) except -1

    cdef void _zero_out_matrix(self)
    cdef _new_unitialized_matrix(self, Py_ssize_t nrows, Py_ssize_t ncols)
    cdef _pickle_version0(self)
    cdef _unpickle_version0(self, data)




