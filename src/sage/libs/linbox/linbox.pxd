include "../../ext/cdefs.pxi"

ctypedef size_t mod_int

cdef class Linbox_modn_dense:
    cdef mod_int** matrix
    cdef mod_int n
    cdef size_t nrows, ncols

    cdef set(self, mod_int n, mod_int** matrix, size_t nrows, size_t ncols)
    cdef int echelonize(self)
    cdef matrix_matrix_multiply(self,
                                mod_int **ans,
                                mod_int **B,
                                size_t B_nr, size_t B_nc)
    cdef unsigned long rank(self) except -1



cdef class Linbox_integer_dense:
    cdef mpz_t** matrix
    cdef size_t nrows, ncols

    cdef set(self, mpz_t** matrix, size_t nrows, size_t ncols)
    cdef matrix_matrix_multiply(self,
                                mpz_t **ans,
                                mpz_t **B,
                                size_t B_nr, size_t B_nc)
    cdef unsigned long rank(self) except -1
