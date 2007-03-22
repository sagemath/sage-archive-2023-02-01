include "../../ext/cdefs.pxi"

cdef extern from "../libs/m4ri/packedmatrix.h":
    ctypedef struct packedmatrix

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


## cdef class Linbox_mod2_dense:
##     cdef packedmatrix *matrix

##     cdef set(self, packedmatrix *matrix)
##     cdef int echelonize(self)
##     cdef matrix_matrix_multiply(self, packedmatrix *ans, packedmatrix *B)
##     cdef unsigned long rank(self) except -1

cdef class Linbox_integer_dense:
    cdef mpz_t** matrix
    cdef size_t nrows, ncols

    cdef set(self, mpz_t** matrix, size_t nrows, size_t ncols)
    cdef matrix_matrix_multiply(self,
                                mpz_t **ans,
                                mpz_t **B,
                                size_t B_nr, size_t B_nc)
    cdef unsigned long rank(self) except -1

