include "../ext/cdefs.pxi"

ctypedef size_t mod_int

cdef class Linbox:
    cdef int modn_dense_echelonize(self, mod_int modulus, mod_int** matrix,
                                   size_t nrows, size_t ncols)

    cdef modn_dense_poly(self, unsigned long modulus, size_t n,
                                mod_int **matrix, minpoly)

