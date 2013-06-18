include "sage/ext/cdefs.pxi"

cimport matrix_dense

cdef class Matrix_integer_2x2(matrix_dense.Matrix_dense):
    cdef mpz_t a
    cdef mpz_t b
    cdef mpz_t c
    cdef mpz_t d
    cdef mpz_t *_entries

    cdef Matrix_integer_2x2 _new_c(self)


