include "../ext/cdefs.pxi"
include '../ext/interrupt.pxi'

cimport matrix_sparse


cdef struct c_vector_modint:
    int *entries
    int *positions
    int p
    int degree
    int num_nonzero


cdef struct mpq_vector:
    mpq_t *entries      # array of nonzero entries
    int   *positions    # positions of those nonzero entries, starting at 0
    int    degree       # the degree of this sparse vector
    int    num_nonzero  # the number of nonzero entries of this vector.


cdef class Matrix_rational_sparse(matrix_sparse.Matrix_sparse):
    cdef mpq_vector* rows
    cdef public int nr, nc
    cdef public is_init

