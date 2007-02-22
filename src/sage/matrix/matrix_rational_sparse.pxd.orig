include "../ext/cdefs.pxi"

cimport matrix_sparse

cdef class Matrix_rational_sparse(matrix_sparse.Matrix_sparse):
    cdef mpq_vector* rows
    cdef public int nr, nc
    cdef public is_init

