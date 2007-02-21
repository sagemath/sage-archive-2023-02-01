include '../ext/cdefs.pxi'
include '../modules/vector_rational_sparse_h.pxi'

cimport matrix_sparse

cdef class Matrix_rational_sparse(matrix_sparse.Matrix_sparse):
    cdef mpq_vector* _matrix
    cdef int _initialized

