cimport matrix_sparse

cdef class Matrix_generic_sparse(matrix_sparse.Matrix_sparse):
    cdef object _entries
    cdef object _zero
    cdef object _base_ring
