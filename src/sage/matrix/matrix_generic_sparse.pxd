cimport matrix

cdef class Matrix_generic_sparse(matrix.Matrix):
    cdef object _entries
    cdef object _zero
    cdef object _base_ring
