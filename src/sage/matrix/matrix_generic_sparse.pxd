cimport matrix_generic

cdef class Matrix_sparse(matrix_generic.Matrix):
    cdef object _entries
    cdef object _base_ring
    cdef object _zero
