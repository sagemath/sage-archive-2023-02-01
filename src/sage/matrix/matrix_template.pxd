cimport matrix

cdef class Matrix_generic_dense(matrix.Matrix):
    cdef object _entries
    cdef int*   _row_indices





