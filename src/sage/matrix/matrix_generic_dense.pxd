cimport matrix_dense

cdef class Matrix_generic_dense(matrix_dense.Matrix_dense):
    cdef object _entries

