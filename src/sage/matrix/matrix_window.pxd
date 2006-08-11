cimport matrix_dense

cdef class MatrixWindow:

    cdef matrix_dense.Matrix_dense _matrix

    cdef int _row
    cdef int _col
    cdef int _nrows
    cdef int _ncols
