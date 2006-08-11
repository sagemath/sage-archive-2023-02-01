cimport sage.matrix.matrix_dense

cdef class MatrixWindow:

    cdef sage.matrix.matrix_dense.Matrix_dense _matrix

    cdef int _row
    cdef int _col
    cdef int _nrows
    cdef int _ncols
