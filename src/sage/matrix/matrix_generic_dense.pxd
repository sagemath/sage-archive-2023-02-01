cimport matrix_dense

cdef class Matrix_generic_dense(matrix_dense.Matrix_dense):
    cdef object _entries

cdef class MatrixWindow(matrix.MatrixWindow):
    cdef Matrix_generic_dense _matrix
    cdef int _row
    cdef int _col
    cdef int _nrows
    cdef int _ncols
