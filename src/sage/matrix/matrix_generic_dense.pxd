cimport matrix_dense

cdef class Matrix_generic_dense(matrix_dense.Matrix_dense):
    cdef object _entries

cimport matrix_window

cdef class MatrixWindow(matrix_window.MatrixWindow):
    cdef Matrix_generic_dense _matrix
    cdef int _row
    cdef int _col
    cdef int _nrows
    cdef int _ncols
