include "../ext/cdefs.pxi"

cimport matrix

cdef class Matrix_generic_dense(matrix.Matrix):
    cdef object _entries

cdef class MatrixWindow(matrix.MatrixWindow):
    cdef Matrix_generic_dense _matrix
    cdef int _row
    cdef int _col
    cdef int _nrows
    cdef int _ncols
