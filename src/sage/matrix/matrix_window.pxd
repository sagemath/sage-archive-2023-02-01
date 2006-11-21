from matrix cimport Matrix

cdef class MatrixWindow:
    cdef Py_ssize_t _row, _col, _nrows, _ncols
    cdef Matrix _matrix
