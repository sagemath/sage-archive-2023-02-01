include "../ext/cdefs.pxi"

cimport matrix_pyx

cdef class Matrix_dense(matrix_pyx.Matrix):
    cdef int _nrows, _ncols
    cdef object _entries
    cdef int* _row_indices


#cdef object multiply_items(object v, int i, object w, int j)

