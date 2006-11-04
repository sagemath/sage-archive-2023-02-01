cimport matrix_dense

ctypedef unsigned int uint

cdef class Matrix_modn_dense(matrix_dense.Matrix_dense):
    cdef uint **matrix
    cdef uint *_entries
    cdef uint p
    cdef uint gather

    #cdef set_matrix(Matrix_modn_dense self, uint **m)
    #cdef uint **get_matrix(Matrix_modn_dense self)
    #cdef uint entry(self, uint i, uint j)

    cdef rescale_row_c(self, Py_ssize_t row, uint multiple, Py_ssize_t start_col)
    cdef rescale_col_c(self, Py_ssize_t col, uint multiple, Py_ssize_t start_row)
    cdef add_multiple_of_row_c(self, uint row_to, uint row_from, uint multiple, \
                            uint start_col)
    cdef add_multiple_of_column_c(self, uint col_to, uint col_from, uint multiple, \
                                  uint start_row)
    cdef swap_rows_c(self, uint row1, uint row2)
    cdef swap_columns_c(self, uint col1, uint col2)



#cdef class MatrixWindow:
#    cdef Matrix_modn_dense _matrix
#    cdef int _row, _col, _nrows, _ncols
