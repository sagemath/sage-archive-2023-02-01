cimport matrix_generic

ctypedef unsigned int uint


cdef class Matrix_modn_dense(matrix_generic.Matrix):
    cdef uint **matrix
    cdef uint _nrows, _ncols, p
    cdef uint gather
    cdef object __pivots
    cdef set_matrix(Matrix_modn_dense self, uint **m)
    cdef uint **get_matrix(Matrix_modn_dense self)
    cdef uint entry(self, uint i, uint j)
    cdef scale_row(self, uint row, uint multiple, uint start_col)
    cdef add_multiple_of_row(self, uint row_from, uint multiple, \
                            uint row_to, uint start_col)
    cdef add_multiple_of_column(self, uint col_from, uint multiple, \
                               uint col_to, uint start_row)
    cdef swap_rows(self, uint row1, uint row2)
    cdef swap_columns(self, uint col1, uint col2)







