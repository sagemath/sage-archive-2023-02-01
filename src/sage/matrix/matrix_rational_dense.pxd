include "../ext/cdefs.pxi"

cimport matrix_field_dense

cdef class Matrix_rational_dense(matrix_field_dense.Matrix_field_dense):

    cdef mpq_t tmp
    cdef mpq_t *_entries
    cdef mpq_t ** _matrix
    cdef object __pivots

    cdef scale_row(self, int row, mpq_t multiple, int start_col)
    cdef add_multiple_of_row(self, int row_from, mpq_t multiple,
                            int row_to, int start_col)
    cdef add_multiple_of_column(self, int col_from, mpq_t multiple,
                               int col_to, int start_row)
    cdef swap_rows(self, int row1, int row2)
    cdef swap_columns(self, int col1, int col2)
    cdef int mpz_denom(self, mpz_t d) except -1
    cdef int mpz_height(self, mpz_t height) except -1
    cdef int _rescale(self, mpq_t a) except -1


cdef class MatrixWindow:
    cdef Matrix_rational_dense _matrix
    cdef int _row, _col, _nrows, _ncols
