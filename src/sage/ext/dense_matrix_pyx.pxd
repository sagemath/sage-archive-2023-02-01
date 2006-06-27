include "cdefs.pxi"

cdef class Matrix_rational:
    """
    Matrix over the rational numbers.
    """
    cdef mpq_t **matrix
    cdef mpq_t tmp
    cdef int _nrows, _ncols
    cdef object __pivots
    cdef int initialized

    cdef set_matrix(Matrix_rational self, mpq_t **m)
    cdef scale_row(self, int row, mpq_t multiple, int start_col)
    cdef add_multiple_of_row(self, int row_from, mpq_t multiple,
                            int row_to, int start_col)
    cdef add_multiple_of_column(self, int col_from, mpq_t multiple,
                               int col_to, int start_row)
    cdef swap_rows(self, int row1, int row2)
    cdef swap_columns(self, int col1, int col2)
    cdef int mpz_denom(self, mpz_t d) except -1
    cdef int mpz_height(self, mpz_t height) except -1
