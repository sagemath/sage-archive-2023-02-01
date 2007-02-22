cimport matrix_sparse

cdef struct c_vector_modint:
    int *entries
    int p
    Py_ssize_t *positions
    Py_ssize_t degree
    Py_ssize_t num_nonzero

cdef class Matrix_modn_sparse(matrix_sparse.Matrix_sparse):
    cdef c_vector_modint* rows
    cdef public int p
    cdef swap_rows_c(self, Py_ssize_t n1, Py_ssize_t n2)
