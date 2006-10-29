cimport matrix_sparse

cdef struct c_vector_modint:
    int *entries
    int *positions
    int p
    int degree
    int num_nonzero


cdef class Matrix_modn_sparse(matrix_pid_sparse.Matrix_pid_sparse):
    cdef c_vector_modint* rows
    cdef public int p
