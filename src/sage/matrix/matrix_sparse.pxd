cimport matrix_generic

cdef class Matrix_sparse(matrix_generic.Matrix):
    cdef object __entries # TODO: move to Matrix?


cdef class Matrix_domain_sparse(Matrix_sparse):
    pass

cdef class Matrix_pid_sparse(Matrix_domain_sparse):
    pass

cdef class Matrix_field_sparse(Matrix_pid_sparse):
    pass
