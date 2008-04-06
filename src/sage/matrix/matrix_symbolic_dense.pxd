cimport matrix_dense

cdef class Matrix_symbolic_dense(matrix_dense.Matrix_dense):
    cdef object _maxima
    cdef set_from_list(self, entries)
    cdef object __variables
    cdef object __number_of_args

    cdef public object _simp
