"""
Sparse matrices over the integers.
"""

cimport matrix_pid_sparse
import  matrix_pid_sparse

cdef class Matrix_integer_sparse(matrix_pid_sparse.Matrix_pid_sparse):

    def __init__(self, parent, coerce_entries=None, copy=True):
        matrix_sparse.Matrix_pid_sparse.__init__(self, parent)
