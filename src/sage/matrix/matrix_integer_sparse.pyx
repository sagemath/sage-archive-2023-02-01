"""nodoctest
Sparse matrices over the integers.
"""

cimport matrix_pid_sparse
import  matrix_pid_sparse

cdef class Matrix_integer_sparse(matrix_pid_sparse.Matrix_pid_sparse):

    def __init__(self, parent, coerce_entries=None, copy=True):
        matrix_sparse.Matrix_pid_sparse.__init__(self, parent)

    def echelon_form(self, *args, **kwds):
        """
        This function works in the same was as echelon form for
        \code{self.dense_matrix()}.  Please see the documentation
        for that function.
        """
        return self.dense_matrix().echelon_form(*args, **kwds)
