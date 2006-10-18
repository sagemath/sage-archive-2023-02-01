"""
Sparse matrices over the integers.
"""

cimport matrix_integer
import matrix_integer

cdef class Matrix_integer_sparse(matrix_integer.Matrix_integer):

    def __init__(self, parent, coerce_entries=None, copy=True):
        matrix_integer.Matrix_integer.__init__(self, parent)
