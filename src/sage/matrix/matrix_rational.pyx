"""
Generic matrices over the rational field.
"""

cimport matrix_field
import matrix_field

cdef class Matrix_rational(matrix_field.Matrix_field):
    def __init__(self, parent):
        matrix_field.Matrix_field.__init__(self, parent)
