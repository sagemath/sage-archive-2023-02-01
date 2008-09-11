from matrix_dense cimport Matrix_dense
from matrix_rational_dense cimport Matrix_rational_dense

cdef class Matrix_cyclo_dense(Matrix_dense):
    # Matrix over ZZ that stores elements
    cdef Matrix_rational_dense _matrix

    # Degree of base cyclotomic field
    cdef int _degree
    cdef int _n

