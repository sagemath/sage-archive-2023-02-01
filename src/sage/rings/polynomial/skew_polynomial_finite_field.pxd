from skew_polynomial_element cimport SkewPolynomial_generic_dense
from sage.matrix.matrix_dense cimport Matrix_dense

cdef class SkewPolynomial_finite_field_dense (SkewPolynomial_generic_dense):
    cdef Matrix_dense _matmul_c(self)
