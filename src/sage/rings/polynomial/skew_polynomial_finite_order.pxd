from sage.rings.polynomial.skew_polynomial_element cimport SkewPolynomial_generic_dense
from sage.matrix.matrix_dense cimport Matrix_dense

cdef class SkewPolynomial_finite_order_dense (SkewPolynomial_generic_dense):
    cdef _norm
    cdef _optbound

    cdef Matrix_dense _matphir_c(self)
    cdef Matrix_dense _matmul_c(self)

