from sage.rings.polynomial.skew_polynomial_finite_order cimport SkewPolynomial_finite_order_dense
from sage.matrix.matrix_dense cimport Matrix_dense

cdef class SkewPolynomial_finite_field_dense (SkewPolynomial_finite_order_dense):
    cdef _norm_factor
    cdef dict _rdivisors
    cdef dict _types
    cdef _factorization

    # Finding divisors
    cdef SkewPolynomial_finite_field_dense _rdivisor_c(P, N)

    # Finding factorizations
    cdef _factor_c(self)
    cdef _factor_uniform_c(self)
