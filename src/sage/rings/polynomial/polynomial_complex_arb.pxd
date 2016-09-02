from sage.libs.arb.acb_poly cimport *
from sage.rings.polynomial.polynomial_element cimport Polynomial_generic_dense

cdef class Polynomial_complex_arb(Polynomial_generic_dense):
    pass
