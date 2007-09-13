from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.libs.ntl.ntl cimport ntl_ZZX

cdef class Polynomial_integer_dense_ntl(Polynomial):
    cdef ntl_ZZX __poly
