include "../../ext/cdefs.pxi"
include "../../libs/flint/fmpz.pxi"
include "../../libs/flint/fmpz_poly.pxi"
include "../../libs/flint/ntl_interface.pxi"

from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.rings.integer cimport Integer

cdef class Polynomial_integer_dense_flint(Polynomial):
    cdef fmpz_poly_t __poly

    cdef Polynomial_integer_dense_flint _new(self)
    cpdef bint is_zero(self)
    cpdef _unsafe_mutate(self, long n, value)
    cpdef Integer content(self)
