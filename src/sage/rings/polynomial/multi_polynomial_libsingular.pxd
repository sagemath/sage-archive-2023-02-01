include "../../libs/singular/singular-cdefs.pxi"

cimport sage.rings.polynomial.multi_polynomial
from sage.rings.polynomial.multi_polynomial_ring_generic cimport MPolynomialRing_generic
from sage.structure.parent cimport Parent

cdef class MPolynomial_libsingular(sage.rings.polynomial.multi_polynomial.MPolynomial):
    cdef object __singular
    cdef poly *_poly
    cdef _repr_short_c(self)
    cdef _singular_init_c(self,singular, have_ring)
    cdef int is_constant_c(self)

cdef class MPolynomialRing_libsingular(MPolynomialRing_generic):
    cdef object __weakref__
    cdef object __singular
    cdef object __macaulay2
    cdef object __m2_set_ring_cache
    cdef ring *_ring
    cdef MPolynomial_libsingular _zero
    cdef int _cmp_c_impl(left, Parent right) except -2
