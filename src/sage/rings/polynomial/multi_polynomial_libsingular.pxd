include "../../libs/singular/singular-cdefs.pxi"

#cimport sage.rings.polynomial.multi_polynomial
cimport multi_polynomial
from sage.rings.polynomial.multi_polynomial_ring_generic cimport MPolynomialRing_generic
from sage.structure.parent cimport Parent

cdef class MPolynomial_libsingular(sage.rings.polynomial.multi_polynomial.MPolynomial):
    cdef poly *_poly
    cpdef _repr_short_(self)
    cpdef is_constant(self)
    cpdef _homogenize(self, int var)

cdef class MPolynomialRing_libsingular(MPolynomialRing_generic):
    cdef object __singular
    cdef object __macaulay2
    cdef object __m2_set_ring_cache
    cdef object __minpoly
    cdef ring *_ring
    cdef int _cmp_c_impl(left, Parent right) except -2

    #cpdef MPolynomial_libsingular _element_constructor_(self, element)

# new polynomials

cdef MPolynomial_libsingular new_MP(MPolynomialRing_libsingular parent, poly *p)
