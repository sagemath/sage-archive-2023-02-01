include "singular-cdefs.pxi"

from sage.rings.rational cimport Rational
from sage.structure.element cimport Element
from sage.rings.integer cimport Integer
from sage.rings.finite_field_givaro cimport FiniteField_givaro,FiniteField_givaroElement
from sage.rings.finite_field_ntl_gf2e cimport FiniteField_ntl_gf2e,FiniteField_ntl_gf2eElement
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular

from sage.rings.number_field.number_field_base cimport NumberField

cdef class Conversion:
    cdef extern Rational si2sa_QQ(self, number (*),ring (*))
    cdef extern Integer si2sa_ZZ(self, number (*),ring (*))
    cdef extern FiniteField_givaroElement si2sa_GFqGivaro(self, number *n, ring *_ring, FiniteField_givaro base)
    cdef extern FiniteField_ntl_gf2eElement si2sa_GFqNTLGF2E(self, number *n, ring *_ring, FiniteField_ntl_gf2e base)
    cdef extern object si2sa_GFqPari(self, number *n, ring *_ring, object base)
    cdef extern object si2sa_NF(self, number *n, ring *_ring, object base)

    cdef extern number *sa2si_QQ(self, Rational ,ring (*))
    cdef extern number *sa2si_ZZ(self, Integer d, ring *_ring)

    cdef extern number *sa2si_GFqGivaro(self, int exp ,ring (*))
    cdef extern number *sa2si_GFqPari(self, object vector, ring *_ring)
    cdef extern number *sa2si_GFqNTLGF2E(self, FiniteField_ntl_gf2eElement elem, ring *_ring)

    cdef extern number *sa2si_NF(self, object element, ring *_ring)

    cdef extern object si2sa(self, number *n, ring *_ring, object base)
    cdef extern number *sa2si(self, Element elem, ring * _ring)


    cdef extern  MPolynomial_libsingular new_MP(self, MPolynomialRing_libsingular parent, poly *p)
