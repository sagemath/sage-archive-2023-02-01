include "singular-cdefs.pxi"

from sage.rings.rational cimport Rational
from sage.structure.element cimport Element
from sage.rings.integer cimport Integer
from sage.rings.integer_mod cimport IntegerMod_abstract
from sage.rings.finite_field_givaro cimport FiniteField_givaro,FiniteField_givaroElement
from sage.rings.finite_field_ntl_gf2e cimport FiniteField_ntl_gf2e,FiniteField_ntl_gf2eElement
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular

from sage.rings.number_field.number_field_base cimport NumberField

# Singular -> Sage

cdef Rational                    si2sa_QQ(number (*),ring (*))
cdef Integer                     si2sa_ZZ(number (*),ring (*))

cdef FiniteField_givaroElement   si2sa_GFqGivaro(number *n, ring *_ring, FiniteField_givaro base)
cdef FiniteField_ntl_gf2eElement si2sa_GFqNTLGF2E(number *n, ring *_ring, FiniteField_ntl_gf2e base)
cdef object                      si2sa_GFqPari(number *n, ring *_ring, object base)
cdef inline object               si2sa_ZZmod(number *n, ring *_ring, object base)

cdef object                       si2sa_NF(number *n, ring *_ring, object base)

cdef object si2sa(number *n, ring *_ring, object base)

# Sage -> Singular

cdef number *sa2si_QQ(Rational ,ring (*))
cdef number *sa2si_ZZ(Integer d, ring *_ring)

cdef number *sa2si_GFqGivaro(int exp ,ring (*))
cdef number *sa2si_GFqPari(object vector, ring *_ring)
cdef number *sa2si_GFqNTLGF2E(FiniteField_ntl_gf2eElement elem, ring *_ring)
cdef inline number *sa2si_ZZmod(IntegerMod_abstract d, ring *_ring)

cdef number *sa2si_NF(object element, ring *_ring)

cdef number *sa2si(Element elem, ring * _ring)
