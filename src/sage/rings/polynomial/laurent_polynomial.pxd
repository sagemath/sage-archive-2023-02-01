
from sage.rings.quotient_ring_element import QuotientRingElement
from sage.structure.element cimport CommutativeAlgebraElement, ModuleElement, RingElement, Element
from sage.rings.polynomial.polydict cimport ETuple, PolyDict
from sage.rings.polynomial.multi_polynomial cimport MPolynomial
from sage.rings.polynomial.polynomial_element cimport Polynomial

cdef class LaurentPolynomial_quotient(QuotientRingElement):
    cdef _new_c(self)

cdef class LaurentPolynomial_mpair(CommutativeAlgebraElement):
    cdef ETuple _mon
    cdef MPolynomial _poly
    cdef PolyDict _prod
    cdef _new_c(self)

#cdef class LaurentPolynomial_upair(CommutativeAlgebraElement):
#    cdef long _mon
#    cdef Polynomial _poly
#    cdef PolyDict _prod
#    cdef _new_c(self)