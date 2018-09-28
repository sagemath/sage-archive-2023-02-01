from sage.structure.element cimport MonoidElement
from sage.structure.element cimport CommutativeAlgebraElement

from sage.rings.polynomial.polydict cimport PolyDict
from sage.rings.polynomial.polydict cimport ETuple

from sage.rings.padics.padic_generic_element cimport pAdicGenericElement


cdef class TateAlgebraTerm(MonoidElement):
    cdef _field
    cdef pAdicGenericElement _coeff
    cdef ETuple _exponent

    cdef TateAlgebraTerm _new_c(self)
    cdef long _valuation_c(self)
    cpdef TateAlgebraTerm monomial(self)
    cpdef TateAlgebraTerm monic(self)
    cdef TateAlgebraTerm _gcd_c(self, TateAlgebraTerm other)
    cdef TateAlgebraTerm _lcm_c(self, TateAlgebraTerm other)
    cdef bint _divides_c(self, TateAlgebraTerm other, bint integral)
    cdef TateAlgebraTerm _floordiv_c(self, TateAlgebraTerm other)


cdef class TateAlgebraElement(CommutativeAlgebraElement):
    cdef _prec
    cdef PolyDict _poly
    cdef list _terms
    cdef bint _is_normalized

    cdef _normalize(self)
    cdef TateAlgebraElement _new_c(self)
    cdef list _terms_c(self)
    cpdef valuation(self)
    cdef TateAlgebraElement _term_mul_c(self, TateAlgebraTerm term)
    cdef TateAlgebraElement _positive_lshift_c(self, n)
    cdef TateAlgebraElement _lshift_c(self, n)
    cpdef TateAlgebraElement monic(self)
    cdef _quo_rem_c(self, list divisors, bint quo, bint rem, bint integral)
    cdef _quo_rem_check(self, divisors, bint quo, bint rem)
    cdef TateAlgebraElement _Spoly_c(self, TateAlgebraElement other)

