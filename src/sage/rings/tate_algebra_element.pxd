from sage.structure.element cimport RingElement
from sage.structure.element cimport MonoidElement
from sage.structure.element cimport CommutativeAlgebraElement

from sage.rings.polynomial.polydict cimport PolyDict
from sage.rings.polynomial.polydict cimport ETuple


cdef class TateAlgebraTerm(MonoidElement):
    cdef _field
    cdef RingElement _coeff
    cdef ETuple _exponent

    cdef TateAlgebraTerm _new_c(self)
    cdef long _valuation_c(self)
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
    cdef _term_mul_c(self, TateAlgebraTerm term)
    cdef _positive_lshift_c(self, n)
    cdef _lshift_c(self, n)
    cpdef monic(self)
    cdef TateAlgebraElement _mod_c(TateAlgebraElement self, list divisors)
    cdef TateAlgebraElement _Spoly_c(self, TateAlgebraElement other)

