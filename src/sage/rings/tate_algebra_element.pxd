from sage.structure.element cimport MonoidElement
from sage.structure.element cimport CommutativeAlgebraElement

from sage.rings.polynomial.polydict cimport ETuple


cdef class TateAlgebraTerm(MonoidElement):
    cdef _field
    cdef _coeff
    cdef ETuple _exponent

cdef class TateAlgebraElement(CommutativeAlgebraElement):
    cdef _poly
    cdef _prec
    cdef list _terms
