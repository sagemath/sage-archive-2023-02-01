from sage.structure.element cimport MonoidElement
from sage.structure.element cimport CommutativeAlgebraElement
#from sage.rings.integer cimport Integer

cdef class TateAlgebraTerm(MonoidElement):
    cdef _field
    cdef _coeff
    cdef tuple _exponent

cdef class TateAlgebraElement(CommutativeAlgebraElement):
    cdef _poly
    cdef _prec
    cdef list _coeffs
