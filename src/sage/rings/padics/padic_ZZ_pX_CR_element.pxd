include "../../libs/ntl/decl.pxi"

from sage.rings.padics.padic_ZZ_pX_element cimport pAdicZZpXElement
from sage.structure.element cimport RingElement, ModuleElement

cdef class pAdicZZpXCRElement(pAdicZZpXElement):
    cdef ZZ_pX_c unit
    cdef long ordp
    cdef long relprec
    cdef bint _normalized
    cdef pAdicZZpXCRElement _new_c(self)
    cdef pAdicZZpXCRElement _lshift_c(self, long n, bint top_zeros)
    cdef pAdicZZpXCRElement _rshift_c(self, long n, bint top_zeros)
    cdef pAdicZZpXCRElement _internal_shift(self, long n, bint top_zeros)
    cdef long _internal_valuation(self)
    cdef pAdicZZpXCRElement _unit_part_c(self, bint top_zeros)
    cdef RingElement _invert_c_impl(self)

