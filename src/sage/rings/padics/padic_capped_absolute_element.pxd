include "../../ext/cdefs.pxi"

from sage.rings.padics.padic_base_generic_element cimport pAdicBaseGenericElement
from sage.structure.element cimport CommutativeRingElement, RingElement, ModuleElement, Element
from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational

cdef class pAdicCappedAbsoluteElement(pAdicBaseGenericElement):
    cdef mpz_t value
    cdef unsigned long absprec

    cdef pAdicCappedAbsoluteElement _new_c(self)
    cpdef RingElement _invert_c_impl(self)
    cpdef ModuleElement _neg_(self)
    cdef pAdicCappedAbsoluteElement _lshift_c(pAdicCappedAbsoluteElement self, long shift)
    cdef pAdicCappedAbsoluteElement _rshift_c(pAdicCappedAbsoluteElement self, long shift)
    cdef object teichmuller_list(pAdicCappedAbsoluteElement self)
    cpdef pAdicCappedAbsoluteElement unit_part(self)
    cdef long valuation_c(self)
    cpdef val_unit(self)
    cpdef Integer lift(self)
