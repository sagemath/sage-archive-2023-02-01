include "../../ext/cdefs.pxi"

cimport sage.rings.padics.padic_base_generic_element
from sage.rings.padics.padic_base_generic_element cimport pAdicBaseGenericElement

cimport sage.structure.element
from sage.structure.element cimport CommutativeRingElement, RingElement, ModuleElement, Element

cimport sage.rings.padics.pow_computer
from sage.rings.padics.pow_computer cimport PowComputer_class

cimport sage.rings.integer
from sage.rings.integer cimport Integer

cdef class pAdicCappedAbsoluteElement(pAdicBaseGenericElement):
    cdef mpz_t value
    cdef unsigned long absprec

    cdef pAdicCappedAbsoluteElement _new_c(self)
    cdef RingElement _invert_c_impl(self)
    cdef ModuleElement _neg_c_impl(self)
    cdef pAdicCappedAbsoluteElement _lshift_c(pAdicCappedAbsoluteElement self, long shift)
    cdef pAdicCappedAbsoluteElement _rshift_c(pAdicCappedAbsoluteElement self, long shift)
    cdef Integer lift_c(pAdicCappedAbsoluteElement self)
    cdef object teichmuller_list(pAdicCappedAbsoluteElement self)
    cdef pAdicCappedAbsoluteElement unit_part_c(pAdicCappedAbsoluteElement self)
    cdef unsigned long valuation_c(self)
    cdef val_unit_c(self)
