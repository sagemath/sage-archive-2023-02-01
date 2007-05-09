include "../../ext/cdefs.pxi"

cimport sage.rings.padics.padic_base_generic_element
from sage.rings.padics.padic_base_generic_element cimport pAdicBaseGenericElement

cimport sage.structure.element
from sage.structure.element cimport CommutativeRingElement, RingElement, ModuleElement, Element

cimport sage.rings.padics.pow_computer
from sage.rings.padics.pow_computer cimport PowComputer_class

cimport sage.rings.integer
from sage.rings.integer cimport Integer

cdef class pAdicRingCappedAbsoluteElement(pAdicBaseGenericElement):
    cdef mpz_t value
    cdef mpz_t p
    cdef mpz_t modulus
    cdef int absprec
    cdef void set_precs(pAdicRingCappedAbsoluteElement self, unsigned int absprec)
    cdef void set_value_from_mpz(pAdicRingCappedAbsoluteElement self, mpz_t value)
    cdef void set_from_Integers(pAdicRingCappedAbsoluteElement self, Integer value, absprec)
    cdef pAdicRingCappedAbsoluteElement _new_c(self)
    cdef RingElement _invert_c_impl(self)
    cdef ModuleElement _neg_c_impl(self)
    cdef ModuleElement _add_c_impl(self, ModuleElement right)
    cdef ModuleElement _sub_c_impl(self, ModuleElement right)
    cdef RingElement _div_c_impl(self, RingElement right)
    cdef RingElement _mul_c_impl(self, RingElement right)
    cdef pAdicRingCappedAbsoluteElement _lshift_c(pAdicRingCappedAbsoluteElement self, int shift)
    cdef pAdicRingCappedAbsoluteElement _rshift_c(pAdicRingCappedAbsoluteElement self, int shift)
    cdef Integer lift_c(pAdicRingCappedAbsoluteElement self)
    cdef object teichmuller_list(pAdicRingCappedAbsoluteElement self)
    cdef pAdicRingCappedAbsoluteElement unit_part_c(pAdicRingCappedAbsoluteElement self)
    cdef int valuation_c(self)
    cdef val_unit_c(self)
    cdef long _hash(self) except -1