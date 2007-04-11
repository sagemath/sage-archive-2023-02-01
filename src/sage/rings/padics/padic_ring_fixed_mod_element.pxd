include "../../ext/cdefs.pxi"

cimport sage.rings.padics.padic_base_generic_element
from sage.rings.padics.padic_base_generic_element cimport pAdicBaseGenericElement

cimport sage.structure.element
from sage.structure.element cimport CommutativeRingElement, RingElement, ModuleElement

cimport sage.rings.integer
from sage.rings.integer cimport Integer

cdef class pAdicRingFixedModElement(pAdicBaseGenericElement):
    cdef mpz_t value
    cdef mpz_t modulus
    cdef mpz_t p
    cdef void set_from_mpz(pAdicRingFixedModElement self, mpz_t value)
    cdef pAdicRingFixedModElement _new_c(self)
    cdef RingElement _invert_c_impl(self)
    cdef pAdicRingFixedModElement _lshift_c(pAdicRingFixedModElement self, int shift)
    cdef pAdicRingFixedModElement _rshift_c(pAdicRingFixedModElement self, int shift)
    cdef ModuleElement _neg_c_impl(self)
    cdef ModuleElement _add_c_impl(self, ModuleElement right)
    cdef ModuleElement _sub_c_impl(self, ModuleElement right)
    cdef RingElement _mul_c_impl(self, RingElement right)
    cdef RingElement _div_c_impl(self, RingElement right)
    cdef Integer lift_c(pAdicRingFixedModElement self)
    cdef object teichmuller_list(pAdicRingFixedModElement self)
    cdef pAdicRingFixedModElement unit_part_c(pAdicRingFixedModElement self)
    cdef int valuation_c(self)
