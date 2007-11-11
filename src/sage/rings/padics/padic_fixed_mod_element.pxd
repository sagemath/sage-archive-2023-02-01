include "../../ext/cdefs.pxi"

cimport sage.rings.padics.padic_base_generic_element
from sage.rings.padics.padic_base_generic_element cimport pAdicBaseGenericElement

cimport sage.structure.element
from sage.structure.element cimport CommutativeRingElement, RingElement, ModuleElement, Element

cimport sage.rings.integer
from sage.rings.integer cimport Integer

cimport sage.rings.padics.pow_computer
from sage.rings.padics.pow_computer cimport PowComputer_class

cdef class pAdicFixedModElement(pAdicBaseGenericElement):
    cdef mpz_t value
    cdef int _set_from_mpz(pAdicFixedModElement self, mpz_t value) except -1
    cdef int _set_from_mpq(pAdicFixedModElement self, mpq_t x) except -1
    cdef int _set_to_mpz(pAdicFixedModElement self, mpz_t dest) except -1
    cdef int _set_to_mpq(pAdicFixedModElement self, mpq_t dest) except -1
    cdef pAdicFixedModElement _new_c(self)
    cdef RingElement _invert_c_impl(self)
    cdef pAdicFixedModElement _lshift_c(pAdicFixedModElement self, long shift)
    cdef pAdicFixedModElement _rshift_c(pAdicFixedModElement self, long shift)
    cdef ModuleElement _neg_c_impl(self)
    cdef Integer lift_c(pAdicFixedModElement self)
    cdef object teichmuller_list(pAdicFixedModElement self)
    cdef pAdicFixedModElement unit_part_c(pAdicFixedModElement self)
    cdef long valuation_c(self)
    cdef val_unit_c(self)
