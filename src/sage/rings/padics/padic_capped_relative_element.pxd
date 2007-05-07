include "../../ext/cdefs.pxi"

cimport sage.rings.padics.padic_base_generic_element
from sage.rings.padics.padic_base_generic_element cimport pAdicBaseGenericElement

cimport sage.structure.element
from sage.structure.element cimport CommutativeRingElement, RingElement, ModuleElement, Element

cimport sage.rings.integer
from sage.rings.integer cimport Integer

cdef class pAdicCappedRelativeElement(pAdicBaseGenericElement):
    cdef mpz_t unit #An exact zero is indicated by unit < 0
    cdef mpz_t ordp
    cdef mpz_t modulus
    cdef mpz_t p
    cdef int in_field
    cdef int relprec
    cdef void set_exact_zero(pAdicCappedRelativeElement self)
    cdef void set_inexact_zero(pAdicCappedRelativeElement self, mpz_t absprec)
    cdef void set_precs(pAdicCappedRelativeElement self, unsigned int relprec)
    cdef void set_from_Integers(pAdicCappedRelativeElement self, Integer ordp, Integer unit, Integer relprec)
    cdef pAdicCappedRelativeElement _new_c(pAdicCappedRelativeElement self)
    cdef ModuleElement _neg_c_impl(self)
    cdef ModuleElement _add_c_impl(self, ModuleElement right)
    cdef RingElement _invert_c_impl(self)
    cdef RingElement _floordiv_c_impl(self, RingElement right)
    cdef pAdicCappedRelativeElement _lshift_c(pAdicCappedRelativeElement self, mpz_t shift)
    cdef pAdicCappedRelativeElement _rshift_c(pAdicCappedRelativeElement self, mpz_t shift)
    cdef RingElement _mul_c_impl(self, RingElement right)
    cdef RingElement _div_c_impl(self, RingElement right)
    cdef pAdicCappedRelativeElement lift_to_precision_c(pAdicCappedRelativeElement self, mpz_t absprec)
    cdef teichmuller_list(pAdicCappedRelativeElement self)
    cdef val_unit_c(self)
    cdef pAdicCappedRelativeElement unit_part_c(self)
    cdef valuation_c(self)
    cdef long _hash(self) except -1

