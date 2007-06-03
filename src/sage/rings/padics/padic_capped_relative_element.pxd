include "../../ext/cdefs.pxi"

cimport sage.rings.padics.padic_base_generic_element
from sage.rings.padics.padic_base_generic_element cimport pAdicBaseGenericElement

cimport sage.structure.element
from sage.structure.element cimport CommutativeRingElement, RingElement, ModuleElement, Element

cimport sage.rings.integer
from sage.rings.integer cimport Integer

cimport sage.rings.padics.pow_computer
from sage.rings.padics.pow_computer cimport PowComputer_class

cdef class pAdicCappedRelativeElement(pAdicBaseGenericElement):
    cdef mpz_t unit    #An exact zero is indicated by unit < 0
    cdef long ordp
    cdef long relprec
    #cdef PowComputer_class prime_pow
    cdef bint _normalized
    cdef void set_exact_zero(pAdicCappedRelativeElement self)
    cdef void set_inexact_zero(pAdicCappedRelativeElement self, long absprec)
    cdef void set_precs(pAdicCappedRelativeElement self, long relprec)
    cdef void set_from_Integers(pAdicCappedRelativeElement self, Integer ordp, Integer unit, Integer relprec)
    cdef pAdicCappedRelativeElement _new_c(pAdicCappedRelativeElement self)
    cdef void _normalize(pAdicCappedRelativeElement self)
    cdef ModuleElement _neg_c_impl(self)
    cdef RingElement _floordiv_c_impl(self, RingElement right)
    cdef pAdicCappedRelativeElement _lshift_c(pAdicCappedRelativeElement self, long shift)
    cdef pAdicCappedRelativeElement _rshift_c(pAdicCappedRelativeElement self, long shift)
    cdef pAdicCappedRelativeElement lift_to_precision_c(pAdicCappedRelativeElement self, long absprec)
    cdef teichmuller_list(pAdicCappedRelativeElement self)
    cdef val_unit_c(self)
    cdef pAdicCappedRelativeElement unit_part_c(self)
    cdef valuation_c(self)
    cdef lift_c(self)

