include "../../ext/cdefs.pxi"

cimport sage.rings.padics.local_generic_element
from sage.rings.padics.local_generic_element cimport LocalGenericElement
cimport sage.structure.element
from sage.structure.element cimport Element
cimport sage.rings.padics.pow_computer
from sage.rings.padics.pow_computer cimport PowComputer_class

cdef class pAdicGenericElement(LocalGenericElement):
    cdef PowComputer_class prime_pow
    cdef int _cmp_c_impl(left, Element right) except -2
    cdef base_p_list(self, mpz_t value, lift_mode)

cdef extern void teichmuller_set_c(mpz_t value, mpz_t p, mpz_t ppow)
