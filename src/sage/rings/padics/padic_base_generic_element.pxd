include "../../ext/cdefs.pxi"

from padic_generic_element cimport pAdicGenericElement
from pow_computer cimport PowComputer_base

cdef class pAdicBaseGenericElement(pAdicGenericElement):
    cdef PowComputer_base prime_pow
    cdef int _set_to_mpz(self, mpz_t dest) except -1
    cdef int _set_to_mpq(self, mpq_t dest) except -1
    cdef int teichmuller_set_c(self, mpz_t value, mpz_t ppow) except -1