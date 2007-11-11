include "../../ext/cdefs.pxi"

from sage.rings.padics.padic_generic_element cimport pAdicGenericElement

cdef class pAdicExtElement(pAdicGenericElement):
    cdef ext_p_list(self, bint pos)
    cdef int _set_from_mpz(self, mpz_t x) except -1
    cdef int _set_from_mpq(self, mpq_t x) except -1