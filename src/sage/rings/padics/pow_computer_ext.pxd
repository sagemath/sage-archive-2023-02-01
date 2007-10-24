include "../../libs/ntl/decl.pxi""

from sage.rings.padics.pow_computer cimport PowComputer_class

cdef class PowComputer_ext(PowComputer_class):
    cdef ZZ_c* small_powers
    cdef unsigned long cache_limit
    cdef ZZ_c top_power
    cdef unsigned long prec_cap

cdef class PowComputer_ZZ_pX(PowComputer_ext):
    cdef ntl_ZZ_pContext get_context(self, unsigned long n)
    cdef ntl_ZZ_pContext get_top_context(self)
    cdef void restore_context(self, unsigned long n)
    cdef void restore_top_context(self)
    cdef ZZ_pX_Modulus_c get_modulus(self, unsigned long n)
    cdef ZZ_pX_Modulus_c get_top_modulus(self)

cdef class PowComputer_ZZ_pX_FM(PowComputer_ZZ_pX):
    cdef ZZ_pX_c poly
    cdef ntl_ZZ_pContext c
    cdef ZZ_pX_Modulus_c mod

cdef class PowComputer_ZZ_pX_FM_Eis(PowComputer_ZZ_pX_FM):
    pass

#cdef class PowComputer_ZZ_pX_small(PowComputer_ZZ_pX):
#    cdef ZZ_pX_c *poly
#    cdef ntl_ZZ_pContext *c
#    cdef ZZ_pX_Modulus_c *mod
#
#cdef class PowComputer_ZZ_pX_small_Eis(PowComputer_ZZ_pX_small):
#    pass
#
#cdef class PowComputer_ZZ_pX_big(PowComputer_ZZ_pX):
#    cdef ZZX_c poly
#    cdef context_dict #currently using a dict, optimize for speed later
#    cdef modulus_dict #currently using a dict, optimize for speed later
#
#cdef class PowComputer_ZZ_pX_big_Eis(PowComputer_ZZ_pX_big):
#    pass
#
#cdef class PowComputer_ZZ_pEX(PowComputer_ext):
#    cdef ntl_ZZ_pEContext get_context(self, unsigned long n)
#    cdef ntl_ZZ_pEContext get_top_context(self)
#    cdef void restore_context(self, unsigned long n)
#    cdef void restore_top_context(self)
#    cdef ZZ_pEX_Modulus_c get_modulus(self, unsigned long n)
#    cdef ZZ_pEX_Modulus_c get_top_modulus(self)
#
#cdef class PowComputer_ZZ_pEX_FM(PowComputer_ZZ_pEX):
#    cdef ZZ_pEX_c poly
#    cdef ntl_ZZ_pEContext c
#    cdef ZZ_pEX_Modulus_c mod
#
#cdef class PowComputer_ZZ_pEX_FM_Eis(PowComputer_ZZ_pEX_FM):
#    pass
#
#cdef class PowComputer_ZZ_pEX_small(PowComputer_ZZ_pEX):
#    cdef ZZ_pEX_c *poly
#    cdef ntl_ZZ_pEContext *c
#    cdef ZZ_pEX_Modulus_c *mod
#
#cdef class PowComputer_ZZ_pEX_small_Eis(PowComputer_ZZ_pEX_small):
#    pass
#
#cdef class PowComputer_ZZ_pEX_big(PowComputer_ZZ_pEX):
#    cdef ZZ_pEX_c poly
#    cdef context_dict #currently using a dict, optimize for speed later
#    cdef modulus_dict #currently using a dict, optimize for speed later
#
#cdef class PowComputer_ZZ_pEX_big_Eis(PowComputer_ZZ_pEX_big):
#    pass

