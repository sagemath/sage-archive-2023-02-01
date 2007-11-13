include "../../ext/cdefs.pxi"
include "../../libs/ntl/decl.pxi"

from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class

cdef class PowComputer_ext(PowComputer_class):
    cdef ZZ_c* small_powers
    cdef ZZ_c top_power
    cdef ZZ_c temp_z

    # Here for Eisenstein pow computers
    cdef int low_length
    cdef int high_length

    # the following three should be set by the subclasses
    cdef long ram_prec_cap # = prec_cap * e
    cdef long deg
    cdef long e
    cdef long f

    cdef ZZ_c* pow_ZZ_tmp(self, unsigned long n)
    cdef ZZ_c* pow_ZZ_top(self)

    cdef void cleanup_ext(self)

cdef class PowComputer_ZZ_pX(PowComputer_ext):
    cdef ntl_ZZ_pContext_class get_context(self, unsigned long n)
    cdef ntl_ZZ_pContext_class get_top_context(self)
    cdef restore_context(self, unsigned long n)
    cdef void restore_top_context(self)
    cdef ZZ_pX_Modulus_c* get_modulus(self, unsigned long n)
    cdef ZZ_pX_Modulus_c* get_top_modulus(self)

cdef class PowComputer_ZZ_pX_FM(PowComputer_ZZ_pX):
    #cdef ZZ_pX_c poly
    cdef ntl_ZZ_pContext_class c
    cdef ZZ_pX_Modulus_c mod

    cdef void cleanup_ZZ_pX_FM(self)

cdef class PowComputer_ZZ_pX_FM_Eis(PowComputer_ZZ_pX_FM):
    cdef ZZ_pX_Multiplier_c* low_shifter
    cdef ZZ_pX_Multiplier_c* high_shifter

    cdef void cleanup_ZZ_pX_FM_Eis(self)

cdef class PowComputer_ZZ_pX_small(PowComputer_ZZ_pX):
    cdef object c # using a python list so that we can store ntl_ZZ_pContext_class objects
    cdef ZZ_pX_Modulus_c *mod

    cdef void cleanup_ZZ_pX_small(self)

cdef class PowComputer_ZZ_pX_small_Eis(PowComputer_ZZ_pX_small):
    cdef ZZ_pX_Multiplier_c* low_shifter
    cdef ZZ_pX_Multiplier_c* high_shifter

    cdef void cleanup_ZZ_pX_small_Eis(self)

cdef class PowComputer_ZZ_pX_big(PowComputer_ZZ_pX):
    cdef object context_list # using a python list so that we can store ntl_ZZ_pContext_class objects
    cdef ZZ_pX_Modulus_c *modulus_list

    cdef ntl_ZZ_pContext_class top_context
    cdef ZZ_pX_Modulus_c top_mod

    cdef object context_dict #currently using a dict, optimize for speed later
    cdef object modulus_dict #currently using a dict, optimize for speed later

    cdef void cleanup_ZZ_pX_big(self)

cdef class PowComputer_ZZ_pX_big_Eis(PowComputer_ZZ_pX_big):
    cdef ZZ_pX_Multiplier_c* low_shifter
    cdef ZZ_pX_Multiplier_c* high_shifter

    cdef void cleanup_ZZ_pX_big_Eis(self)

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

#cdef int Eis_init(PowComputer_ZZ_pX prime_pow, ZZ_pX_Multiplier_c* low_shifter, ZZ_pX_Multiplier_c *high_shifter) except -1
