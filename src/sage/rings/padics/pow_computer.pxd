include "../../ext/cdefs.pxi"
include "../../libs/ntl/decl.pxi"

cimport sage.structure.sage_object
from sage.structure.sage_object cimport SageObject
cimport sage.rings.integer
from sage.rings.integer cimport Integer

cdef class PowComputer_class(SageObject):
    cdef Integer prime
    cdef bint in_field
    cdef int _initialized

    cdef unsigned long cache_limit
    cdef unsigned long prec_cap

    cdef Integer pow_Integer(self, unsigned long n)
    cdef mpz_t* pow_mpz_top(self)
    cdef mpz_t* pow_mpz_t(self, unsigned long n)
    cdef mpz_t* pow_mpz_t_tmp(self, unsigned long n)
    cdef ZZ_c pow_ZZ(self, unsigned long n)

cdef class PowComputer_base(PowComputer_class):
    cdef mpz_t* small_powers
    cdef mpz_t top_power
    cdef mpz_t temp
    cdef object __weakref__