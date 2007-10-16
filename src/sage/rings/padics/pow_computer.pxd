include "../../ext/cdefs.pxi"

cimport sage.structure.sage_object
from sage.structure.sage_object cimport SageObject
cimport sage.rings.integer
from sage.rings.integer cimport Integer

cdef class PowComputer_class(SageObject):
    cdef Integer prime
    cdef bint in_field
    cdef int _initialized
    cdef object __weakref__

    cdef Integer pow_Integer(self, unsigned long n)
    cdef mpz_t pow_mpz_t(self, unsigned long n)
    cdef ZZ_c pow_ZZ(self, unsigned long n)

cdef class PowComputer_base(PowComputer_class):
    cdef mpz_t* small_powers
    cdef unsigned long cache_limit
    cdef mpz_t top_power
    cdef unsigned long prec_cap
    cdef mpz_t temp
