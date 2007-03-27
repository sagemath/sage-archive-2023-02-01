include "../ext/cdefs.pxi"
include "../libs/ntl/decl.pxi"

from ring cimport PrincipalIdealDomain
from integer cimport Integer

cdef class IntegerRing_class(PrincipalIdealDomain):
    cdef Integer _coerce_ZZ(self, ntl_c_ZZ *z)
    cdef int _randomize_mpz(self, mpz_t value, x, y, distribution) except -1
