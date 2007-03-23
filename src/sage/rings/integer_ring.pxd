include "../ext/cdefs.pxi"

from ring cimport PrincipalIdealDomain

cdef class IntegerRing_class(PrincipalIdealDomain):
    cdef int _randomize_mpz(self, mpz_t value, x, y, distribution) except -1
