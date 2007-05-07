include "../../ext/cdefs.pxi"

cimport sage.structure.sage_object
from sage.structure.sage_object cimport SageObject
cimport sage.rings.integer
from sage.rings.integer cimport Integer

cdef class PowComputer_class(SageObject):
    cdef Integer basex
    cdef unsigned int log_of_dense_limit
    cdef unsigned int dense_mask
    cdef mpz_t dense_mask_mpz
    cdef mpz_t* dense_list
    cdef object dense_list_Integer
    cdef int _initialized
    cdef Integer zero
    cdef object _cache
    cdef object __weakref__
    cdef void pow_mpz_ui(self, mpz_t ans, unsigned int n)
    cdef void pow_mpz_mpz(self, mpz_t ans, mpz_t n)
