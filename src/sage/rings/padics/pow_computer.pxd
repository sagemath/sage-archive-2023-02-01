include "../../ext/cdefs.pxi"

cimport sage.structure.sage_object
from sage.structure.sage_object cimport SageObject
cimport sage.rings.integer
from sage.rings.integer cimport Integer

cdef class PowComputer_class(SageObject):
    cdef Integer prime
    cdef unsigned long _cache_limit
    cdef mpz_t modulus
    cdef int in_field
    cdef mpz_t* dense_list
    cdef object dense_list_Integer
    cdef int _initialized
    cdef object __weakref__