from sage.libs.flint.types cimport fmpz, fmpz_t, fmpz_poly_t, flint_rand_t

ctypedef fmpz_t cdigit
ctypedef fmpz* cdigit_ptr
ctypedef fmpz_poly_t celement
ctypedef flint_rand_t randgen


include "lazy_template_header.pxi"


cdef class pAdicLazyElement(LazyElement):
    pass

cdef class pAdicLazyElement_zero(LazyElement_zero):
    pass

cdef class pAdicLazyElement_one(LazyElement_one):
    pass

cdef class pAdicLazyElement_bound(LazyElement_bound):
    pass

cdef class pAdicLazyElement_value(LazyElement_value):
    pass

cdef class pAdicLazyElement_random(LazyElement_random):
    pass

cdef class pAdicLazyElement_slice(LazyElement_slice):
    pass

cdef class pAdicLazyElement_add(LazyElement_add):
    pass

cdef class pAdicLazyElement_sub(LazyElement_sub):
    pass

cdef class pAdicLazyElement_mul(LazyElement_mul):
    pass

cdef class pAdicLazyElement_muldigit(LazyElement_muldigit):
    pass

cdef class pAdicLazyElement_div(LazyElement_div):
    pass

cdef class pAdicLazyElement_sqrt(LazyElement_sqrt):
    pass

cdef class pAdicLazyElement_teichmuller(LazyElement_teichmuller):
    pass

cdef class pAdicLazyElement_selfref(LazyElement_selfref):
    pass
