from sage.libs.flint.types cimport slong
from sage.libs.flint.types cimport flint_rand_t
from sage.libs.flint.types cimport fmpz, fmpz_t, fmpz_poly_t

from sage.rings.integer cimport Integer
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement


cdef class pAdicLazyElement(pAdicGenericElement):
    cdef fmpz_t _prime
    cdef fmpz_poly_t _digits
    cdef slong _valuation
    cdef slong _precrel

    cpdef jump(self, prec)
    cdef bint next(self) except -1
    cdef Integer _digit(self, i)

    cdef bint _is_equal(self, pAdicLazyElement right, slong prec)



cdef class pAdicLazyElement_zero(pAdicLazyElement):
    pass

cdef class pAdicLazyElement_one(pAdicLazyElement):
    pass

cdef class pAdicLazyElement_value(pAdicLazyElement):
    cdef _maxprec

cdef class pAdicLazyElement_random(pAdicLazyElement):
    cdef flint_rand_t _randstate



cdef class pAdicLazyElement_shift(pAdicLazyElement):
    cdef pAdicLazyElement _x
    cdef slong _shift

cdef class pAdicLazyElement_add(pAdicLazyElement):
    cdef list _summands
    cdef list _signs

cdef class pAdicLazyElement_mul(pAdicLazyElement):
    cdef pAdicLazyElement _x
    cdef pAdicLazyElement _y
    cdef fmpz_poly_t _tmp
    
#cdef class pAdicLazyElement_div(pAdicLazyElement)
#cdef class pAdicLazyElement_sqrt(pAdicLazyElement)



cdef class pAdicLazyElement_selfref(pAdicLazyElement):
    cdef pAdicLazyElement _definition
    cdef bint _next
