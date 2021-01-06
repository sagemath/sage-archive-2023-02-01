from sage.libs.flint.types cimport slong
from sage.libs.flint.types cimport flint_rand_t
from sage.libs.flint.types cimport fmpz, fmpz_t, fmpz_poly_t

from sage.rings.integer cimport Integer
from sage.rings.padics.pow_computer_flint cimport PowComputer_flint
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement



cpdef lazy_sum(parent, summands)


cdef class pAdicLazyElement(pAdicGenericElement):
    cdef PowComputer_flint prime_pow

    cdef fmpz_poly_t _digits
    cdef slong _valuation
    cdef slong _precrel

    #cdef pAdicLazyElement _new_c(self, type cls)

    cdef int _jump_c(self, slong prec)
    cdef int _next_c(self)
    cdef Integer _digit(self, slong i)

    cdef bint _is_equal(self, pAdicLazyElement right, slong prec, bint permissive) except -1



cdef class pAdicLazyElement_zero(pAdicLazyElement):
    pass

cdef class pAdicLazyElement_one(pAdicLazyElement):
    pass

cdef class pAdicLazyElement_value(pAdicLazyElement):
    cdef slong _maxprec
    cdef slong _shift
    cdef bint _finished

cdef class pAdicLazyElement_random(pAdicLazyElement):
    pass



cdef class pAdicLazyElement_shift(pAdicLazyElement):
    cdef pAdicLazyElement _x
    cdef slong _shift

cdef class pAdicLazyElement_add(pAdicLazyElement):
    cdef list _summands
    cdef list _signs

cdef class pAdicLazyElement_mul(pAdicLazyElement):
    cdef pAdicLazyElement _x
    cdef pAdicLazyElement _y

cdef class pAdicLazyElement_muldigit(pAdicLazyElement):
    cdef fmpz* _x
    cdef pAdicLazyElement _y
    
cdef class pAdicLazyElement_div(pAdicLazyElement):
    cdef slong _maxprec
    cdef fmpz_t _inverse
    cdef pAdicLazyElement _num
    cdef pAdicLazyElement _denom
    cdef pAdicLazyElement _definition
    cdef int _bootstrap_c(self)

cdef class pAdicLazyElement_sqrt(pAdicLazyElement):
    cdef slong _maxprec
    cdef pAdicLazyElement _x
    cdef pAdicLazyElement _definition
    cdef int _bootstrap_c(self)


cdef class pAdicLazyElement_selfref(pAdicLazyElement):
    cdef pAdicLazyElement _definition
    cdef bint _next

    cpdef set(self, pAdicLazyElement definition)
