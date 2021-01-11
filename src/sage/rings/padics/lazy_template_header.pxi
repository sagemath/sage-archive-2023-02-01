from sage.libs.flint.types cimport slong

from sage.rings.integer cimport Integer
from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement


cdef class LazyElement(pAdicGenericElement):
    cdef celement _digits
    cdef slong _valuation
    cdef slong _precrel
    cdef PowComputer_class prime_pow

    cdef int _jump_c(self, slong prec)
    cdef int _next_c(self)
    cdef Integer _digit(self, slong i)
    cdef bint _is_equal(self, LazyElement right, slong prec, bint permissive) except -1


cdef class LazyElement_zero(LazyElement):
    pass

cdef class LazyElement_one(LazyElement):
    pass

cdef class LazyElement_value(LazyElement):
    cdef slong _maxprec
    cdef slong _shift
    cdef bint _finished

cdef class LazyElement_random(LazyElement):
    pass



cdef class LazyElement_shift(LazyElement):
    cdef LazyElement _x
    cdef slong _shift

cdef class LazyElement_add(LazyElement):
    cdef list _summands
    cdef list _signs

cdef class LazyElement_mul(LazyElement):
    cdef LazyElement _x
    cdef cdigit _lastdigit_x
    cdef LazyElement _y
    cdef cdigit _lastdigit_y
    cdef int _update_last_digit(self)

cdef class LazyElement_muldigit(LazyElement):
    cdef cdigit_ptr _x
    cdef LazyElement _y
    cdef void _erase_first_digit(self)
    
cdef class LazyElement_div(LazyElement):
    cdef slong _maxprec
    cdef cdigit _inverse
    cdef LazyElement _num
    cdef LazyElement _denom
    cdef LazyElement _definition
    cdef int _bootstrap_c(self)

cdef class LazyElement_sqrt(LazyElement):
    cdef slong _maxprec
    cdef LazyElement _x
    cdef LazyElement _definition
    cdef int _bootstrap_c(self)

cdef class LazyElement_teichmuller(LazyElement):
    cdef bint _trivial
    cdef list _xns
    cdef LazyElement _xbar
    cdef LazyElement _xp


cdef class LazyElement_selfref(LazyElement):
    cdef LazyElement _definition
    cdef bint _next
    cpdef set(self, LazyElement definition)
