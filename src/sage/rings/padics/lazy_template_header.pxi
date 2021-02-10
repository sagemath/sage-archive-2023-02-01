from sage.rings.integer cimport Integer
from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement


cdef class LazyElement(pAdicGenericElement):
    cdef long _valuation
    cdef long _precrel
    cdef long _precbound
    cdef PowComputer_class prime_pow

    cdef cdigit_ptr _getdigit_relative(self, long i)
    cdef cdigit_ptr _getdigit_absolute(self, long i)
    cdef void _getslice_relative(self, celement slice, long start, long length)

    cdef int _init_jump(self) except -1
    cdef int _jump_c(self, long prec)
    cdef int _jump_relative_c(self, long prec, long halt)
    cdef int _next_c(self)
    cdef bint _is_equal(self, LazyElement right, long prec, bint permissive) except -1

cdef class LazyElement_abandon(LazyElement):
    pass
cdef lazyelement_abandon

cdef class LazyElement_init(LazyElement):
    cdef celement _digits


# Assignment

cdef class LazyElement_zero(LazyElement):
    pass

cdef class LazyElement_one(LazyElement_init):
    pass

cdef class LazyElement_bound(LazyElement):
    cdef LazyElement _x

cdef class LazyElement_value(LazyElement_init):
    cdef long _valuebound
    cdef long _shift
    cdef _value

cdef class LazyElement_random(LazyElement_init):
    cdef randgen _generator
    # for pickling
    cdef long _initialvaluation
    cdef long _seed


# Operations

cdef class LazyElement_slice(LazyElement):
    cdef LazyElement _x
    cdef long _start
    cdef long _stop
    cdef long _shift

cdef class LazyElement_add(LazyElement_init):
    cdef LazyElement _x
    cdef LazyElement _y

cdef class LazyElement_sub(LazyElement_init):
    cdef LazyElement _x
    cdef LazyElement _y

cdef class LazyElement_mul(LazyElement_init):
    cdef LazyElement _x
    cdef cdigit _lastdigit_x
    cdef LazyElement _y
    cdef cdigit _lastdigit_y
    cdef int _update_last_digit(self)

cdef class LazyElement_muldigit(LazyElement_init):
    cdef cdigit_ptr _x
    cdef LazyElement _y
    cdef void _erase_first_digit(self)
    
cdef class LazyElement_div(LazyElement_init):
    cdef long _maxprec
    cdef cdigit _inverse
    cdef LazyElement _num
    cdef LazyElement _denom
    cdef LazyElement _definition
    cdef int _bootstrap_c(self)
    cdef bint _bootstraping

cdef class LazyElement_sqrt(LazyElement_init):
    cdef LazyElement _x
    cdef LazyElement _definition
    cdef int _bootstrap_c(self)

cdef class LazyElement_teichmuller(LazyElement_init):
    cdef bint _ready
    cdef bint _trivial
    cdef list _xns
    cdef LazyElement _xbar
    cdef LazyElement _xp

# Self-referent numbers

cdef class LazyElement_selfref(LazyElement_init):
    cdef LazyElement _definition
    cdef long _next
    cpdef set(self, LazyElement definition)
    # for pickling
    cdef long _initialvaluation
    cdef long _initialprecrel
