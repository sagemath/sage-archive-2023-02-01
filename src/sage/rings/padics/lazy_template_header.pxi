"""
This file provides the declaration for the LazyElement class,
which collects common functionality for the different lazy p-adic
template classes.

It is included in padic_lazy_element.pxd and should be included
in any pxd file implementing lazy `p`-adics.

AUTHORS:

- Xavier Caruso (2021-02) -- initial version
"""

#*****************************************************************************
#       Copyright (C) 2021 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer cimport Integer
from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement

cdef enum expansion_mode:
    simple_mode, smallest_mode, teichmuller_mode

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

cdef class LazyElementWithDigits(LazyElement):
    cdef celement _digits


# Assignment

cdef class LazyElement_zero(LazyElement):
    pass

cdef class LazyElement_one(LazyElementWithDigits):
    pass

cdef class LazyElement_bound(LazyElement):
    cdef LazyElement _x

cdef class LazyElement_value(LazyElementWithDigits):
    cdef long _valuebound
    cdef long _shift
    cdef _value

cdef class LazyElement_random(LazyElementWithDigits):
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

cdef class LazyElement_add(LazyElementWithDigits):
    cdef LazyElement _x
    cdef LazyElement _y

cdef class LazyElement_sub(LazyElementWithDigits):
    cdef LazyElement _x
    cdef LazyElement _y

cdef class LazyElement_mul(LazyElementWithDigits):
    cdef LazyElement _x
    cdef cdigit _lastdigit_x
    cdef LazyElement _y
    cdef cdigit _lastdigit_y
    cdef int _update_last_digit(self)

cdef class LazyElement_muldigit(LazyElementWithDigits):
    cdef cdigit_ptr _x
    cdef LazyElement _y
    
cdef class LazyElement_div(LazyElementWithDigits):
    cdef long _maxprec
    cdef cdigit _inverse
    cdef LazyElement _num
    cdef LazyElement _denom
    cdef LazyElement _definition
    cdef int _bootstrap_c(self)
    cdef bint _bootstraping

cdef class LazyElement_sqrt(LazyElementWithDigits):
    cdef LazyElement _x
    cdef LazyElement _definition
    cdef int _bootstrap_c(self)

cdef class LazyElement_teichmuller(LazyElementWithDigits):
    cdef bint _ready
    cdef bint _trivial
    cdef list _xns
    cdef LazyElement _xbar
    cdef LazyElement _xp

# Self-referent numbers

cdef class LazyElement_selfref(LazyElementWithDigits):
    cdef LazyElement _definition
    cdef long _next
    cpdef set(self, LazyElement definition)
    # for pickling
    cdef long _initialvaluation
    cdef long _initialprecrel

# Expansion

cdef class ExpansionIter(object):
    cdef LazyElement elt
    cdef expansion_mode mode
    cdef long start
    cdef long stop
    cdef long current
    cdef cdigit digit
    cdef cdigit carry

    cdef void _next_simple(self)
    cdef void _next_smallest(self)
    cdef void _next_teichmuller(self)
