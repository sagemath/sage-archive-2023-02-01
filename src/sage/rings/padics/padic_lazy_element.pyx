# ****************************************************************************
#       Copyright (C) 2021 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************


from sage.libs.flint.types cimport *
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_poly cimport *

cdef extern from "sage/rings/padics/padic_lazy_element_helper.c":
    cdef void flint_randinit(flint_rand_t state)
    cdef void flint_randclear(flint_rand_t state)
    cdef fmpz* get_coeff(fmpz_poly_t poly, slong i)
    cdef void get_slice(fmpz_poly_t slice, fmpz_poly_t poly, slong start, slong length)
    cdef void iadd_coeff(fmpz_poly_t poly, fmpz_t summand, slong i)
    cdef void isub_coeff(fmpz_poly_t poly, fmpz_t summand, slong i)
    cdef void iadd_shifted(fmpz_poly_t poly, fmpz_poly_t summand, slong shift)
    cdef void reduce_coeff(fmpz_poly_t poly, slong i, fmpz_t modulus)

from sage.ext.stdsage cimport PY_NEW
from sage.structure.element import coerce_binop

from sage.rings.all import ZZ
from sage.rings.integer cimport Integer
from sage.rings.infinity import Infinity

from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.rings.padics.precision_error import PrecisionError

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1
MAXORDP = ZZ(maxordp)


cdef class pAdicLazyElement(pAdicGenericElement):
    def __cinit__(self):
        fmpz_poly_init(self._digits)

    def __dealloc__(self):
        fmpz_poly_clear(self._digits)

    def __init__(self, parent):
        pAdicGenericElement.__init__(self, parent)
        cdef Integer p = self._parent.prime()
        fmpz_set_mpz(self._prime, p.value)
        self._valuation = 0
        self._precrel = 0

    def prime(self):
        return self._parent.prime()

    cdef bint next(self) except -1:
        raise NotImplementedError("must be implemented in subclasses")

    cpdef jump(self, prec):
        cdef slong i, stop
        stop = min(prec, maxordp) - self._precrel - self._valuation
        for i in range(stop):
            if not self.next():
                return False
        return prec < maxordp

    def expansion(self, n=None):
        if n is None:
            return ExpansionIter(self, self._valuation, self._valuation + self._precrel)
        if isinstance(n, slice):
            if n.step is not None:
                raise NotImplementedError
            return ExpansionIter(self, n.start, n.stop)
        else:
            return self._digit(n)

    cdef Integer _digit(self, i):
        if not self.jump(i+1):
            raise PrecisionError
        cdef Integer ans = PY_NEW(Integer)
        cdef fmpz* coeff = get_coeff(self._digits, i - self._valuation)
        fmpz_get_mpz(ans.value, coeff)
        return ans

    def _repr_(self):
        self.jump(self._parent.default_prec())
        return pAdicGenericElement._repr_(self)

    cdef bint _is_equal(self, pAdicLazyElement right, slong prec):
        cdef slong i
        self.jump(prec)
        right.jump(prec)
        if self._valuation >= prec:
            return right._valuation >= prec
        if right._valuation >= prec:
            return False
        if self._valuation != right._valuation:
            return False
        for i in range(prec - self._valuation):
            if not fmpz_equal(get_coeff(self._digits, i), get_coeff(right._digits, i)):
                return False
        return True

    @coerce_binop
    def is_equal_at_precision(self, other, prec=None):
        if prec is None:
            prec = min(self.precision_absolute(), other.precision_absolute())
            prec = max(prec, self._parent.default_prec())
        return self._is_equal(other, long(prec))

    def __eq__(self, other):
        return self.is_equal_at_precision(other)

    cpdef bint _is_exact_zero(self) except -1:
        return isinstance(self, pAdicLazyElement_zero)

    cpdef bint _is_inexact_zero(self) except -1:
        return self._precrel == 0

    def is_exact(self):
        return self._is_exact_zero()

    def precision_absolute(self):
        if self._is_exact_zero():
            return Infinity
        return ZZ(self._precrel + self._valuation)

    def precision_relative(self):
        return ZZ(self._precrel)

    def valuation(self):
        if self._is_exact_zero():
            return Infinity
        return ZZ(self._valuation)

    def unit_part(self):
        if self._precrel == 0:
            raise PrecisionError
        return self >> self._valuation

    def __lshift__(self, s):
        cdef slong shift = long(s)
        if shift:
            return pAdicLazyElement_shift((<pAdicLazyElement>self)._parent, self, -shift)
        else:
            return self

    def __rshift__(self, s):
        cdef slong shift = long(s)
        if shift:
            return pAdicLazyElement_shift((<pAdicLazyElement>self)._parent, self, shift)
        else:
            return self

    cpdef _add_(self, other):
        cdef list summands
        cdef list signs
        if isinstance(self, pAdicLazyElement_zero):
            return other
        if isinstance(other, pAdicLazyElement_zero):
            return self
        if isinstance(self, pAdicLazyElement_add):
            summands = list((<pAdicLazyElement_add>self)._summands)
            signs = list((<pAdicLazyElement_add>self)._signs)
        else:
            summands = [self]
            signs = [True]
        if isinstance(other, pAdicLazyElement_add):
            summands.extend((<pAdicLazyElement_add>other)._summands)
            signs.extend((<pAdicLazyElement_add>other)._signs)
        else:
            summands.append(other)
            signs.append(True)
        return pAdicLazyElement_add(self._parent, summands, signs)

    cpdef _sub_(self, other):
        cdef list summands
        cdef list signs
        if isinstance(self, pAdicLazyElement_zero):
            return -other
        if isinstance(other, pAdicLazyElement_zero):
            return self
        if isinstance(self, pAdicLazyElement_add):
            summands = list((<pAdicLazyElement_add>self)._summands)
            signs = list((<pAdicLazyElement_add>self)._signs)
        else:
            summands = [self]
            signs = [True]
        if isinstance(other, pAdicLazyElement_add):
            summands.extend((<pAdicLazyElement_add>other)._summands)
            signs.extend([ not sign for sign in (<pAdicLazyElement_add>other)._signs ])
        else:
            summands.append(other)
            signs.append(False)
        return pAdicLazyElement_add(self._parent, summands, signs)

    cpdef _neg_(self):
        cdef list summands
        cdef list signs
        if isinstance(self, pAdicLazyElement_add):
            summands = list((<pAdicLazyElement_add>self)._summands)
            signs = [ not sign for sign in (<pAdicLazyElement_add>self)._signs ]
        else:
            summands = [self]
            signs = [False]
        return pAdicLazyElement_add(self._parent, summands, signs)

    cpdef _mul_(self, other):
        if isinstance(self, pAdicLazyElement_zero) or isinstance(other, pAdicLazyElement_one):
            return self
        if isinstance(self, pAdicLazyElement_one) or isinstance(other, pAdicLazyElement_zero):
            return other
        return pAdicLazyElement_mul(self._parent, self, <pAdicLazyElement>other)


# Assignations
##############

# Zero

cdef class pAdicLazyElement_zero(pAdicLazyElement):
    def __init__(self, parent):
        pAdicLazyElement.__init__(self, parent)
        self._valuation = maxordp

    cpdef jump(self, prec):
        return True

    cdef bint next(self) except -1:
        return True


# One

cdef class pAdicLazyElement_one(pAdicLazyElement):
    def __init__(self, parent):
        pAdicLazyElement.__init__(self, parent)
        fmpz_poly_set_ui(self._digits, 1)

    cpdef jump(self, prec):
        if self._precrel < prec:
            self._precrel = prec
        return True

    cdef bint next(self) except -1:
        self._precrel += 1
        return True


# Value

cdef class pAdicLazyElement_value(pAdicLazyElement):
    def __init__(self, parent, value, prec=None, maxprec=None):
        pAdicLazyElement.__init__(self, parent)
        cdef Integer x = ZZ(value)
        fmpz_poly_set_mpz(self._digits, x.value)
        self._maxprec = maxprec

    cdef bint next(self) except -1:
        if (self._maxprec is not None) and (self._valuation + self._precrel >= self._maxprec):
            return False
        reduce_coeff(self._digits, self._precrel, self._prime)
        if self._precrel == 0 and fmpz_is_zero(get_coeff(self._digits, 0)):
            self._valuation += 1
            fmpz_poly_shift_right(self._digits, self._digits, 1)
        else:
            self._precrel += 1
        return True


# Random

cdef class pAdicLazyElement_random(pAdicLazyElement):
    def __cinit__(self):
        flint_randinit(self._randstate)

    def __dealloc__(self):
        flint_randclear(self._randstate)

    def __init__(self, parent):
        pAdicLazyElement.__init__(self, parent)

    cdef bint next(self) except -1:
        cdef fmpz_t r;
        fmpz_randm(r, self._randstate, self._prime)
        if self._precrel == 0 and fmpz_is_zero(r):
            self._valuation += 1
        else:
            fmpz_poly_set_coeff_fmpz(self._digits, self._precrel, r);
            self._precrel += 1
        return True


# Operations
############

# Shift

cdef class pAdicLazyElement_shift(pAdicLazyElement):
    def __init__(self, parent, pAdicLazyElement x, slong s):
        pAdicLazyElement.__init__(self, parent)
        self._x = x
        self._shift = s
        self._valuation = x._valuation - s
        self._precrel = x._precrel
        fmpz_poly_set(self._digits, x._digits)

    cdef bint next(self) except -1:
        cdef slong n = self._valuation + self._precrel
        cdef pAdicLazyElement x = self._x
        cdef fmpz* digit

        if not x.jump(n + self._shift + 1):
            return False
        digit = get_coeff(x._digits, n + self._shift - x._valuation)
        fmpz_poly_set_coeff_fmpz(self._digits, self._precrel, digit)
        if self._precrel == 0 and fmpz_is_zero(get_coeff(self._digits, 0)):
            self._valuation += 1
            fmpz_poly_shift_right(self._digits, self._digits, 1)
        else:
            self._precrel += 1
        return True


# Addition

cdef class pAdicLazyElement_add(pAdicLazyElement):
    def __init__(self, parent, summands, signs):
        pAdicLazyElement.__init__(self, parent)
        self._summands = summands
        self._signs = signs

    cdef bint next(self) except -1:
        cdef pAdicLazyElement summand
        cdef slong n = self._valuation + self._precrel
        cdef fmpz* coeff
        for summand in self._summands:
            if not summand.jump(n+1):
                return False
        cdef slong i
        for i in range(len(self._summands)):
            summand = self._summands[i]
            coeff = get_coeff(summand._digits, n - summand._valuation)
            if self._signs[i]:
                iadd_coeff(self._digits, coeff, n - self._valuation)
            else:
                isub_coeff(self._digits, coeff, n - self._valuation)
        reduce_coeff(self._digits, self._precrel, self._prime)
        if self._precrel == 0 and fmpz_is_zero(get_coeff(self._digits, 0)):
            self._valuation += 1
            fmpz_poly_shift_right(self._digits, self._digits, 1)
        else:
            self._precrel += 1
        return True


# Multiplication

cdef class pAdicLazyElement_mul(pAdicLazyElement):
    def __cinit__(self):
        fmpz_poly_init(self._tmp)

    def __dealloc__(self):
        fmpz_poly_clear(self._tmp)

    def __init__(self, parent, x, y):
        pAdicLazyElement.__init__(self, parent)
        self._x = <pAdicLazyElement?>x
        self._y = <pAdicLazyElement?>y
        self._valuation = self._x._valuation + self._y._valuation

    cdef bint next(self) except -1:
        cdef pAdicLazyElement x = self._x
        cdef pAdicLazyElement y = self._y
        cdef slong n = self._valuation + self._precrel

        cdef bint success = (x.jump(n - y._valuation + 1) 
                         and y.jump(n - x._valuation + 1))
        if self._precrel == 0:
            self._valuation = x._valuation + y._valuation
            if self._valuation > n:
                return True
            if self._valuation < n or x._precrel == 0 or y._precrel == 0:
                return False
        elif not success:
            return False

        n = self._precrel
        cdef slong m = n + 2
        cdef slong len = 1
        cdef fmpz_poly_t slicex, slicey
        while m > 1:
            get_slice(slicex, x._digits, len - 1, len)
            get_slice(slicey, y._digits, (m-1)*len - 1, len)
            fmpz_poly_mul(self._tmp, slicex, slicey)
            iadd_shifted(self._digits, self._tmp, n)
            if m > 2:
                get_slice(slicex, x._digits, (m-1)*len - 1, len)
                get_slice(slicey, y._digits, len - 1, len)
                fmpz_poly_mul(self._tmp, slicex, slicey)
                iadd_shifted(self._digits, self._tmp, n)
            if m & 1:
                break
            m >>= 1
            len <<= 1
        reduce_coeff(self._digits, self._precrel, self._prime)
        self._precrel += 1
        return True


# Division

#cdef class pAdicLazyElement_div(pAdicLazyElement):


# Square root

#cdef class pAdicLazyElement_sqrt(pAdicLazyElement):



# Self-referent definitions
###########################

cdef class pAdicLazyElement_selfref(pAdicLazyElement):
    def __init__(self, parent):
        pAdicLazyElement.__init__(self, parent)
        self._definition = None
        self._next = False

    def set(self, definition):
        if self._definition is not None:
            raise ValueError("this self-referent number is already defined")
        self._definition = <pAdicLazyElement?>self._parent(definition)
        try:
            self.next()
        except RecursionError:
            self._definition = None
            self._valuation = 0
            self._precrel = 0
            raise

    def __eq__(self, other):
        if self._definition is None:
            self.set(other)
            return True
        return pAdicLazyElement.__eq__(self, other)

    cdef bint next(self) except -1:
        cdef pAdicLazyElement definition = self._definition
        cdef fmpz* digit

        if definition is None:
            return False
        if self._next:
            raise RecursionError("definition looks circular")

        self._next = True
        success = definition.jump(self._valuation + self._precrel + 1)
        if success:
            if definition._valuation > self._valuation:
                self._valuation = definition._valuation
            else:
                digit = get_coeff(definition._digits, self._precrel)
                fmpz_poly_set_coeff_fmpz(self._digits, self._precrel, digit);
                self._precrel += 1
        self._next = False
        return success


# Expansion
###########

cdef class ExpansionIter(object):
    cdef pAdicLazyElement elt
    cdef slong start
    cdef slong stop
    cdef slong current

    def __init__(self, pAdicLazyElement elt, start, stop):
        self.elt = elt
        if start is None:
            self.start = 0
        else:
            self.start = long(start)
        if stop is None:
            self.stop = self.start - 1
        else: 
            self.stop = long(stop)
            if self.stop <= self.start:
                self.stop = self.start
        self.current = self.start

    def __repr__(self):
        return "%s-adic expansion of %s" % (self.elt._parent.prime(), self.elt)

    def __len__(self):
        if self.stop < self.start:
            raise NotImplementedError("infinite sequence")
        return ZZ(self.stop - self.start)

    def __iter__(self):
        return self

    def __next__(self):
        if self.stop < self.start or self.current < self.stop:
            digit = self.elt._digit(self.current)
            self.current += 1
            return digit
        else:
            raise StopIteration
