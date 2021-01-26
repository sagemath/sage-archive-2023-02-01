# ****************************************************************************
#       Copyright (C) 2021 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

import re

from libc.stdlib cimport malloc, free
from sage.libs.flint.types cimport slong
from sage.libs.gmp.mpz cimport mpz_sizeinbase, mpz_tstbit

from sage.structure.element import coerce_binop
from sage.structure.element cimport have_same_parent
from sage.structure.coerce cimport coercion_model

from sage.rings.all import ZZ
from sage.rings.integer cimport Integer
from sage.rings.infinity import Infinity

from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.padics.padic_lazy_errors cimport *
from sage.rings.padics.padic_lazy_errors import raise_error

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1


cdef cdigit tmp_digit
cdef celement tmp_poly
digit_init(tmp_digit)
element_init(tmp_poly)


cdef class LazyElement(pAdicGenericElement):
    def __init__(self, parent):
        pAdicGenericElement.__init__(self, parent)
        self.prime_pow = <PowComputer_class>self._parent.prime_pow
        self._valuation = 0
        self._precrel = 0
        self._precbound = maxordp

    cdef int _init_jump(self) except -1:
        cdef slong default_prec = self._parent.default_prec()
        cdef int error = 0
        if self._precbound < maxordp:
            error = self._jump_c(self._precbound)
            raise_error(error, True)
        return error

    cpdef bint _is_base_elt(self, p) except -1:
        return True

    def prime(self):
        return self._parent.prime()

    cdef cdigit_ptr _getdigit_relative(self, slong i):
        pass

    cdef cdigit_ptr _getdigit_absolute(self, slong i):
        pass

    cdef void _getslice_relative(self, celement slice, slong start, slong length):
        pass

    cdef int _next_c(self):
        raise NotImplementedError("must be implemented in subclasses")

    cdef int _jump_c(self, slong prec):
        cdef int error
        if prec > maxordp:
            return ERROR_OVERFLOW
        cdef slong pr = min(prec, self._precbound)
        while self._precrel + self._valuation < pr:
            error = self._next_c()
            if error:
                return error
        if prec > self._precbound:
            return ERROR_PRECISION
        return 0

    cdef int _jump_relative_c(self, slong prec, slong halt):
        if self._valuation >= maxordp:
            return 0
        cdef int error = 0
        if self._valuation <= -maxordp:
            error = self._next_c()
        cdef slong valhalt = min(halt, self._precbound)
        while not error and self._valuation < valhalt and self._precrel == 0:
            error = self._next_c()
        if self._valuation >= self._precbound:
            error |= ERROR_PRECISION
        elif self._valuation >= halt:
            error |= ERROR_ABANDON
        if not error:
            error = self._jump_c(self._valuation + prec)
        return error

    def digit(self, slong i):
        cdef int error = self._jump_c(i+1)
        if error:
            raise_error(error)
        cdef cdigit_ptr coeff = self._getdigit_absolute(i)
        return digit_get_sage(coeff)

    def expansion(self, n=None):
        if n is None:
            return ExpansionIter(self, self._valuation, self._valuation + self._precrel)
        if isinstance(n, slice):
            if n.step is not None:
                raise NotImplementedError("step is not allowed")
            if n.start is None:
                start = 0
            else:
                start = Integer(n.start)
            if n.stop is None:
                stop = None
            else:
                stop = Integer(n.stop)
            return ExpansionIter(self, start, stop)
        else:
            n = Integer(n)
            if n >= maxordp:
                raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
            return self.digit(n)

    def slice(self, start=None, stop=None):
        if start is None:
            start = maxordp
        elif start >= maxordp:
            raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
        if stop is None:
            stop = maxordp
        elif stop >= maxordp:
            raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
        return element_class_slice(self._parent, self, start, stop, 0)

    def _repr_(self):
        # This code should be integrated to the p-adic printer
        if self._valuation <= -maxordp:
            return "valuation not known"
        if self._precbound >= maxordp:
            unbounded = True
            try:
                if self.prime_pow.in_field:
                    x = self.at_precision_relative(permissive=False)
                else:
                    x = self.at_precision_absolute(permissive=False)
            except (ValueError, PrecisionError, RecursionError):
                unbounded = False
                x = element_class_bound(self._parent, self, self._valuation + self._precrel)
        else:
            unbounded = False
            x = self
        s = pAdicGenericElement._repr_(x)
        mode = self._parent._printer.dict()['mode']
        if unbounded:
            s = re.sub(r'O\(.*\)$', '...', s)
        elif mode == "digits":
            s = re.sub(r'^\.\.\.\??', '...?', s)
        elif mode == "bars":
            s = re.sub(r'^\.\.\.\|?', '...?|', s)
        if s == "...":
            s = "0 + ..."
        return s

    cdef bint _is_equal(self, LazyElement right, slong prec, bint permissive) except -1:
        cdef int error
        cdef slong i
        if self._valuation <= -maxordp:
            error = self._next_c()
            raise_error(error)
        if right._valuation <= -maxordp:
            error = right._next_c()
            raise_error(error)
        for i in range(min(self._valuation, right._valuation), prec):
            while self._valuation + self._precrel <= i:
                error = self._next_c()
                raise_error(error, permissive)
                if error:   # not enough precision
                    return True
            while right._valuation + right._precrel <= i:
                error = right._next_c()
                raise_error(error, permissive)
                if error:   # not enough precision
                    return True
            if not digit_equal(self._getdigit_absolute(i), right._getdigit_absolute(i)):
                return False
        return True

    @coerce_binop
    def is_equal_at_precision(self, other, prec):
        return self._is_equal(other, min(prec, maxordp), False)

    def __eq__(self, other):
        if not have_same_parent(self, other):
            a, b = coercion_model.canonical_coercion(self, other)
            return a == b
        cdef LazyElement right = <LazyElement?>other
        prec = min(self._valuation + self._precrel, right._valuation + right._precrel)
        c = self._is_equal(right, prec, True)
        if not c:
            return False
        prec = min(self._precbound, right._precbound)
        if prec < maxordp:
            return self._is_equal(right, prec, True)
        else:
            raise PrecisionError("unable to decide equality; try to bound precision")

    cpdef bint _is_exact_zero(self) except -1:
        return self._valuation >= maxordp

    cpdef bint _is_inexact_zero(self) except -1:
        return self._precrel == 0

    def is_exact(self):
        return self._is_exact_zero()

    def precision_absolute(self):
        if self._is_exact_zero():
            return Infinity
        if self._valuation <= -maxordp:
            raise PrecisionError("no lower bound on the valuation is known")
        if self._precbound >= maxordp:
            return Infinity
        return Integer(self._precrel + self._valuation)

    def precision_relative(self):
        if self._valuation <= -maxordp:
            raise PrecisionError("no lower bound on the valuation is known")
        if self._precbound >= maxordp:
            return Infinity
        return Integer(self._precrel)

    def at_precision_absolute(self, prec=None, permissive=None):
        if prec is None:
            if permissive is None:
                permissive = True
            prec = self._parent.default_prec()
        else:
            prec = Integer(prec)
            if not self.prime_pow.in_field and prec < 0:
                raise ValueError("precision must be nonnegative")
            if prec > maxordp:
                raise OverflowError("beyond maximal precision (which is %s)" % maxordp)
        error = self._jump_c(prec)
        if permissive is None:
            permissive = False
        raise_error(error, permissive)
        return element_class_bound((<LazyElement>self)._parent, self, prec)

    def add_bigoh(self, prec):
        return self.at_precision_absolute(prec, True)

    def at_precision_relative(self, prec=None, halt=True, permissive=None):
        default_prec = self._parent.default_prec()
        if prec is None:
            if permissive is None:
                permissive = True
            prec = default_prec
        else:
            prec = Integer(prec)
            if prec < 0:
                raise ValueError("precision must be nonnegative")
            if prec > maxordp:
                raise OverflowError("beyond maximal precision (which is %s)" % maxordp)
        if halt is True:
            halt = self._valuation + default_prec
        elif halt is False:
            halt = maxordp
        else:
            halt = min(maxordp, halt)
        error = self._jump_relative_c(prec, halt)
        if permissive is None:
            permissive = False
        raise_error(error, permissive)
        return element_class_bound((<LazyElement>self)._parent, self, self._valuation +  prec)

    def __matmul__(self, prec):
        if (<LazyElement>self).prime_pow.in_field:
            return self.at_precision_relative(prec)
        else:
            return self.at_precision_absolute(prec)

    def lift_to_precision(self, absprec=None):
        if self._precbound >= maxordp:
            return self
        cdef slong prec
        cdef slong default_prec = self._parent.default_prec()
        if absprec is None:
            if self.prime_pow.in_field:
                self._jump_relative_c(1, self._precbound)
                if self._precrel == 0:
                    return self._parent.zero()
                prec = self._valuation + default_prec
            else:
                prec = default_prec
        else:
            if absprec > maxordp:
                raise OverflowError("beyond the maximal precision (which is %s)" % maxordp)
            prec = absprec
        if prec <= self._precbound:
            return self
        cdef LazyElement ans = element_class_slice(self._parent, self, -maxordp, self._precbound, 0)
        ans._precbound = prec
        ans._init_jump()
        return ans

    cdef long valuation_c(self):
        return self._valuation

    def valuation(self, halt=True):
        if self._is_exact_zero():
            return Infinity
        if self._valuation <= -maxordp:
            raise PrecisionError("no lower bound on the valuation is known")
        if halt is True:
            halt = self._valuation + self._parent.default_prec()
        elif halt is False:
            halt = maxordp
        else:
            halt = min(maxordp, halt)
        self._jump_relative_c(1, halt)
        if self._precbound >= maxordp and self._precrel == 0:
            raise PrecisionError("cannot determine the valuation; try to increase the halting precision")
        return Integer(self._valuation)

    def unit_part(self, halt=True):
        val = self.valuation(halt)
        return self >> val

    def val_unit(self, halt=True):
        val = self.valuation(halt)
        return val, self >> val

    def residue(self, slong prec=1, field=None, check_prec=True):
        if prec >= maxordp:
            raise OverflowError
        if prec < 0:
            raise ValueError("cannot reduce modulo a negative power of p")
        error = self._jump_c(prec)
        if error:
            raise_error(error, not check_prec)
        if self._valuation < 0:
            raise ValueError("element must have non-negative valuation in order to compute residue")
        cdef celement digits
        self._getslice_relative(digits, 0, min(self._precrel, prec - self._valuation))
        cdef Integer ans = element_get_sage(digits, self.prime_pow)
        if field and prec == 1:
            return self._parent.residue_class_field()(ans)
        else:
            return self._parent.residue_ring(prec)(ans)

    def lift(self, prec=None):
        if prec is None:
            if self._precbound < maxordp:
                prec = self._precbound
            else:
                raise ValueError("you must specify a precision for unbounded elements")
        else:
            prec = Integer(prec)
        cdef int error = self._jump_c(prec)
        if error:
            raise_error(error)
        if prec < self._valuation:
            return Integer(0)
        cdef celement digits
        self._getslice_relative(digits, 0, prec - self._valuation)
        cdef Integer ans = element_get_sage(digits, self.prime_pow)
        if self._valuation:
            ans *= self._parent.prime() ** self._valuation
        return ans

    def __rshift__(self, s):
        cdef slong start
        cdef slong shift = long(s)
        if shift:
            if (<LazyElement>self)._parent.is_field():
                start = -maxordp
            else:
                start = shift
            return element_class_slice((<pAdicLazyElement>self)._parent, self, start, maxordp, shift)
        else:
            return self

    def __lshift__(self, s):
        return self.__rshift__(-s)

    cpdef _add_(self, other):
        cdef list summands
        cdef list signs
        if isinstance(self, LazyElement_zero):
            return other
        if isinstance(other, LazyElement_zero):
            return self
        if isinstance(self, LazyElement_add):
            summands = list((<LazyElement_add>self)._summands)
            signs = list((<LazyElement_add>self)._signs)
        else:
            summands = [self]
            signs = [True]
        if isinstance(other, LazyElement_add):
            summands.extend((<LazyElement_add>other)._summands)
            signs.extend((<LazyElement_add>other)._signs)
        else:
            summands.append(other)
            signs.append(True)
        return element_class_add(self._parent, summands, signs)

    cpdef _sub_(self, other):
        cdef list summands
        cdef list signs
        if isinstance(self, LazyElement_zero):
            return -other
        if isinstance(other, LazyElement_zero):
            return self
        if isinstance(self, LazyElement_add):
            summands = list((<LazyElement_add>self)._summands)
            signs = list((<LazyElement_add>self)._signs)
        else:
            summands = [self]
            signs = [True]
        if isinstance(other, LazyElement_add):
            summands.extend((<LazyElement_add>other)._summands)
            signs.extend([ not sign for sign in (<LazyElement_add>other)._signs ])
        else:
            summands.append(other)
            signs.append(False)
        return element_class_add(self._parent, summands, signs)

    cpdef _neg_(self):
        cdef list summands
        cdef list signs
        if isinstance(self, LazyElement_add):
            summands = list((<LazyElement_add>self)._summands)
            signs = [ not sign for sign in (<LazyElement_add>self)._signs ]
        else:
            summands = [self]
            signs = [False]
        return element_class_add(self._parent, summands, signs)

    cpdef _mul_(self, other):
        if isinstance(self, LazyElement_zero) or isinstance(other, LazyElement_one):
            return self
        if isinstance(self, LazyElement_one) or isinstance(other, LazyElement_zero):
            return other
        return element_class_mul(self._parent, self, <LazyElement>other)

    cpdef _div_(self, other):
        if isinstance(other, LazyElement_zero):
            return ZeroDivisionError("cannot divide by zero")
        if isinstance(other, LazyElement_one):
            return self
        return element_class_div(self._parent.fraction_field(), self, <LazyElement>other)

    def __invert__(self):
        if isinstance(self, LazyElement_zero):
            return ZeroDivisionError("cannot divide by zero")
        if isinstance(self, LazyElement_one):
            return self
        return element_class_div(self._parent.fraction_field(), self._parent.one(), self)

    def sqrt(self):
        return element_class_sqrt(self._parent, self)


cdef class LazyElement_abandon(LazyElement):
    def __init__(self):
        self._valuation = -maxordp

    cdef int _next_c(self):
        return ERROR_ABANDON

cdef lazyelement_abandon = LazyElement_abandon()


cdef class LazyElement_init(LazyElement):
    def __cinit__(self):
        element_init(self._digits)

    def __dealloc__(self):
        element_clear(self._digits)

    cdef cdigit_ptr _getdigit_relative(self, slong i):
        return element_get_digit(self._digits, i)

    cdef cdigit_ptr _getdigit_absolute(self, slong i):
        return element_get_digit(self._digits, i - self._valuation)

    cdef void _getslice_relative(self, celement slice, slong start, slong length):
        element_get_slice(slice, self._digits, start, length)



# Assignations
##############

# Zero

cdef class LazyElement_zero(LazyElement):
    def __init__(self, parent):
        LazyElement.__init__(self, parent)
        self._valuation = maxordp

    cdef cdigit_ptr _getdigit_relative(self, slong i):
        return digit_zero

    cdef cdigit_ptr _getdigit_absolute(self, slong i):
        return digit_zero

    cdef void _getslice_relative(self, celement slice, slong start, slong length):
        element_init(slice)

    cdef int _jump_c(self, slong prec):
        return 0

    cdef int _next_c(self):
        return 0


# One

cdef class LazyElement_one(LazyElement_init):
    def __init__(self, parent):
        LazyElement.__init__(self, parent)
        element_set_digit_ui(self._digits, 1, 0)
        self._precrel = self._parent.default_prec()

    cdef int _jump_c(self, slong prec):
        if self._precrel < prec:
            self._precrel = prec
        return 0

    cdef int _next_c(self):
        self._precrel += 1
        return 0


# Bound

cdef class LazyElement_bound(LazyElement):
    def __init__(self, parent, LazyElement x, precbound=None):
        LazyElement.__init__(self, parent)
        self._x = x
        if precbound is None:
            self._precbound = x._precbound
        else:
            self._precbound = min(x._precbound, precbound)
        self._valuation = min(x._valuation, self._precbound)
        self._precrel = min(x._precrel, self._precbound - self._valuation)
        self._init_jump()

    cdef cdigit_ptr _getdigit_relative(self, slong i):
        return self._x._getdigit_relative(i)

    cdef cdigit_ptr _getdigit_absolute(self, slong i):
        return self._x._getdigit_absolute(i)

    cdef void _getslice_relative(self, celement slice, slong start, slong length):
        self._x._getslice_relative(slice, start, length)

    cdef int _jump_c(self, slong prec):
        cdef LazyElement x = self._x
        cdef int error
        if prec > self._precbound:
            error = ERROR_PRECISION | x._jump_c(self._precbound)
        else:
            error = x._jump_c(prec)
        self._precbound = min(self._precbound, x._precbound)
        self._valuation = min(x._valuation, self._precbound)
        self._precrel = min(x._precrel, self._precbound - self._valuation)
        return error

    cdef int _next_c(self):
        cdef LazyElement x = self._x
        cdef int error = x._next_c()
        self._precbound = min(self._precbound, x._precbound)
        self._valuation = min(x._valuation, self._precbound)
        self._precrel = min(x._precrel, self._precbound - self._valuation)
        return error


# Value

cdef class LazyElement_value(LazyElement_init):
    def __init__(self, parent, value, slong shift=0, precbound=None):
        LazyElement.__init__(self, parent)
        element_set_digit_sage(self._digits, value, 0)
        self._valuation = -shift
        self._finished = False
        if precbound is not None and precbound is not Infinity:
            self._precbound = min(maxordp, precbound)
        self._init_jump()

    cdef int _jump_c(self, slong prec):
        if not self._finished:
            return LazyElement._jump_c(self, prec)
        if (self._precbound is not None) and (prec > self._precbound):
            self._precrel = self._precbound - self._valuation
            return ERROR_PRECISION
        cdef slong precrel = min(prec, maxordp) - self._valuation
        if precrel > self._precrel:
            self._precrel = precrel
        if prec < maxordp:
            return 0
        return ERROR_OVERFLOW

    cdef int _next_c(self):
        # The algorithm is not optimal (quadratic in the precision),
        # but it sounds okay
        element_reduce_digit(self._digits, self._precrel, self.prime_pow)
        if self._precrel == 0 and digit_is_zero(self._getdigit_relative(0)):
            self._valuation += 1
            element_shift_right(self._digits)
        else:
            self._precrel += 1
        self._finished = digit_is_zero(self._getdigit_relative(self._precrel))
        return 0


# Random

cdef class LazyElement_random(LazyElement_init):
    def __init__(self, parent, integral, precbound):
        LazyElement.__init__(self, parent)
        if not integral:
            self._valuation = ZZ.random_element()
        if precbound is not None:
            self._precbound = min(maxordp, precbound)
        self._init_jump()

    cdef int _next_c(self):
        cdef cdigit r
        digit_random(r, self.prime_pow)
        if self._precrel == 0 and digit_is_zero(r):
            self._valuation += 1
            digit_clear(r)
        else:
            element_set_digit(self._digits, r, self._precrel)
            self._precrel += 1
        return 0


# Operations
############

# Shift

cdef class LazyElement_slice(LazyElement):
    def __init__(self, parent, LazyElement x, slong start, slong stop, slong shift):
        # self[0] = x[shift]
        LazyElement.__init__(self, parent)
        self._x = x
        self._start = start
        self._shift = shift
        self._valuation = max(x._valuation, start) - shift
        self._precrel = min(x._precrel + x._valuation, stop) - self._valuation - shift
        if x._precbound < stop:
            self._precbound = min(maxordp, x._precbound - shift)
        self._stop = stop
        if self._precrel < 0:
            self._precrel = 0
        while self._precrel > 0 and digit_is_zero(self._getdigit_relative(0)):
            self._precrel -= 1
            self._valuation += 1
        self._init_jump()

    cdef cdigit_ptr _getdigit_relative(self, slong i):
        return self._getdigit_absolute(i + self._valuation)

    cdef cdigit_ptr _getdigit_absolute(self, slong i):
        cdef slong j = i + self._shift
        if j < self._start or j >= self._stop:
            return digit_zero
        else:
            return self._x._getdigit_absolute(j)

    cdef void _getslice_relative(self, celement slice, slong start, slong length):
        cdef LazyElement x = self._x
        cdef slong s = start + self._valuation + self._shift
        cdef slong start_absolute = max(self._start, s)
        cdef slong stop_absolute = min(self._stop, s + length)
        x._getslice_relative(slice, start_absolute - x._valuation, stop_absolute - start_absolute)

    cdef int _jump_c(self, slong prec):
        cdef LazyElement x = self._x
        cdef int error = 0
        cdef slong pr
        if prec <= self._valuation + self._precrel:
            return 0
        if prec > self._precbound:
            prec = self._precbound
            error = ERROR_PRECISION
        cdef int errorx = x._jump_c(min(prec + self._shift, self._stop))
        pr = max(self._valuation, x._valuation + x._precrel - self._shift)
        if self._precrel == 0:
            while self._valuation < pr and digit_is_zero(self._getdigit_relative(0)):
                self._valuation += 1
        self._precrel = prec - self._valuation
        if errorx:
            return errorx
        else:
            return error

    cdef int _next_c(self):
        cdef int error
        cdef slong n = self._precrel + self._valuation
        n += self._shift
        if n <= self._stop:
            error = self._x._jump_c(n+1)
            if error:
                return error
        if self._precrel == 0 and (n > self._stop or digit_is_zero(self._getdigit_relative(0))):
            self._valuation += 1
        else:
            self._precrel += 1
        return 0


# Addition

cdef class LazyElement_add(LazyElement_init):
    def __init__(self, parent, summands, signs):
        LazyElement.__init__(self, parent)
        self._summands = summands
        self._signs = signs
        self._valuation = min((<LazyElement>summand)._valuation for summand in summands)
        self._precbound = min((<LazyElement>summand)._precbound for summand in summands)
        self._init_jump()

    cdef int _next_c(self):
        cdef LazyElement summand
        cdef slong n = self._valuation + self._precrel
        cdef cdigit_ptr coeff
        cdef int error

        for summand in self._summands:
            error = summand._jump_c(n+1)
            if error:
                return error
        cdef slong i
        for i in range(len(self._summands)):
            summand = self._summands[i]
            coeff = summand._getdigit_absolute(n)
            if self._signs[i]:
                element_iadd_digit(self._digits, coeff, n - self._valuation)
            else:
                element_isub_digit(self._digits, coeff, n - self._valuation)
        element_reduce_digit(self._digits, self._precrel, self.prime_pow)
        if self._precrel == 0 and digit_is_zero(self._getdigit_relative(0)):
            self._valuation += 1
            element_shift_right(self._digits)
        else:
            self._precrel += 1
        return 0



# Multiplication

cdef class LazyElement_mul(LazyElement_init):
    def __cinit__(self):
        digit_init(self._lastdigit_x)
        digit_init(self._lastdigit_y)

    def __dealloc__(self):
        digit_clear(self._lastdigit_x)
        digit_clear(self._lastdigit_y)

    def __init__(self, parent, LazyElement x, LazyElement y):
        LazyElement.__init__(self, parent)
        self._x = x
        self._y = y
        self._valuation = x._valuation + y._valuation
        if x._precbound < maxordp:
            self._precbound = min(self._precbound, y._valuation + x._precbound)
        if y._precbound < maxordp:
            self._precbound = min(self._precbound, x._valuation + y._precbound)
        self._init_jump()

    cdef int _next_c(self):
        global tmp_digit, tmp_poly
        cdef LazyElement x = self._x
        cdef LazyElement y = self._y
        cdef slong n = self._valuation + self._precrel

        cdef int errorx = x._jump_c(n - y._valuation + 1)
        cdef int errory = y._jump_c(n - x._valuation + 1)
        cdef int error = errorx | errory
        if self._precrel == 0:
            self._valuation = x._valuation + y._valuation
            if self._valuation > n:
                return 0
            if self._valuation < n or x._precrel == 0 or y._precrel == 0:
                return error | ERROR_PRECISION
        elif error:
            return error

        n = self._precrel
        digit_set(self._lastdigit_x, x._getdigit_relative(n))
        digit_set(self._lastdigit_y, y._getdigit_relative(n))
        digit_mul(tmp_digit, x._getdigit_relative(0), self._lastdigit_y)
        element_iadd_digit(self._digits, tmp_digit, n)
        if n:
            digit_mul(tmp_digit, self._lastdigit_x, y._getdigit_relative(0))
            element_iadd_digit(self._digits, tmp_digit, n)

        cdef slong m = n + 2
        cdef slong len = 1
        cdef celement slicex, slicey
        while (m & 1 == 0) and (m > 3):
            m >>= 1
            len <<= 1
            x._getslice_relative(slicex, len - 1, len)
            y._getslice_relative(slicey, (m-1)*len - 1, len)
            element_mul(tmp_poly, slicex, slicey)
            element_iadd_slice(self._digits, tmp_poly, n)
            if m > 2:
                x._getslice_relative(slicex, (m-1)*len - 1, len)
                y._getslice_relative(slicey, len - 1, len)
                element_mul(tmp_poly, slicex, slicey)
                element_iadd_slice(self._digits, tmp_poly, n)

        element_reduce_digit(self._digits, n, self.prime_pow)
        self._precrel += 1
        return 0

    cdef int _update_last_digit(self):
        if self._precrel == 0:
            return ERROR_UNEXPECTED
        cdef LazyElement x = self._x
        cdef LazyElement y = self._y
        cdef slong n = self._precrel - 1
        cdef slong m = n + 2
        cdef slong len = 2
        cdef celement slice
        while (m & 1 == 0) and (m > 3):
            m >>= 1
            len <<= 1
        len -= 1

        digit_sub(tmp_digit, x._getdigit_relative(n), self._lastdigit_x)
        y._getslice_relative(slice, 0, len)
        element_scalarmul(tmp_poly, slice, tmp_digit)
        element_iadd_slice(self._digits, tmp_poly, n)
        if m == 2:
            len -= 1
        digit_sub(tmp_digit, y._getdigit_relative(n), self._lastdigit_y)
        x._getslice_relative(slice, 0, len)
        element_scalarmul(tmp_poly, slice, tmp_digit)
        element_iadd_slice(self._digits, tmp_poly, n)
        if m == 2:
            digit_mul(tmp_digit, tmp_digit, self._lastdigit_x)
            element_iadd_digit(self._digits, tmp_digit, 2*len)
        element_reduce_digit(self._digits, n, self.prime_pow)

        digit_set(self._lastdigit_x, x._getdigit_relative(n))
        digit_set(self._lastdigit_y, y._getdigit_relative(n))
        return 0


cdef class LazyElement_muldigit(LazyElement_init):
    # compute x*y where x is assumed to be in (0,p)
    def __init__(self, parent, LazyElement_div x, LazyElement y):
        LazyElement.__init__(self, parent)
        self._x = <cdigit_ptr>x._inverse
        self._y = y
        self._valuation = y._valuation
        self._init_jump()

    cdef int _next_c(self):
        cdef slong n = self._valuation + self._precrel
        cdef int error = self._y._jump_c(n+1)
        if error:
            return error
        digit_mul(tmp_digit, self._x, self._y._getdigit_absolute(n))
        element_iadd_digit(self._digits, tmp_digit, self._precrel)
        element_reduce_digit(self._digits, self._precrel, self.prime_pow)
        self._precrel += 1
        return 0

    cdef void _erase_first_digit(self):
        # This is a awful hack,
        # but it's very useful for division
        if self._precrel:
            self._valuation += 1
            self._precrel -= 1
            element_shift_right(self._digits)


# Division

cdef class LazyElement_div(LazyElement_init):
    def __cinit__(self):
        digit_init(self._inverse)

    def __dealloc__(self):
        digit_clear(self._inverse)

    def __init__(self, parent, LazyElement num, LazyElement denom, precbound=None):
        LazyElement.__init__(self, parent)
        self._num = num
        self._denom = denom
        if denom._valuation <= -maxordp:
            self._maxprec = maxordp + 1
        else:
            self._maxprec = denom._valuation + max(1, self._parent.default_prec())
        if precbound is not None:
            self._precbound = min(maxordp, precbound)
        cdef int error = self._bootstrap_c()
        if num._precbound < maxordp:
            self._precbound = min(self._precbound, num._precbound - denom._valuation)
        if denom._precbound < maxordp:
            self._precbound = min(self._precbound, denom._precbound + num._valuation - 2*denom._valuation)
        if error:
            self._valuation = -maxordp
        else:
            self._init_jump()

    cdef int _bootstrap_c(self):
        cdef int error
        cdef LazyElement num = self._num
        cdef LazyElement denom = self._denom

        while denom._valuation < self._maxprec and denom._precrel == 0:
            error = denom._next_c()
            if error:
                if error & ERROR_PRECISION:
                    error |= ERROR_DIVISION
                return error
            if self._maxprec > maxordp and denom._valuation > -maxordp:
                self._maxprec = denom._valuation + max(1, self._parent.default_prec())
        if denom._precrel == 0:
            return ERROR_ABANDON

        self._valuation = num._valuation - denom._valuation
        digit_inv(self._inverse, denom._getdigit_relative(0), self.prime_pow)
        self._definition = lazyelement_abandon
        cdef parent = self._parent
        cdef LazyElement a = element_class_muldigit(parent, self, num)
        cdef LazyElement b = element_class_muldigit(parent, self, denom)
        cdef LazyElement c = element_class_slice(parent, b, denom._valuation + 1, maxordp, 0)
        cdef LazyElement d = element_class_mul(parent, c, self)
        self._definition = element_class_add(parent, [a,d], [True,False])
        return 0

    cdef int _next_c(self):
        cdef LazyElement definition = self._definition
        cdef slong val
        if definition is None:
            error = self._bootstrap_c()
            if error:
                return error
        else:
            val = self._valuation + self._denom._valuation
            error = definition._jump_c(val + self._precrel + 1)
            if error:
                return error
            if definition._valuation > val:
                self._valuation = definition._valuation - self._denom._valuation
                self._precbound = definition._precbound - self._denom._valuation
            else:
                digit = definition._getdigit_relative(self._precrel)
                element_set_digit(self._digits, digit, self._precrel)
                self._precrel += 1
        return 0


# Square root

cdef class LazyElement_sqrt(LazyElement_init):
    def __init__(self, parent, LazyElement x):
        LazyElement.__init__(self, parent)
        self._x = x
        if x._valuation <= -maxordp:
            self._valuation = -maxordp
            self._maxprec = maxordp + 1
        else:
            self._valuation = x._valuation >> 1
            self._maxprec = x._valuation + 2*self._parent.default_prec()
        cdef int error = self._bootstrap_c()
        if error & ERROR_NOTSQUARE:
            raise ValueError("not a square")
        if not error:
            if x._precbound < maxordp:
                self._precbound = x._precbound - x._valuation / 2
                if self._parent.prime() == 2:
                    self._precbound -= 1
                self._precbound = min(maxordp, self._precbound)
            self._init_jump()

    cdef int _bootstrap_c(self):
        cdef LazyElement x = self._x
        while x._valuation < self._maxprec and x._precrel == 0:
            error = x._next_c()
            if error:
                return error
            if self._maxprec > maxordp and x._valuation > -maxordp:
                self._maxprec = x._valuation + 2*self._parent.default_prec()
        if x._valuation > -maxordp:
            if x._valuation & 1:
                self._valuation = (x._valuation + 1) >> 1
            else:
                self._valuation = x._valuation >> 1
        if x._precrel == 0:
            return ERROR_ABANDON
        if x._valuation & 1 != 0:
            return ERROR_NOTSQUARE

        cdef parent = self._parent
        cdef slong val = self._valuation
        cdef cdigit digit
        cdef Integer zd, p = self.prime_pow.prime
        cdef LazyElement u, y, c, d

        if p == 2:
            element_set_digit_ui(self._digits, 1, 0)
            self._precrel = 1
            if x._precrel == 1:
                error = x._next_c()
                if error:
                    return error
            if not digit_equal_ui(x._getdigit_relative(1), 0):
                return ERROR_NOTSQUARE
            if x._precrel == 2:
                error = x._next_c()
                if error:
                    return error
            if not digit_equal_ui(x._getdigit_relative(2), 0):
                return ERROR_NOTSQUARE
            zd = Integer(1)
            self._definition = lazyelement_abandon
            u = element_class_slice(parent, self, val + 2, maxordp, val)
            y = element_class_slice(parent, x, -maxordp, maxordp, val + 1)
            c = element_class_value(parent, zd, shift=-val+1) 
            d = element_class_slice(parent, u*u, -maxordp, maxordp, -val + 1)
            self._definition = y + c - d
        else:
            digit_init(digit)
            if digit_sqrt(digit, x._getdigit_relative(0), self.prime_pow):
                digit_clear(digit)
                return ERROR_NOTSQUARE
            element_set_digit(self._digits, digit, 0)
            self._precrel = 1
            zd = digit_get_sage(digit)
            self._definition = lazyelement_abandon
            u = element_class_slice(parent, self, val + 1, maxordp, val)
            y = element_class_slice(parent, x, -maxordp, maxordp, 2*val)
            c = element_class_value(parent, zd*zd)
            d = element_class_value(parent, 2*zd, shift=val)
            self._definition = (y + c - u*u) / d
        return 0

    cdef int _next_c(self):
        cdef LazyElement x = self._x
        cdef LazyElement definition = self._definition
        cdef slong n = self._valuation + self._precrel
        cdef int error
        if definition is None:
            error = self._bootstrap_c()
            if error:
                return error
        else:
            error = definition._jump_c(n+1)
            if error:
                return error
            element_set_digit(self._digits, definition._getdigit_relative(self._precrel), self._precrel)
            self._precrel += 1
        return 0

    def next(self):
        return self._next_c()


# Teichmüller lifts

cdef class LazyElement_teichmuller(LazyElement_init):
    def __init__(self, parent, Integer xbar):
        LazyElement.__init__(self, parent)
        cdef cdigit digit
        digit_init(digit)
        digit_set_sage(digit, xbar)
        digit_mod(digit, digit, self.prime_pow)
        if digit_equal_ui(digit, 0):
            digit_clear(digit)
            self._trivial = True
            self._valuation = maxordp
            return
        element_set_digit(self._digits, digit, 0)
        self._trivial = digit_equal_ui(digit, 1)
        self._precrel += 1

        cdef LazyElement xn = self
        cdef Integer p = self.prime_pow.prime
        cdef int size = mpz_sizeinbase(p.value, 2)
        cdef int i = size - 2
        self._xns = [ ]
        while i >= 0:
            xn = xn * xn
            self._xns.append(xn)
            if mpz_tstbit(p.value, i):
                xn = xn * self
                self._xns.append(xn)
            i -= 1
        self._xp = xn
        self._ready = True
        self._init_jump()

    cdef int _jump_c(self, slong prec):
        if not self._ready:
            return ERROR_ABANDON
        if self._trivial:
            if self._valuation == 0 and self._precrel < prec:
                self._precrel = prec
            return 0
        return LazyElement._jump_c(<LazyElement>self, prec)

    cdef int _next_c(self):
        if self._trivial:
            if self._valuation:
                self._precrel += 1
            return 0
        cdef int error
        cdef LazyElement xp = self._xp
        cdef LazyElement_mul xn
        self._precrel += 1
        xp._next_c()
        element_set_digit(self._digits, xp._getdigit_relative(self._precrel - 1), self._precrel - 1)
        for xn in self._xns:
            error = xn._update_last_digit()
            if error:
                return error | ERROR_UNEXPECTED
        return 0


# Self-referent definitions
###########################

cdef class LazyElement_selfref(LazyElement_init):
    def __init__(self, parent, slong valuation, digits=None):
        LazyElement.__init__(self, parent)
        if valuation >= maxordp:
            raise OverflowError("valuation is too large (maximum is %s)" % maxordp)
        self._valuation = valuation
        self._definition = None
        self._next = maxordp
        cdef cdigit digit
        if digits is not None:
            digits = [ Integer(d) for d in digits ]
            for d in digits:
                digit_set_sage(digit, d)
                element_iadd_digit(self._digits, digit, self._precrel)
                element_reduce_digit(self._digits, self._precrel, self.prime_pow)
                if self._precrel == 0 and digit_is_zero(element_get_digit(self._digits, 0)):
                    self._valuation += 1
                    element_shift_right(self._digits)
                else:
                    self._precrel += 1

    cpdef set(self, LazyElement definition):
        if self._definition is not None:
            raise ValueError("this self-referent number is already defined")
        self._definition = definition
        self._precbound = max(self._valuation + self._precrel, definition._precbound)
        self._init_jump()

    cdef int _next_c(self):
        cdef LazyElement definition = self._definition
        cdef cdigit_ptr digit
        cdef slong n = self._valuation + self._precrel
        cdef slong diffval
        cdef int error

        if definition is None:
            return ERROR_NOTDEFINED
        if n >= self._next:
            return ERROR_CIRCULAR

        cdef slong svenext = self._next
        self._next = n
        error = definition._jump_c(n+1)
        if not error:
            digit = definition._getdigit_absolute(n)
            if self._precrel == 0 and digit_is_zero(digit):
                self._valuation += 1
            else:
                element_set_digit(self._digits, digit, self._precrel)
                self._precrel += 1
        self._next = svenext
        return error


# Expansion
###########

cdef class ExpansionIter(object):
    cdef LazyElement elt
    cdef slong start
    cdef slong stop
    cdef slong current

    def __init__(self, LazyElement elt, start, stop):
        self.elt = elt
        if start is None:
            self.start = 0
        else:
            if start >= maxordp:
                raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
            self.start = long(start)
        if stop is None:
            self.stop = self.start - 1
        else: 
            if stop >= maxordp:
                raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
            self.stop = long(stop)
            if self.stop <= self.start:
                self.stop = self.start
        self.current = self.start

    def __repr__(self):
        return "%s-adic expansion of %s" % (self.elt._parent.prime(), self.elt)

    def __len__(self):
        if self.stop < self.start:
            raise NotImplementedError("infinite sequence")
        return Integer(self.stop - self.start)

    def __iter__(self):
        return self

    def __next__(self):
        if self.stop < self.start or self.current < self.stop:
            digit = self.elt.digit(self.current)
            self.current += 1
            return digit
        else:
            raise StopIteration
