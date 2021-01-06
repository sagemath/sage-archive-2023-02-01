# ****************************************************************************
#       Copyright (C) 2021 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.misc import walltime

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
from sage.structure.element cimport have_same_parent
from sage.structure.coerce cimport coercion_model

from sage.rings.all import ZZ
from sage.rings.integer cimport Integer
from sage.rings.infinity import Infinity
from sage.rings.finite_rings.finite_field_constructor import GF

from sage.rings.padics.pow_computer_flint cimport PowComputer_flint
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.padics.padic_lazy_errors cimport *
from sage.rings.padics.padic_lazy_errors import raise_error, error_to_str

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) - 1
MAXORDP = ZZ(maxordp)


cpdef lazy_sum(parent, summands):
    if not summands:
        return parent.zero()
    n = len(summands)
    signs = n * [True]
    return pAdicLazyElement_add(parent, summands, signs)



cdef class pAdicLazyElement(pAdicGenericElement):
    def __cinit__(self):
        fmpz_poly_init(self._digits)

    def __dealloc__(self):
        fmpz_poly_clear(self._digits)

    def __init__(self, parent):
        pAdicGenericElement.__init__(self, parent)
        self.prime_pow = <PowComputer_flint>self._parent.prime_pow
        self._valuation = 0
        self._precrel = 0

    cpdef bint _is_base_elt(self, p) except -1:
        return True

    def prime(self):
        return self._parent.prime()

    cdef int _next_c(self):
        raise NotImplementedError("must be implemented in subclasses")

    cdef int _jump_c(self, slong prec):
        cdef int error
        cdef slong i, stop
        prec = min(prec, maxordp)
        while self._precrel + self._valuation < prec:
            error = self._next_c()
            if error:
                return error
        if prec < maxordp:
            return 0
        return ERROR_OVERFLOW

    def jump(self, prec=None, quiet=False):
        if prec is None:
            permissive = True
            default_prec = self._parent.default_prec()
            error = 0
            if self._valuation <= -maxordp:
                error = self._next_c()
            if not error:
                if self.prime_pow.in_field:
                    prec = self._valuation + default_prec
                    error = self._jump_c(prec)
                    if not error and self._valuation < prec:
                        error = self._jump_c(self._valuation + default_prec)
                else:
                    error = self._jump_c(default_prec)
        elif prec not in ZZ:
            raise ValueError("precision must be an integer")
        else:
            permissive = False
            prec = ZZ(prec)
            if prec >= maxordp:
                raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
            if prec < -maxordp:
                prec = -maxordp
            error = self._jump_c(prec)
        if quiet:
            return error
        else:
            raise_error(error, permissive)

    def expansion(self, n=None):
        if n is None:
            return ExpansionIter(self, self._valuation, self._valuation + self._precrel)
        if isinstance(n, slice):
            if n.step is not None:
                raise NotImplementedError("step is not allowed")
            if n.start is None:
                start = 0
            elif n.start in ZZ:
                start = ZZ(n.start)
            else:
                raise ValueError("invalid slice")
            if n.stop is None:
                stop = None
            elif n.stop in ZZ:
                stop = ZZ(n.stop)
            else:
                raise ValueError("invalid slice")
            return ExpansionIter(self, start, stop)
        else:
            if n not in ZZ:
                raise IndexError("index must be an integer")
            n = ZZ(n)
            if n >= maxordp:
                raise OverflowError("beyond maximum precision (which is %s)" % maxordp)
            return self._digit(long(n))

    cdef Integer _digit(self, slong i):
        self.jump(i+1)
        cdef Integer ans = PY_NEW(Integer)
        cdef fmpz* coeff = get_coeff(self._digits, i - self._valuation)
        fmpz_get_mpz(ans.value, coeff)
        return ans

    def _repr_(self):
        error = self.jump(quiet=True)
        s = error_to_str(error, permissive=True)
        if s is not None:
            return s
        if self._valuation <= -maxordp:
            return "valuation not known"
        return pAdicGenericElement._repr_(self)

    cdef bint _is_equal(self, pAdicLazyElement right, slong prec) except -1:
        if self._valuation <= -maxordp:
            error = self._next_c()
            if error:
                self._error(error)
        if right._valuation <= -maxordp:
            error = right._next_c()
            if error:
                self._error(error)
        cdef slong i = 0
        cdef slong j = self._valuation - right._valuation
        if j > 0:
            i = -j
            j = 0
        prec -= self._valuation
        while i < prec:
            while self._precrel <= i:
                error = self._next_c()
                if error:
                    self._error(error)
            while right._precrel <= j:
                error = right._next_c()
                if error:
                    self._error(error)
            if not fmpz_equal(get_coeff(self._digits, i), get_coeff(right._digits, j)):
                return False
            i += 1
            j += 1
        return True

    @coerce_binop
    def is_equal_at_precision(self, other, prec):
        return self._is_equal(other, min(prec, maxordp))

    def __eq__(self, other):
        if not have_same_parent(self, other):
            a, b = coercion_model.canonical_coercion(self, other)
            return a == b
        cdef pAdicLazyElement right = <pAdicLazyElement?>other
        if self._valuation <= -maxordp:
            error = self._next_c()           
            if error:
                self._error(error)
        if right._valuation <= -maxordp:
            error = right._next_c()           
            if error:
                right._error(error)
        minprec = min(self.precision_absolute(), right.precision_absolute())
        default_prec = self._parent.default_prec()
        if self.prime_pow.in_field:
            prec = self._valuation + default_prec
            if not self._is_equal(right, max(minprec, prec)):
                return False
            if self._valuation < prec:
                return self._is_equal(right, max(minprec, self._valuation + default_prec))
        else:
            return self._is_equal(right, max(minprec, default_prec))

    cpdef bint _is_exact_zero(self) except -1:
        return isinstance(self, pAdicLazyElement_zero)

    cpdef bint _is_inexact_zero(self) except -1:
        return self._precrel == 0

    def is_exact(self):
        return self._is_exact_zero()

    def precision_absolute(self):
        if self._is_exact_zero():
            return Infinity
        if self._valuation <= -maxordp:
            raise PrecisionError("no lower bound on the valuation is known")
        return ZZ(self._precrel + self._valuation)

    def precision_relative(self):
        if self._valuation <= -maxordp:
            raise PrecisionError("no lower bound on the valuation is known")
        return ZZ(self._precrel)

    def valuation(self, secure=False):
        if self._is_exact_zero():
            return Infinity
        if self._valuation <= -maxordp:
            raise PrecisionError("no lower bound on the valuation is known")
        if secure and self._precrel == 0:
            raise PrecisionError("cannot determine the valuation; try to increase the precision")
        return ZZ(self._valuation)

    def unit_part(self):
        if self._precrel == 0:
            raise PrecisionError("cannot determine the valuation; try to increase the precision")
        return self >> self._valuation

    def lift(self):
        if self._precrel == 0:
            return ZZ(0)
        cdef fmpz_t fans
        fmpz_init(fans)
        cdef fmpz_poly_t digits
        get_slice(digits, self._digits, 0, self._precrel)
        fmpz_poly_evaluate_fmpz(fans, digits, self.prime_pow.fprime)
        cdef Integer ans = PY_NEW(Integer)
        fmpz_get_mpz(ans.value, fans)
        fmpz_clear(fans)
        if self._valuation:
            ans *= self._parent.prime() ** self._valuation
        return ans

    def __rshift__(self, s):
        cdef slong shift = long(s)
        if shift:
            return pAdicLazyElement_shift((<pAdicLazyElement>self)._parent, self, shift, not (<pAdicLazyElement>self)._parent.is_field())
        else:
            return self

    def __lshift__(self, s):
        return self.__rshift__(-s)

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

    cpdef _div_(self, other):
        if isinstance(other, pAdicLazyElement_zero):
            return ZeroDivisionError("cannot divide by zero")
        if isinstance(other, pAdicLazyElement_one):
            return self
        return pAdicLazyElement_div(self._parent.fraction_field(), self, <pAdicLazyElement>other)

    def __invert__(self):
        if isinstance(self, pAdicLazyElement_zero):
            return ZeroDivisionError("cannot divide by zero")
        if isinstance(self, pAdicLazyElement_one):
            return self
        return pAdicLazyElement_div(self._parent.fraction_field(), self._parent.one(), self)

    def sqrt(self):
        return pAdicLazyElement_sqrt(self._parent, self)


# Assignations
##############

# Zero

cdef class pAdicLazyElement_zero(pAdicLazyElement):
    def __init__(self, parent):
        pAdicLazyElement.__init__(self, parent)
        self._valuation = maxordp

    cdef int _jump_c(self, slong prec):
        return 0

    cdef int _next_c(self):
        return 0


# One

cdef class pAdicLazyElement_one(pAdicLazyElement):
    def __init__(self, parent):
        pAdicLazyElement.__init__(self, parent)
        fmpz_poly_set_ui(self._digits, 1)

    cdef int _jump_c(self, slong prec):
        if self._precrel < prec:
            self._precrel = prec
        return 0

    cdef int _next_c(self):
        self._precrel += 1
        return 0


# Value

cdef class pAdicLazyElement_value(pAdicLazyElement):
    def __init__(self, parent, Integer value, slong shift=0, maxprec=None):
        pAdicLazyElement.__init__(self, parent)
        cdef Integer x = value
        fmpz_poly_set_mpz(self._digits, x.value)
        if maxprec is None:
            self._maxprec = maxordp
        else:
            self._maxprec = maxprec
        self._shift = self._valuation = shift
        self._finished = False

    cdef int _jump_c(self, slong prec):
        if not self._finished:
            return pAdicLazyElement._jump_c(self, prec)
        if (self._maxprec is not None) and (prec > self._maxprec):
            self._precrel = self._maxprec - self._valuation
            return ERROR_PRECISION
        cdef slong precrel = min(prec, maxordp) - self._valuation
        if precrel > self._precrel:
            self._precrel = precrel
        if prec < maxordp:
            return 0
        return ERROR_OVERFLOW

    cdef int _next_c(self):
        if (self._maxprec is not None) and (self._valuation + self._precrel >= self._maxprec):
            return ERROR_PRECISION
        reduce_coeff(self._digits, self._precrel, self.prime_pow.fprime)
        if self._precrel == 0 and fmpz_is_zero(get_coeff(self._digits, 0)):
            self._valuation += 1
            fmpz_poly_shift_right(self._digits, self._digits, 1)
        else:
            self._precrel += 1
        self._finished = fmpz_is_zero(get_coeff(self._digits, self._precrel))
        return 0


# Random

cdef flint_rand_t flint_randstate
flint_randinit(flint_randstate)

cdef class pAdicLazyElement_random(pAdicLazyElement):
    def __init__(self, parent, integral):
        pAdicLazyElement.__init__(self, parent)
        if not integral:
            self._valuation = ZZ.random_element()

    cdef int _next_c(self):
        cdef fmpz_t r
        fmpz_randm(r, flint_randstate, self.prime_pow.fprime)
        if self._precrel == 0 and fmpz_is_zero(r):
            self._valuation += 1
        else:
            fmpz_poly_set_coeff_fmpz(self._digits, self._precrel, r);
            self._precrel += 1
        return 0


# Operations
############

# Shift

cdef class pAdicLazyElement_shift(pAdicLazyElement):
    def __init__(self, parent, pAdicLazyElement x, slong s, bint truncate):
        cdef slong start = 0
        cdef fmpz_poly_t digits;
        pAdicLazyElement.__init__(self, parent)
        self._x = x
        self._shift = s
        if truncate and x._valuation < s:
            start = s - x._valuation
            self._valuation = 0
            self._precrel = x._precrel - start
        else:
            self._valuation = x._valuation - s
            self._precrel = x._precrel
        if self._precrel < 0:
            self._precrel = 0
        else:
            get_slice(digits, x._digits, start, self._precrel)
            fmpz_poly_set(self._digits, digits)

    cdef int _next_c(self):
        cdef slong n = self._valuation + self._precrel
        cdef pAdicLazyElement x = self._x
        cdef fmpz* digit
        cdef int error

        error = x._jump_c(n + self._shift + 1)
        if error:
            return error
        digit = get_coeff(x._digits, n + self._shift - x._valuation)
        fmpz_poly_set_coeff_fmpz(self._digits, self._precrel, digit)
        if self._precrel == 0 and fmpz_is_zero(get_coeff(self._digits, 0)):
            self._valuation += 1
            fmpz_poly_shift_right(self._digits, self._digits, 1)
        else:
            self._precrel += 1
        return 0


# Addition

cdef class pAdicLazyElement_add(pAdicLazyElement):
    def __init__(self, parent, summands, signs):
        pAdicLazyElement.__init__(self, parent)
        self._valuation = min((<pAdicLazyElement>summand)._valuation for summand in summands)
        self._summands = summands
        self._signs = signs

    cdef int _next_c(self):
        cdef pAdicLazyElement summand
        cdef slong n = self._valuation + self._precrel
        cdef fmpz* coeff
        cdef int error

        for summand in self._summands:
            error = summand._jump_c(n+1)
            if error:
                return error
        cdef slong i
        for i in range(len(self._summands)):
            summand = self._summands[i]
            coeff = get_coeff(summand._digits, n - summand._valuation)
            if self._signs[i]:
                iadd_coeff(self._digits, coeff, n - self._valuation)
            else:
                isub_coeff(self._digits, coeff, n - self._valuation)
        reduce_coeff(self._digits, self._precrel, self.prime_pow.fprime)
        if self._precrel == 0 and fmpz_is_zero(get_coeff(self._digits, 0)):
            self._valuation += 1
            fmpz_poly_shift_right(self._digits, self._digits, 1)
        else:
            self._precrel += 1
        return 0



# Multiplication

cdef fmpz_t tmp_coeff
cdef fmpz_poly_t tmp_poly
fmpz_init(tmp_coeff)
fmpz_poly_init(tmp_poly)

cdef class pAdicLazyElement_mul(pAdicLazyElement):
    def __init__(self, parent, x, y):
        pAdicLazyElement.__init__(self, parent)
        self._x = <pAdicLazyElement?>x
        self._y = <pAdicLazyElement?>y
        self._valuation = self._x._valuation + self._y._valuation

    cdef int _next_c(self):
        global tmp_coeff, tmp_poly
        cdef pAdicLazyElement x = self._x
        cdef pAdicLazyElement y = self._y
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
        fmpz_mul(tmp_coeff, get_coeff(x._digits, 0), get_coeff(y._digits, n))
        iadd_coeff(self._digits, tmp_coeff, n)
        if n:
            fmpz_mul(tmp_coeff, get_coeff(x._digits, n), get_coeff(y._digits, 0))
            iadd_coeff(self._digits, tmp_coeff, n)

        cdef slong m = n + 2
        cdef slong len = 1
        cdef fmpz_poly_t slicex, slicey
        while (m & 1 == 0) and (m > 3):
            m >>= 1
            len <<= 1
            get_slice(slicex, x._digits, len - 1, len)
            get_slice(slicey, y._digits, (m-1)*len - 1, len)
            fmpz_poly_mul(tmp_poly, slicex, slicey)
            iadd_shifted(self._digits, tmp_poly, n)
            if m > 2:
                get_slice(slicex, x._digits, (m-1)*len - 1, len)
                get_slice(slicey, y._digits, len - 1, len)
                fmpz_poly_mul(tmp_poly, slicex, slicey)
                iadd_shifted(self._digits, tmp_poly, n)

        reduce_coeff(self._digits, self._precrel, self.prime_pow.fprime)
        self._precrel += 1
        return 0


cdef class pAdicLazyElement_muldigit(pAdicLazyElement):
    # compute x*y where x is assumed to be in (0,p)
    def __init__(self, parent, pAdicLazyElement_div x, pAdicLazyElement y):
        pAdicLazyElement.__init__(self, parent)
        self._x = <fmpz*>x._inverse
        self._y = y
        self._valuation = y._valuation

    cdef int _next_c(self):
        cdef slong n = self._valuation + self._precrel
        cdef int error = self._y._jump_c(n+1)
        if error:
            return error
        fmpz_mul(tmp_coeff, self._x, get_coeff(self._y._digits, n - self._y._valuation))
        iadd_coeff(self._digits, tmp_coeff, self._precrel)
        reduce_coeff(self._digits, self._precrel, self.prime_pow.fprime)
        self._precrel += 1
        return 0


# Division

cdef class pAdicLazyElement_div(pAdicLazyElement):
    def __cinit__(self):
        fmpz_init(self._inverse)

    def __dealloc__(self):
        fmpz_clear(self._inverse)

    def __init__(self, parent, pAdicLazyElement num, pAdicLazyElement denom):
        pAdicLazyElement.__init__(self, parent)
        self._num = num
        self._denom = denom
        if denom._valuation <= -maxordp:
            self._maxprec = maxordp + 1
        else:
            self._maxprec = denom._valuation + self._parent.default_prec()
        cdef int error = self._bootstrap_c()
        if error & ERROR_DIVISION:
            raise ZeroDivisionError("cannot divide by something indistinguishable from zero")
        if error:
            self._valuation = -maxordp

    cdef int _bootstrap_c(self):
        cdef int error
        cdef pAdicLazyElement num = self._num
        cdef pAdicLazyElement denom = self._denom
        cdef fmpz_t gcd

        while denom._valuation < self._maxprec and denom._precrel == 0:
            error = denom._next_c()
            if error:
                if error & ERROR_PRECISION:
                    error |= ERROR_DIVISION
                return error
            if self._maxprec > maxordp and denom._valuation > -maxordp:
                self._maxprec = denom._valuation + self._parent.default_prec()
        if denom._precrel == 0:
            return ERROR_ABANDON

        self._valuation = num._valuation - denom._valuation
        fmpz_gcdinv(gcd, self._inverse, get_coeff(denom._digits, 0), self.prime_pow.fprime)
        cdef pAdicLazyElement a = pAdicLazyElement_muldigit(self._parent, self, num)
        cdef pAdicLazyElement b = pAdicLazyElement_muldigit(self._parent, self, denom)
        error = b._next_c()
        if error:
            return ERROR_UNKNOWN
        b._valuation += 1
        b._precrel = 0
        fmpz_poly_shift_right(b._digits, b._digits, 1)
        self._definition = a - b * self
        return 0

    cdef int _next_c(self):
        cdef pAdicLazyElement definition = self._definition
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
            else:
                digit = get_coeff(definition._digits, self._precrel)
                fmpz_poly_set_coeff_fmpz(self._digits, self._precrel, digit);
                self._precrel += 1
        return 0


# Square root

cdef bint fmpz_modp_sqrt(fmpz* ans, fmpz* x, fmpz* p):
    # Need to do something better
    cdef Integer zx = PY_NEW(Integer)
    fmpz_get_mpz(zx.value, x)
    cdef Integer zp = PY_NEW(Integer)
    fmpz_get_mpz(zp.value, p)
    cdef Integer zans
    try:
        k = GF(zp)
        zans = ZZ(k(zx).sqrt(extend=False))
    except ValueError:
        return 1
    fmpz_set_mpz(ans, zans.value)
    return 0


cdef class pAdicLazyElement_sqrt(pAdicLazyElement):
    def __init__(self, parent, pAdicLazyElement x):
        pAdicLazyElement.__init__(self, parent)
        if parent.prime() == 2:
            raise NotImplementedError
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

    cdef int _bootstrap_c(self):
        cdef pAdicLazyElement x = self._x
        while x._valuation < self._maxprec and x._precrel == 0:
            error = x._next_c()
            if error:
                return error
            if self._maxprec > maxordp and x._valuation > -maxordp:
                self._maxprec = x._valuation + 2*self._parent.default_prec()
        if x._valuation > -maxordp:
            self._valuation = x._valuation >> 1
        if x._precrel == 0:
            return ERROR_ABANDON

        if x._valuation & 1 != 0:
            return ERROR_NOTSQUARE
        cdef fmpz_t digit
        fmpz_init(digit)
        if fmpz_modp_sqrt(digit, get_coeff(x._digits, 0), self.prime_pow.fprime):
            return ERROR_NOTSQUARE

        fmpz_poly_set_coeff_fmpz(self._digits, 0, digit)
        cdef parent = self._parent
        cdef slong val = self._valuation
        cdef Integer zd = PY_NEW(Integer)
        fmpz_get_mpz(zd.value, digit)
        cdef pAdicLazyElement u = pAdicLazyElement_shift(parent, self, val + 1, True)
        cdef pAdicLazyElement y = pAdicLazyElement_shift(parent, x, 2*val + 2, False)
        cdef pAdicLazyElement c = pAdicLazyElement_value(parent, zd*zd, shift=-2)
        cdef pAdicLazyElement d = pAdicLazyElement_value(parent, 2*zd, shift=-val-2)
        self._definition = (y + c - u*u) / d
        self._precrel = 1
        return 0

    cdef int _next_c(self):
        cdef pAdicLazyElement x = self._x
        cdef pAdicLazyElement definition = self._definition
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
            fmpz_poly_set_coeff_fmpz(self._digits, self._precrel, get_coeff(definition._digits, self._precrel))
            self._precrel += 1
        return 0



# Self-referent definitions
###########################

cdef class pAdicLazyElement_selfref(pAdicLazyElement):
    def __init__(self, parent, valuation):
        pAdicLazyElement.__init__(self, parent)
        self._valuation = valuation
        self._definition = None
        self._next = False

    cpdef set(self, pAdicLazyElement definition):
        if self._definition is not None:
            raise ValueError("this self-referent number is already defined")
        self._definition = definition
        cdef slong valsve = self._valuation
        cdef int error = self._jump_c(self._valuation + self._parent.default_prec())
        if error & ERROR_CIRCULAR:
            self._definition = None
            self._valuation = valsve
            self._precrel = 0
            self._next = False
        raise_error(error, permissive=True)

    def __eq__(self, other):
        if self._definition is None:
            self.set(other)
        return pAdicLazyElement.__eq__(self, other)

    cdef int _next_c(self):
        cdef pAdicLazyElement definition = self._definition
        cdef fmpz* digit
        cdef slong diffval
        cdef int error

        if definition is None:
            return ERROR_NOTDEFINED
        if self._next:
            return ERROR_CIRCULAR

        self._next = True
        error = definition._jump_c(self._valuation + self._precrel + 1)
        if not error:
            diffval = self._valuation - definition._valuation
            if diffval < 0:
                self._valuation = definition._valuation
            else:
                digit = get_coeff(definition._digits, self._precrel + diffval)
                if self._precrel == 0 and fmpz_is_zero(digit):
                    self._valuation += 1
                else:
                    fmpz_poly_set_coeff_fmpz(self._digits, self._precrel, digit);
                    self._precrel += 1
        self._next = False
        return error


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
