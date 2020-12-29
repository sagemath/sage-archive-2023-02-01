# ****************************************************************************
#       Copyright (C) 2021 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************


from sage.structure.element import coerce_binop

from sage.rings.all import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.padics.padic_generic_element import pAdicGenericElement

from sage.rings.padics.precision_error import PrecisionError


class pAdicLazyElement(pAdicGenericElement):
    def __init__(self, parent):
        pAdicGenericElement.__init__(self, parent)
        self._digits = [ ]
        self._precision = 0

    def prime(self):
        return self.parent().prime()

    def next(self):
        raise NotImplementedError("must be implemented in subclasses")

    def jump(self, prec):
        while self._precision < prec:
            self.next()

    def _repr_(self):
        try:
            self.jump(self.parent().precision_cap())
        except PrecisionError:
            pass
        if self._precision == 0:
            return "O(1)"
        s = "..."
        for i in range(self._precision - 1, -1, -1):
            if self.prime() > 9:
                s += "|"
            s += str(self._digits[i])
        return s

    def __getitem__(self, i):
        if isinstance(i, slice):
            if i.step is not None or i.stop is None:
                raise IndexError
            start = i.start
            stop = i.stop
        else:
            start = i
            stop = i + 1
        if start < 0:
            start = 0
        if stop < 0 or start >= stop:
            return ZZ(0)
        self.jump(stop)
        if isinstance(i, slice):
            return self._digits[start:stop]
        else:
            return self._digits[start]

    @coerce_binop
    def equal_at_current_precision(self, other):
        for i in range(min(self._precision, other._precision)):
            if self[i] != other[i]:
                return False
        return True

    def __eq__(self, other):
        return self.equal_at_current_precision(other)

    def __lshift__(self, s):
        if s == 0:
            return self
        return pAdicLazyElement_shift(self.parent(), self, ZZ(s))

    def __rshift__(self, s):
        if s == 0:
            return self
        return pAdicLazyElement_shift(self.parent(), self, ZZ(-s))

    def _add_(self, other):
        summands = [ ]
        coeffs = [ ]
        if isinstance(self, pAdicLazyElement_add):
            summands.extend(self._summands)
            coeffs.extend(self._coeffs)
        else:
            summands.append(self)
            coeffs.append(1)
        if isinstance(other, pAdicLazyElement_add):
            summands.extend(other._summands)
            coeffs.extend(other._coeffs)
        else:
            summands.append(other)
            coeffs.append(1)
        return pAdicLazyElement_add(self.parent(), summands, coeffs)

    def _sub_(self, other):
        summands = [ ]
        coeffs = [ ]
        if isinstance(self, pAdicLazyElement_add):
            summands.extend(self._summands)
            coeffs.extend(self._coeffs)
        else:
            summands.append(self)
            coeffs.append(1)
        if isinstance(other, pAdicLazyElement_add):
            summands.extend(other._summands)
            coeffs.extend([ -c for c in other._coeffs ])
        else:
            summands.append(other)
            coeffs.append(-1)
        return pAdicLazyElement_add(self.parent(), summands, coeffs)

    def _neg_(self):
        return pAdicLazyElement_add(self.parent(), self, [-1])

    def valuation_at_current_precision(self):
        for i in range(self._precision):
            if self._digits[i] != 0:
                return i
        return self._precision

    def _mul_(self, other):
        vs = self.valuation_at_current_precision()
        vo = other.valuation_at_current_precision()
        return pAdicLazyElement_mul(self.parent(), self >> vs, other >> vo) << (vs + vo)

    def _div_(self, other):
        v = self.valuation_at_current_precision()
        return pAdicLazyElement_div(self.parent(), self >> v, other) << v

    def __invert__(self):
        return pAdicLazyElement_div(self.parent(), None, self)

    def sqrt(self):
        return pAdicLazyElement_sqrt(self.parent(), self)

    def convert(self, ring, prec=None):
        if prec is None:
            prec = self.parent().prec()
        digits = self[:prec]
        value = 0
        for d in reversed(digits):
            value *= self.prime()
            value += d
        return ring(value)


# Value

class pAdicLazyElement_value(pAdicLazyElement):
    def __init__(self, parent, value, prec=None, maxprec=None):
        pAdicLazyElement.__init__(self, parent)
        self._digits = [ ZZ(value) ]
        self._maxprec = maxprec
        if prec is None:
            if maxprec is None:
                prec = parent.precision_cap()
            else:
                prec = maxprec
        try:
            self.jump(prec)
        except PrecisionError:
            pass

    def next(self):
        if self._maxprec is not None and self._precision >= self._maxprec:
            raise PrecisionError("not enough precision")
        carry, digit = self._digits[-1].quo_rem(self.prime())
        self._digits[-1] = digit
        self._digits.append(carry)
        self._precision += 1

# Random

class pAdicLazyElement_random(pAdicLazyElement):
    def next(self):
        d = ZZ.random_element(0, self.prime())
        self._digits.append(d)
        self._precision += 1


# Shift

class pAdicLazyElement_shift(pAdicLazyElement):
    def __init__(self, parent, x, s):
        pAdicLazyElement.__init__(self, parent)
        self._x = x
        self._s = s
        if s > 0:
            self._digits = s * [ZZ(0)]
            self._precision = s

    def next(self):
        n = self._precision
        self._digits.append(self._x[n - self._s])
        self._precision += 1


# Addition

class pAdicLazyElement_add(pAdicLazyElement):
    def __init__(self, parent, summands, coeffs):
        pAdicLazyElement.__init__(self, parent)
        self._summands = summands
        self._coeffs = coeffs
        self._digits = [ ZZ(0) ]

    def next(self):
        n = self._precision
        v = self._digits[-1]
        for i in range(len(self._summands)):
            v += self._coeffs[i] * self._summands[i][n]
        carry, digit = v.quo_rem(self.prime())
        self._digits[-1] = digit
        self._digits.append(carry)
        self._precision += 1


# Multiplication

class pAdicLazyElement_mul(pAdicLazyElement):
    _Zt = PolynomialRing(ZZ, 't')

    def __init__(self, parent, x, y):
        pAdicLazyElement.__init__(self, parent)
        self._x = x
        self._y = y
        self._digits = [ ZZ(0) ]

    def next(self):
        x = self._x
        y = self._y
        Zt = self._Zt
        n = self._precision
        m = n + 2
        ell = 0
        S = Zt(0)
        while m > 1:
            P = Zt(x[2**ell - 1 : 2**(ell+1) - 1])
            Q = Zt(y[(m-1)*2**ell - 1 : m*2**ell - 1])
            S += P * Q
            if m > 2:
                P = Zt(y[2**ell - 1 : 2**(ell+1) - 1])
                Q = Zt(x[(m-1)*2**ell - 1 : m*2**ell - 1])
                S += P * Q
            if m % 2 == 1:
                break
            m = m // 2
            ell += 1
        i = 0
        j = n
        digits = self._digits
        while i <= S.degree() and j < len(digits):
            digits[j] += S[i]
            i += 1
            j += 1
        while i <= S.degree():
            digits.append(S[i])
            i += 1
        carry, digit = digits[n].quo_rem(self.prime())
        digits[n] = digit
        if len(digits) < n+2:
            digits.append(carry)
        else:
            digits[n+1] += carry
        self._precision += 1


# Division

class pAdicLazyElement_div(pAdicLazyElement):
    def __init__(self, parent, x, y):
        pAdicLazyElement.__init__(self, parent)
        self._x = x
        self._y = y
        try:
            self._bootstrap()
            self._ready = True
        except PrecisionError:
            self._ready = False

    def _bootstrap(self):
        parent = self.parent()
        d, u, _ = self._y[0].xgcd(parent.prime())
        if d > 1:
            raise ValueError("divisor is not a unit")
        denom = pAdicLazyElement_add(parent, [self._y], [u])
        if self._x is None:
            num = pAdicLazyElement_value(parent, u)
        else:
            num = pAdicLazyElement_add(parent, [self._x], [u])
        self._definition = num - (denom >> 1) * (self << 1)

    def next(self):
        if not self._ready:
            self._bootstrap()
            self._ready = True
        self._digits.append(self._definition[self._precision])
        self._precision += 1


# Square root

class pAdicLazyElement_sqrt(pAdicLazyElement):
    def __init__(self, parent, x):
        if parent.prime() == 2:
            raise NotImplementedError
        pAdicLazyElement.__init__(self, parent)
        self._x = x
        try:
            self._bootstrap()
            self._ready = True
        except PrecisionError:
            self._ready = False

    def _bootstrap(self):
        parent = self.parent()
        a = parent.residue_field()(self._x[0])
        if a == 0:
            raise NotImplementedError
        b = ZZ(a.sqrt(extend=False))
        c = b + (self._x - ZZ(a)) / (2*b)
        s1 = self >> 1
        self._definition = c - ((s1**2 / (2*b)) << 2)

    def next(self):
        if not self._ready:
            self._bootstrap()
            self._ready = True
        self._digits.append(self._definition[self._precision])
        self._precision += 1


# Self-referent definition

class pAdicLazyElement_selfref(pAdicLazyElement):
    def __init__(self, parent):
        pAdicLazyElement.__init__(self, parent)
        self._definition = None
        self._next = False

    def set(self, definition):
        if self._definition is not None:
            raise ValueError("this self-referent number is already defined")
        self._definition = definition
        try:
            self.next()
        except PrecisionError:
            pass
        except RecursionError:
            self._definition = None
            raise

    def __eq__(self, other):
        if self._definition is None:
            self.set(other)
            return True
        return pAdicLazyElement.__eq__(self, other)

    def next(self):
        if self._definition is None:
            raise PrecisionError("this self-referent number is not defined")
        if self._next:
            raise RecursionError("definition looks circular")
        self._next = True
        try:
            self._digits.append(self._definition[self._precision])
            self._precision += 1
        finally:
            self._next = False
