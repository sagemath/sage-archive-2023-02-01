r"""
This module provides basic support for lazy p-adic integers.

    sage: R = ZpL(5, print_mode="digits")
    sage: R
    5-adic Ring with lazy precision

One creates elements as usual::

    sage: a = R(17/42)
    sage: a
    ...00244200244200244201

    sage: R.random_element()  # random
    ...21013213133412431402

By default, 20 digits of precision are computed (and printed).
If more (or less) digits are needed, one can specify it as follows::

    sage: b = R(42/17, prec=30)
    sage: b
    ...104201213402432310420121340301

Alternatively, one can increase the precision using the method meth:`jump`::

    sage: a.jump(30)
    sage: a
    ...244200244200244200244200244201

Standard operations are implemented::

    sage: a + b
    ...03232011214322140002
    sage: a - b
    ...42311334324023403400

    sage: a * b
    ...00000000000000000001
    sage: a / b
    ...12442142113021233401

We observe again that only 20 digits are computed.
If we need more, we have to create a new variable::

    sage: c = a / b
    sage: c.jump(40)
    sage: c
    ...4230030104200433321312442142113021233401

Note that this automatically increases the precision on `a` and `b`::

    sage: a
    ...4200244200244200244200244200244200244201
    sage: b
    ...2134024323104201213402432310420121340301

Equality test works but equality is only checked up to the
*minimal* current precision of the elements::

    sage: c == R(289/1764, prec=100)
    True
    sage: c == R(289/1764 + 5^50, prec=100)
    True

    sage: c.jump(100)
    sage: c == R(289/1764 + 5^50, prec=100)
    False
    sage: c == R(289/1764 + 5^50)
    True

A quite interesting feature with lazy p-adics is the possibility to
create (in somes cases) self-referent numbers. Here is an example.
We first declare a new variable as follows::

    sage: x = R()
    sage: x
    O(1)

We then write down an equation satisfied by `x`::

    sage: x == 1 + 5*x^2
    True

The variable `x` now contains the unique solution of the above equation::

    sage: x
    ...04222412141121000211

This works because the `n`-th digit of the right hand size of the
defining equation only involves the `i`-th digits of `x` with `i < n`
(this is due to the factor `5`).

As a comparison, the following produces an error::

    sage: y = R()
    sage: y == 1 + 3*y^2
    Traceback (most recent call last):
    ...
    RecursionError: definition looks circular

Self-referent definitions also work with systems of equations::

    sage: u = R(); v = R(); w = R()

    sage: u == 1 + 2*v + 3*w^2 + 5*u*v*w
    True
    sage: v == 2 + 4/w^3 + 10*(u + v + w)^2
    True
    sage: w == 3 + 25*(u*v + v*w + u*w)
    True

    sage: u
    ...21111231040212114001
    sage: v
    ...04313122332303332234
    sage: w
    ...30020001022124410403

"""

# ****************************************************************************
#       Copyright (C) 2021 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************


from sage.rings.all import ZZ, QQ
from sage.rings.infinity import Infinity
from sage.rings.padics.generic_nodes import pAdicRingBaseGeneric, pAdicFieldBaseGeneric
from sage.rings.padics.padic_generic_element import pAdicGenericElement
from sage.rings.padics.pow_computer_flint import PowComputer_flint

from sage.rings.padics.padic_lazy_element import *


class pAdicRingLazy(pAdicRingBaseGeneric):
    def __init__(self, p, prec, print_mode, names):
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, None)
        self._zero_element = pAdicLazyElement_zero(self)
        self._default_prec = ZZ(prec)
        self.prime_pow = PowComputer_flint(p, 1, 1, 1, False)

    def _prec_type(self):
        return "lazy"

    def is_lazy(self):
        return True

    def default_prec(self):
        return self._default_prec

    def precision_cap(self):
        return Infinity

    def _element_constructor_(self, x):
        if isinstance(x, pAdicLazyElement):
            return x
        if isinstance(x, pAdicGenericElement):
            return pAdicLazyElement_value(self, ZZ(x), maxprec=x.precision_absolute())
        if x in ZZ:
            if x == 0:
                return self._zero_element
            elif x == 1:
                return pAdicLazyElement_one(self)
            else:
                return pAdicLazyElement_value(self, ZZ(x))
        raise TypeError("unable to convert '%s' to a lazy %s-adic number" % (x, self.prime()))

    def selfref(self, valuation=0):
        if valuation not in ZZ or valuation < 0:
            raise ValueError("valuation must be a nonnegative integer")
        valuation = ZZ(valuation)
        if valuation >= MAXORDP:
            raise OverflowError("valuation is too large (maximum is %s)" % MAXORDP)
        return pAdicLazyElement_selfref(self, valuation)

    def _an_element_(self):
        return self._zero_element

    def random_element(self):
        return pAdicLazyElement_random(self, True)


class pAdicFieldLazy(pAdicFieldBaseGeneric):
    def __init__(self, p, prec, print_mode, names):
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, None)
        self._zero_element = pAdicLazyElement_zero(self)
        self._default_prec = ZZ(prec)
        self.prime_pow = PowComputer_flint(p, 1, 1, 1, True)

    def _prec_type(self):
        return "lazy"

    def is_lazy(self):
        return True

    def default_prec(self):
        return self._default_prec

    def precision_cap(self):
        return Infinity

    def _element_constructor_(self, x):
        if isinstance(x, pAdicLazyElement):
            if x.parent() is self:
                return x
            return pAdicLazyElement_shift(self, x, 0, False)
        if isinstance(x, pAdicGenericElement):
            return pAdicLazyElement_value(self, ZZ(x), maxprec=x.precision_absolute())
        if x in ZZ:
            if x == 0:
                return self._zero_element
            elif x == 1:
                return pAdicLazyElement_one(self)
            else:
                return pAdicLazyElement_value(self, ZZ(x))
        if x in QQ:
            num = x.numerator()
            denom = x.denominator()
            if denom % self.prime() == 0:
                raise ValueError("denominator is not a unit")
            num = pAdicLazyElement_value(self, num)
            denom = pAdicLazyElement_value(self, denom)
            return num / denom
        raise TypeError("unable to convert '%s' to a lazy %s-adic number" % (x, self.prime()))

    def selfref(self, valuation=0):
        if valuation not in ZZ:
            raise ValueError("valuation must be an integer")
        valuation = ZZ(valuation)
        if valuation >= MAXORDP:
            raise OverflowError("valuation is too large (maximum is %s)" % MAXORDP)
        return pAdicLazyElement_selfref(self, valuation)

    def _an_element_(self):
        return self._zero_element

    def random_element(self, integral=True):
        return pAdicLazyElement_random(self, integral)

    def sum(self, summands):
        n = len(summands)
        signs = n * [True]
        return pAdicLazyElement_add(self, summands, signs)
