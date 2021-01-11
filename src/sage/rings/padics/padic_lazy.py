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

By default, 20 digits of precision are printed.
If more digits are needed, one can increase the precision by using the
meth:`jump`::

    sage: a.jump(30)
    sage: a
    ...244200244200244200244200244201

Standard operations are implemented::

    sage: b = R(42/17)

    sage: a + b
    ...03232011214322140002
    sage: a - b
    ...42311334324023403400

    sage: a * b
    ...00000000000000000001
    sage: a / b
    ...12442142113021233401

We observe again that only 20 digits are printed, even if the precision
on the operands is higher::

    sage: b.jump(30)
    sage: a
    ...244200244200244200244200244201
    sage: b
    ...104201213402432310420121340301
    sage: a / b
    ...12442142113021233401

If more digits are needed, we have to create a new variable::

    sage: c = a / b
    sage: c.jump(40)
    sage: c
    ...4230030104200433321312442142113021233401

Note that this automatically increases the precision on `a` and `b`::

    sage: a
    ...4200244200244200244200244200244200244201
    sage: b
    ...2134024323104201213402432310420121340301

::

A quite interesting feature with lazy p-adics is the possibility to
create (in somes cases) self-referent numbers. Here is an example.
We first declare a new variable as follows::

    sage: x = R.selfref()
    sage: x
    ...?.0

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

    sage: y = R.selfref()
    sage: y == 1 + 3*y^2
    Traceback (most recent call last):
    ...
    RecursionError: definition looks circular

Self-referent definitions also work with systems of equations::

    sage: u = R.selfref()
    sage: v = R.selfref()
    sage: w = R.selfref()

    sage: u == 1 + 2*v + 3*w^2 + 5*u*v*w
    True
    sage: v == 2 + 4*w + sqrt(1 + 5*u + 10*v + 15*w)
    True
    sage: w == 3 + 25*(u*v + v*w + u*w)
    True

    sage: u
    ...31203130103131131433
    sage: v
    ...33441043031103114240
    sage: w
    ...30212422041102444403

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

from sage.rings.padics.padic_lazy_element import pAdicLazyElement
from sage.rings.padics.padic_lazy_element import pAdicLazyElement_zero
from sage.rings.padics.padic_lazy_element import pAdicLazyElement_random
from sage.rings.padics.padic_lazy_element import pAdicLazyElement_selfref
from sage.rings.padics.padic_lazy_element import pAdicLazyElement_teichmuller
from sage.rings.padics.padic_lazy_element import element_constructor


class pAdicRingLazy(pAdicRingBaseGeneric):
    def __init__(self, p, prec, print_mode, names):
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicLazyElement)
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

    def _an_element_(self):
        return self._zero_element

    def _element_constructor_(self, x):
        return element_constructor(self, x)

    def selfref(self, start_val=0):
        if start_val not in ZZ or start_val < 0:
            raise ValueError("valuation must be a nonnegative integer")
        start_val = ZZ(start_val)
        if start_val >= MAXORDP:
            raise OverflowError("valuation is too large (maximum is %s)" % MAXORDP)
        return pAdicLazyElement_selfref(self, start_val)

    def random_element(self):
        return pAdicLazyElement_random(self, True)

    def teichmuller(self, x):
        return pAdicLazyElement_teichmuller(self, ZZ(x))

    def teichmuller_system(self):
        R = self.residue_class_field()
        return [ self.teichmuller(ZZ(i)) for i in R if i != 0 ]


class pAdicFieldLazy(pAdicFieldBaseGeneric):
    def __init__(self, p, prec, print_mode, names):
        pAdicRingBaseGeneric.__init__(self, p, prec, print_mode, names, pAdicLazyElement)
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

    def _coerce_map_from_(self, R):
        if isinstance(R, pAdicRingLazy) and self is R.fraction_field():
            return True

    def _an_element_(self):
        return self._zero_element

    def _element_constructor_(self, x):
        return element_constructor(self, x)

    def selfref(self, start_val=0):
        if start_val not in ZZ:
            raise ValueError("valuation must be an integer")
        start_val = ZZ(start_val)
        if start_val >= MAXORDP:
            raise OverflowError("valuation is too large (maximum is %s)" % MAXORDP)
        return pAdicLazyElement_selfref(self, start_val)

    def random_element(self, integral=False):
        return pAdicLazyElement_random(self, integral)

    def teichmuller(self, x):
        return pAdicLazyElement_teichmuller(self, ZZ(x))

    def teichmuller_system(self):
        R = self.residue_class_field()
        return [ self.teichmuller(ZZ(i)) for i in R if i != 0 ]
