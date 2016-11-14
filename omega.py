
# *****************************************************************************
# Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# http://www.gnu.org/licenses/
# *****************************************************************************

from __future__ import print_function
from __future__ import absolute_import

import operator
from sage.groups.indexed_free_group import IndexedFreeAbelianGroup
from six import iteritems, itervalues


def HomogenousSymmetricFunction(j, x):
    r"""
    EXAMPLES::

        sage: P = PolynomialRing(ZZ, 'X', 3)
        sage: HomogenousSymmetricFunction(0, P.gens())
        1
        sage: HomogenousSymmetricFunction(1, P.gens())
        X0 + X1 + X2
        sage: HomogenousSymmetricFunction(2, P.gens())
        X0^2 + X0*X1 + X1^2 + X0*X2 + X1*X2 + X2^2
        sage: HomogenousSymmetricFunction(3, P.gens())
        X0^3 + X0^2*X1 + X0*X1^2 + X1^3 + X0^2*X2 +
        X0*X1*X2 + X1^2*X2 + X0*X2^2 + X1*X2^2 + X2^3
    """
    from sage.combinat.integer_vector import IntegerVectors
    return sum(prod(xx**pp for xx, pp in zip(x, p))
               for p in IntegerVectors(j, length=len(x)))


def Omega_numerator(a, n, m):
    r"""
    EXAMPLES::

        sage: L.<x0, x1, x2, y0, y1, y2> = LaurentPolynomialRing(QQ)
        sage: Omega_numerator(0, 1, 1)
        1
        sage: Omega_numerator(0, 2, 1)
        -x0*x1*y0 + 1
        sage: Omega_numerator(0, 1, 2)
        1
        sage: Omega_numerator(0, 3, 1)
        x0*x1*x2*y0^2 + x0*x1*x2*y0 - x0*x1*y0 - x0*x2*y0 - x1*x2*y0 + 1
        sage: Omega_numerator(0, 2, 2)
        x0^2*x1*y0*y1 + x0*x1^2*y0*y1 - x0*x1*y0*y1 - x0*x1*y0 - x0*x1*y1 + 1

        sage: Omega_numerator(-2, 1, 1)
        x0^2
        sage: Omega_numerator(-1, 1, 1)
        x0
        sage: Omega_numerator(1, 1, 1)
        -x0*y0 + y0 + 1
        sage: Omega_numerator(2, 1, 1)
        -x0*y0^2 - x0*y0 + y0^2 + y0 + 1
    """
    Y = LaurentPolynomialRing(
        QQ, ', '.join('y{}'.format(mm) for mm in range(m)))
    y = Y.gens()

    def P(n):
        if n == 1:
            L = LaurentPolynomialRing(Y, 'x0')
            x0 = L.gen()
            return x0**(-a) + \
                (prod(1 - x0*yy for yy in y) *
                 sum(HomogenousSymmetricFunction(j, y) * (1-x0**(j-a))
                     for j in srange(a))
                 if a > 0 else 0)
        else:
            Pprev = P(n-1)
            L = LaurentPolynomialRing(Pprev.parent(), 'x{}'.format(n-1))
            x1 = L.gen()
            x2 = Pprev.parent().gen()
            p1 = L(Pprev.subs({x2: x1}))
            p2 = L(Pprev)
            x2 = L({0: x2})
            q, r = (x1 * (1-x2) * prod(1 - x2*yy for yy in y) * p1 - \
                    x2 * (1-x1) * prod(1 - x1*yy for yy in y) * p2).quo_rem(x1 - x2)
            assert r == 0
            return q

    from itertools import chain
    XY = LaurentPolynomialRing(QQ, ', '.join(chain(
        iter('x{}'.format(nn) for nn in range(n)),
        iter('y{}'.format(mm) for mm in range(m)))))

    if m == 0:
        return XY(1 - (prod(factors_denominator) *
                       sum(HomogenousSymmetricFunction(j, XY.gens())
                           for j in srange(-a))
                       if a < 0 else 0))
    elif n == 0:
        return XY(sum(HomogenousSymmetricFunction(j, XY.gens())
                      for j in srange(a+1)))
    else:
        return XY(P(n))


def Omega_fundamental(a, x, y, group_factors=False):
    r"""
    Return `\Omega_{\ge}` of the expression specified by the input.

    .. MATH::

        \Omega_{\ge} \frac{\lambda^a}{
        (1 - x_1 \lambda) \dots (1 - x_n \lambda)
        (1 - y_1 / \lambda) \dots (1 - y_m / \lambda)

    INPUT:

    - ``a`` -- an integer.

    - ``x`` and ``y`` -- lists of laurent polynomials

    OUTPUT:

    A pair representing a quotient as follows: Its first component is the
    numerator as a laurent polynomial, its second component a factorization
    of the denominator as a list of laurent polynomials.

    EXAMPLES::

        sage: L.<x, y, z, w> = LaurentPolynomialRing(QQ)
        sage: Omega_fundamental(0, [x], [y])
        (1, (-x + 1, -x*y + 1))
        sage: Omega_fundamental(0, [x], [y, z])
        (1, (-x + 1, -x*y + 1, -x*z + 1))
        sage: Omega_fundamental(0, [x, y], [z])
        (-x*y*z + 1, (-x + 1, -y + 1, -x*z + 1, -y*z + 1))
        sage: Omega_fundamental(0, [x, y, z], [w])
        (x*y*z*w^2 + x*y*z*w - x*y*w - x*z*w - y*z*w + 1,
         (-x + 1, -y + 1, -z + 1, -x*w + 1, -y*w + 1, -z*w + 1))
        sage: Omega_fundamental(0, [x, y], [z, w])
        (x^2*y*z*w + x*y^2*z*w - x*y*z*w - x*y*z - x*y*w + 1,
         (-x + 1, -y + 1, -x*z + 1, -x*w + 1, -y*z + 1, -y*w + 1))

        sage: Omega_fundamental(-2, [x], [y])
        (x^2, (-x + 1, -x*y + 1))
        sage: Omega_fundamental(-1, [x], [y])
        (x, (-x + 1, -x*y + 1))
        sage: Omega_fundamental(1, [x], [y])
        (-x*y + y + 1, (-x + 1, -x*y + 1))
        sage: Omega_fundamental(2, [x], [y])
        (-x*y^2 - x*y + y^2 + y + 1, (-x + 1, -x*y + 1))
    """
    def flatten(z):
        return sum((tuple(zz) for zz in z), tuple())

    if group_factors:
        flat_x = flatten(x)
        flat_y = flatten(y)
    else:
        flat_x = x
        flat_y = y
        x = tuple((xx,) for xx in x)
        y = tuple((yy,) for yy in y)

    numerator = Omega_numerator(a, len(flat_x), len(flat_y))
    numerator = numerator.subs(
        {xi: xj for xi, xj in zip(numerator.parent().gens(), flat_x+flat_y)})

    Factors_denominator = \
        tuple(tuple(1 - xx for xx in gx) for gx in x) + \
        tuple(tuple(1 - xx*yy for xx in gx for yy in gy)
              for gx in x for gy in y)
    if not group_factors:
        factors_denominator = flatten(factors_denominator)

    return numerator, factors_denominator


def Omega_higher(a, z):
    r"""
    Return `\Omega_{\ge}` of the expression specified by the input.

    .. MATH::

        \Omega_{\ge} \frac{\lambda^a}{
        (1 - z_1 \lambda^{e_1}) \dots (1 - z_n \lambda^{e_n})

    INPUT:

    - ``a`` -- an integer.

    - ``z`` and ``y`` -- a lists with each entry either
      a laurent polynomial (implicit exponent `1`) or
      a pair of a laurent polynomial and an integer exponent.

    OUTPUT:

    A pair representing a quotient as follows: Its first component is the
    numerator as a laurent polynomial, its second component a factorization
    of the denominator as a list of laurent polynomials.

    EXAMPLES::

        sage: L.<x, y, z, w> = LaurentPolynomialRing(QQ)

        sage: Omega_higher(0, [(x, 1), (y, -2)])
        (1, (-x + 1, -x^2*y + 1))
        sage: Omega_higher(0, [(x, 1), (y, -3)])
        (1, (-x + 1, -x^3*y + 1))
        sage: Omega_higher(0, [(x, 1), (y, -4)])
        (1, (-x + 1, -x^4*y + 1))

        sage: Omega_higher(0, [(x, 2), (y, -1)])
        (x*y + 1, (-x + 1, -x*y^2 + 1))
        sage: Omega_higher(0, [(x, 3), (y, -1)])
        (x*y^2 + x*y + 1, (-x + 1, -x*y^3 + 1))
        sage: Omega_higher(0, [(x, 4), (y, -1)])
        (x*y^3 + x*y^2 + x*y + 1, (-x + 1, -x*y^4 + 1))

        sage: Omega_higher(0, [(x, 1), (y, 1), (z, -2)])
        (-x^2*y*z - x*y^2*z + x*y*z + 1,
         (-x + 1, -y + 1, -x^2*z + 1, -y^2*z + 1))
        sage: Omega_higher(0, [(x, 2), (y, -1), (z, -1)])
        (x*y*z + x*y + x*z + 1, (-x + 1, -x*y^2 + 1, -x*z^2 + 1))
        sage: Omega_higher(0, [(x, 2), (y, 1), (z, -1)])
        (-x*y*z^2 - x*y*z + x*z + 1, (-x + 1, -y + 1, -x*z^2 + 1, -y*z + 1))
    """
    class Factor(object):
        def __init__(self, zz):
            if isinstance(zz, (tuple, list)):
                self.value, self.exponent = zz
            else:
                self.value = zz
                self.exponent = 1
        def is_higher(self):
            return abs(self.exponent) > 1
        def z(self, positive=True):
            if (self.exponent > 0) != positive:
                return tuple()
            elif self.is_higher():
                rho = powers[abs(self.exponent)]
                w = self.var
                return tuple(L_high(rho**j * w)
                             for j in srange(abs(self.exponent)))
            else:
                return (L_high(self.value),)
        def x(self):
            return self.z(positive=True)
        def y(self):
            return self.z(positive=False)

    z = list(Factor(zz) for zz in z)

    # -2. create new (larger) laurent polynomial ring
    L_orig = z[0].value.parent()
    B_orig = L_orig.base_ring()
    exponents_pre = sorted(set(abs(factor.exponent)
                               for factor in z if factor.is_higher()))
    B_high = B_orig.extension(
        list(cyclotomic_polynomial(r) for r in exponents_pre),
        tuple('rho{}'.format(i) for i in srange(len(exponents_pre))))
    powers = dict(zip(exponents_pre, B_high.gens()))

    nv = len(tuple(None for factor in z if factor.is_higher()))
    variable_names_high = \
        tuple('Omega{}'.format(i) for i in srange(nv)) + \
        L_orig.variable_names()
    L_high = LaurentPolynomialRing(B_high, variable_names_high)

    # -1. rewrite factors with higher powers
    v = iter(L_high.gens())
    Omega_map = dict()
    for factor in z:
        if factor.is_higher():
            factor.var = next(v)
            Omega_map[factor.var] = factor

    # 0. apply Omega
    def nonempty(T):
        return tuple(t for t in T if t)
    numerator, factors_denominator = Omega_fundamental(
        a,
        nonempty(factor.x() for factor in z),
        nonempty(factor.y() for factor in z),
        group_factors=True)

    # 1. multiply grouped factors of denominator
    factors_denominator = tuple(prod(factor for factor in factors)
                                for factors in factors_denominator)

    # 2. substitute helper variable with actual value
    def subs_power(expression, var, exponent, value):
        r"""
        Substitute ``var^exponent`` by ``value`` in ``expression``.
        """
        p = tuple(var.dict().popitem()[0]).index(1)
        def subs_e(e):
            e = list(e)
            assert e[p] % exponent == 0
            e[p] = e[p] // exponent
            return tuple(e)
        parent = expression.parent()
        result = parent({subs_e(e): c for e, c in iteritems(expression.dict())})
        return result.subs({var: value})

    def subs_Omega(expression):
        r"""
        Substitute all helper variables by their actual values.
        """
        for var, factor in iteritems(Omega_map):
            expression = subs_power(expression, var,
                                    abs(factor.exponent), factor.value)
        return expression

    factors_denominator = tuple(subs_Omega(factor)
                                for factor in factors_denominator)

    from sage.rings.fraction_field import FractionField_generic
    if isinstance(numerator.parent(), FractionField_generic):
        numerator = subs_Omega(L_high(numerator.numerator())) / \
                    subs_Omega(L_high(numerator.denominator()))
    else:
        numerator = subs_Omega(numerator)

    return numerator, factors_denominator


def Omega(var, numerator, factors_denominator, operator=operator.ge):
    pass

