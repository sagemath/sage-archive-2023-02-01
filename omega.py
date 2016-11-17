
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


def _laurent_polynomial_ring_(n, m):
    from itertools import chain
    return LaurentPolynomialRing(QQ, ', '.join(chain(
        iter('x{}'.format(nn) for nn in range(n)),
        iter('y{}'.format(mm) for mm in range(m)))))


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

    XY = _laurent_polynomial_ring_(n, m)

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


def Omega_factors_denominator(n, m):
    r"""
    EXAMPLES::

        sage: Omega_factors_denominator(1, 1)
        ((-x0 + 1,), (-x0*y0 + 1,))
        sage: Omega_factors_denominator(1, 2)
        ((-x0 + 1,), (-x0*y0 + 1,), (-x0*y1 + 1,))
        sage: Omega_factors_denominator(2, 1)
        ((-x0 + 1,), (-x1 + 1,), (-x0*y0 + 1,), (-x1*y0 + 1,))
        sage: Omega_factors_denominator(3, 1)
        ((-x0 + 1,), (-x1 + 1,), (-x2 + 1,),
         (-x0*y0 + 1,), (-x1*y0 + 1,), (-x2*y0 + 1,))
        sage: Omega_factors_denominator(2, 2)
        ((-x0 + 1,), (-x1 + 1,), (-x0*y0 + 1,),
         (-x0*y1 + 1,), (-x1*y0 + 1,), (-x1*y1 + 1,))
    """
    if isinstance(n, tuple):
        x = n
        n = sum(x)
    else:
        x = tuple(1 for _ in range(n))
    if isinstance(m, tuple):
        y = m
        m = sum(y)
    else:
        y = tuple(1 for _ in range(m))

    XY = _laurent_polynomial_ring_(n, m)
    ixy = iter(XY.gens())
    x = tuple(tuple(next(ixy) for _ in range(nx)) for nx in x)
    y = tuple(tuple(next(ixy) for _ in range(my)) for my in y)

    return tuple(tuple(1 - xx for xx in gx) for gx in x) + \
           tuple(tuple(1 - xx*yy for xx in gx for yy in gy)
                 for gx in x for gy in y)


def Omega_higher(a, z):
    r"""
    EXAMPLES::

        sage: Omega_higher(0, [1, -2])
        (1, (-z0 + 1, -z0^2*z1 + 1))
        sage: Omega_higher(0, [1, -3])
        (1, (-z0 + 1, -z0^3*z1 + 1))
        sage: Omega_higher(0, [1, -4])
        (1, (-z0 + 1, -z0^4*z1 + 1))

        sage: Omega_higher(0, [2, -1])
        (z0*z1 + 1, (-z0 + 1, -z0*z1^2 + 1))
        sage: Omega_higher(0, [3, -1])
        (z0*z1^2 + z0*z1 + 1, (-z0 + 1, -z0*z1^3 + 1))
        sage: Omega_higher(0, [4, -1])
        (z0*z1^3 + z0*z1^2 + z0*z1 + 1, (-z0 + 1, -z0*z1^4 + 1))

        sage: Omega_higher(0, [1, 1, -2])
        (-z0^2*z1*z2 - z0*z1^2*z2 + z0*z1*z2 + 1,
         (-z0 + 1, -z1 + 1, -z0^2*z2 + 1, -z1^2*z2 + 1))
        sage: Omega_higher(0, [2, -1, -1])
        (z0*z1*z2 + z0*z1 + z0*z2 + 1, (-z0 + 1, -z0*z1^2 + 1, -z0*z2^2 + 1))
        sage: Omega_higher(0, [2, 1, -1])
        (-z0*z1*z2^2 - z0*z1*z2 + z0*z2 + 1, (-z0 + 1, -z1 + 1, -z0*z2^2 + 1, -z1*z2 + 1))
    """
    if not z or any(zz == 0 for zz in z):
        raise NotImplementedError

    x = tuple(zz for zz in z if zz > 0)
    y = tuple(-zz for zz in z if zz < 0)
    xy = x + y
    n = sum(x)
    m = sum(y)

    exponents = sorted(set(zz for zz in xy) - set([1]))
    B = QQ.extension(
        list(cyclotomic_polynomial(r) for r in exponents),
        tuple('rho{}'.format(i) for i in range(len(exponents))))
    L = LaurentPolynomialRing(B, ', '.join('z{}'.format(nn)
                                           for nn in range(len(z))))
    powers = dict(zip(exponents, iter(L(g) for g in B.gens())))
    powers[2] = L(-1)
    powers[1] = L(1)

    def subs_power(expression, var, exponent, value=None):
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
        if value is None:
            return result
        return result.subs({var: value})

    Z = L.change_ring(QQ)

    def de_power(expression):
        for zz, var in zip(z, L.gens()):
            if abs(zz) == 1:
                continue
            expression = subs_power(expression, var, abs(zz))
        return Z(expression)

    xy_var = _laurent_polynomial_ring_(n, m).gens()
    x_var = iter(xy_var[:n])
    y_var = iter(xy_var[n:])
    rules = {next(x_var) if zz > 0 else next(y_var):
             powers[abs(zz)]**j * var
             for zz, var in zip(z, L.gens()) for j in range(abs(zz))}
    factors_denominator = tuple(de_power(prod(f.subs(rules) for f in factors))
                                for factors in Omega_factors_denominator(x, y))

    numerator = de_power(Omega_numerator(a, n, m).subs(rules))

    return numerator, factors_denominator

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


def Omega(var, numerator, factors_denominator, operator=operator.ge):
    pass

