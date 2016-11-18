
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
    if n + m == 0:
        return QQ, tuple()
    L = LaurentPolynomialRing(QQ, ', '.join(chain(
        iter('x{}'.format(nn) for nn in range(n)),
        iter('y{}'.format(mm) for mm in range(m)))))
    return L, L.gens()


@cached_function
def Omega_numerator(a, n, m):
    r"""
    Return the numerator of `\Omega_{\ge}` of the expression
    specified by the input.

    To be more precise, this calculates

    .. MATH::

        \Omega_{\ge} \frac{\lambda^a}{
        (1 - x_1 \lambda) \dots (1 - x_n \lambda)
        (1 - y_1 / \lambda) \dots (1 - y_m / \lambda)

    and returns its numerator.

    INPUT:

    - ``a`` -- an integer.

    - ``n`` and ``m`` -- nonnegative integers.

    OUTPUT:

    A laurent polynomial.

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

    TESTS::

        sage: Omega_factors_denominator(0, 0)
        ()
        sage: Omega_numerator(0, 0, 0)
        1
        sage: Omega_numerator(+2, 0, 0)
        1
        sage: Omega_numerator(-2, 0, 0)
        0

        sage: Omega_factors_denominator(1, 0)
        ((1 - x0,),)
        sage: Omega_numerator(0, 1, 0)
        1
        sage: Omega_numerator(+2, 1, 0)
        1
        sage: Omega_numerator(-2, 1, 0)
        x0^2

        sage: Omega_factors_denominator(0, 1)
        ()
        sage: Omega_numerator(0, 0, 1)
        1
        sage: Omega_numerator(+2, 0, 1)
        1 + y0 + y0^2
        sage: Omega_numerator(-2, 0, 1)
        0
    """
    if m == 0:
        Y = QQ
        y = tuple()
    else:
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

    XY, xy_vars = _laurent_polynomial_ring_(n, m)

    if m == 0:
        return XY(1 - (prod(prod(f) for f in Omega_factors_denominator(n, m)) *
                       sum(HomogenousSymmetricFunction(j, xy_vars)
                           for j in srange(-a))
                       if a < 0 else 0))
    elif n == 0:
        return XY(sum(HomogenousSymmetricFunction(j, xy_vars)
                      for j in srange(a+1)))
    else:
        return XY(P(n))


@cached_function
def Omega_factors_denominator(n, m):
    r"""
    Return the denominator of `\Omega_{\ge}` of the expression
    specified by the input.

    To be more precise, this calculates

    .. MATH::

        \Omega_{\ge} \frac{1}{
        (1 - x_1 \lambda) \dots (1 - x_n \lambda)
        (1 - y_1 / \lambda) \dots (1 - y_m / \lambda)

    and returns a factorization of its denominator.

    INPUT:

    - ``n`` and ``m`` -- nonnegative integers or
      tuples of nonnegative integers. The latter specifys how the factors
      of the result are grouped. An integer (former case) corresponds to
      a tuple consisting of that many `1`s.

    OUTPUT:

    A factorization of the denominator as
    a tuple of tuples of laurent polynomials.

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

    TESTS::

        sage: Omega_factors_denominator(0, 0)
        ()
        sage: Omega_factors_denominator(1, 0)
        ((1 - x0,),)
        sage: Omega_factors_denominator(0, 1)
        ()
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

    ixy = iter(_laurent_polynomial_ring_(n, m)[1])
    x = tuple(tuple(next(ixy) for _ in range(nx)) for nx in x)
    y = tuple(tuple(next(ixy) for _ in range(my)) for my in y)

    return tuple(tuple(1 - xx for xx in gx) for gx in x) + \
           tuple(tuple(1 - xx*yy for xx in gx for yy in gy)
                 for gx in x for gy in y)


@cached_function
def Omega_higher(a, exponents):
    r"""
    Return `\Omega_{\ge}` of the expression specified by the input.

    To be more precise, this calculates

    .. MATH::

        \Omega_{\ge} \frac{\lambda^a}{
        (1 - z_1 \lambda^{e_1}) \dots (1 - z_n \lambda^{e_n})

    and returns its numerator and a factorization of its denominator.

    INPUT:

    - ``a`` -- an integer.

    - ``exponents`` -- a tuple of integers.

    OUTPUT:

    A pair representing a quotient as follows: Its first component is the
    numerator as a laurent polynomial, its second component a factorization
    of the denominator as a tuple of laurent polynomials.

    EXAMPLES::

        sage: Omega_higher(0, (1, -2))
        (1, (-z0 + 1, -z0^2*z1 + 1))
        sage: Omega_higher(0, (1, -3))
        (1, (-z0 + 1, -z0^3*z1 + 1))
        sage: Omega_higher(0, (1, -4))
        (1, (-z0 + 1, -z0^4*z1 + 1))

        sage: Omega_higher(0, (2, -1))
        (z0*z1 + 1, (-z0 + 1, -z0*z1^2 + 1))
        sage: Omega_higher(0, (3, -1))
        (z0*z1^2 + z0*z1 + 1, (-z0 + 1, -z0*z1^3 + 1))
        sage: Omega_higher(0, (4, -1))
        (z0*z1^3 + z0*z1^2 + z0*z1 + 1, (-z0 + 1, -z0*z1^4 + 1))

        sage: Omega_higher(0, (1, 1, -2))
        (-z0^2*z1*z2 - z0*z1^2*z2 + z0*z1*z2 + 1,
         (-z0 + 1, -z1 + 1, -z0^2*z2 + 1, -z1^2*z2 + 1))
        sage: Omega_higher(0, (2, -1, -1))
        (z0*z1*z2 + z0*z1 + z0*z2 + 1, (-z0 + 1, -z0*z1^2 + 1, -z0*z2^2 + 1))
        sage: Omega_higher(0, (2, 1, -1))
        (-z0*z1*z2^2 - z0*z1*z2 + z0*z2 + 1,
         (-z0 + 1, -z1 + 1, -z0*z2^2 + 1, -z1*z2 + 1))
    """
    if not exponents or any(e == 0 for e in exponents):
        raise NotImplementedError

    x = tuple(e for e in exponents if e > 0)
    y = tuple(-e for e in exponents if e < 0)
    n = sum(x)
    m = sum(y)

    xy = sorted(set(x + y) - set([1]))
    B = QQ.extension(
        list(cyclotomic_polynomial(r) for r in xy),
        tuple('rho{}'.format(i) for i in range(len(xy))))
    L = LaurentPolynomialRing(B, ', '.join('z{}'.format(nn)
                                           for nn in range(len(exponents))))
    powers = dict(zip(xy, iter(L(g) for g in B.gens())))
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
        for e, var in zip(exponents, L.gens()):
            if abs(e) == 1:
                continue
            expression = subs_power(expression, var, abs(e))
        return Z(expression)

    xy_vars = _laurent_polynomial_ring_(n, m)[1]
    x_vars = iter(xy_vars[:n])
    y_vars = iter(xy_vars[n:])
    rules = {next(x_vars) if e > 0 else next(y_vars):
             powers[abs(e)]**j * var
             for e, var in zip(exponents, L.gens()) for j in range(abs(e))}
    factors_denominator = tuple(de_power(prod(f.subs(rules) for f in factors))
                                for factors in Omega_factors_denominator(x, y))

    numerator = de_power(Omega_numerator(a, n, m).subs(rules))

    return numerator, factors_denominator


def Omega(var, expression, denominator=None, op=operator.ge):
    r"""
    Return `\Omega_{\mathrm{op}}` of ``expression`` with respect to ``var``.

    To be more precise, this calculates

    .. MATH::

        \Omega_{\mathrm{op}} \frac{n}{d_1 \dots d_n}

    for the numerator `n` and the factors `d_1`, ..., `d_n` of
    the denominator, all of which are laurent polynomials in ``var``
    and returns a (partial) factorization of the result.

    INPUT:

    - ``var`` -- a variable or a representation string of a variable.

    - ``expression`` -- an element of the quotient field of some
      laurent polynomials. If ``denominator`` is specified, then
      this laurent polynomial is interpreted as the numerator of the
      expression.

    - ``denominator`` -- a laurent polynomial or a
      :class:`~sage.structure.factorization.Factorization` (consisting
      of laurent polynomial factors) or a tuple/list of factors (laurent
      polynomials).

    - ``op`` -- (default: ``operator.ge``) an operator.

    OUTPUT:

    A (partial) :class:`~sage.structure.factorization.Factorization`
    of the result whose factors are laurent polynomials.

    .. NOTE::

        The numerator of the result may not be factored.

    EXAMPLES::

        sage: L.<mu, x, y, z, w> = LaurentPolynomialRing(QQ)

        sage: Omega(mu, 1, [1 - x*mu, 1 - y/mu])
        1 * (-x + 1)^-1 * (-x*y + 1)^-1

        sage: Omega(mu, 1, [1 - x*mu, 1 - y/mu, 1 - z/mu])
        1 * (-x + 1)^-1 * (-x*y + 1)^-1 * (-x*z + 1)^-1
        sage: Omega(mu, 1, [1 - x*mu, 1 - y*mu, 1 - z/mu])
        (-x*y*z + 1) * (-x + 1)^-1 * (-y + 1)^-1 * (-x*z + 1)^-1 * (-y*z + 1)^-1
        sage: Omega(mu, 1, [1 - x*mu, 1 - y/mu^2])
        1 * (-x + 1)^-1 * (-x^2*y + 1)^-1
        sage: Omega(mu, 1, [1 - x*mu^2, 1 - y/mu])
        (x*y + 1) * (-x + 1)^-1 * (-x*y^2 + 1)^-1

        sage: Omega(mu, 1, [1 - x*mu, 1 - y*mu, 1 - z/mu^2])
        (-x^2*y*z - x*y^2*z + x*y*z + 1) *
        (-x + 1)^-1 * (-y + 1)^-1 * (-x^2*z + 1)^-1 * (-y^2*z + 1)^-1

        sage: Omega(mu, 1, [1 - x*mu, 1 - y/mu^3])
        1 * (-x + 1)^-1 * (-x^3*y + 1)^-1
        sage: Omega(mu, 1, [1 - x*mu, 1 - y/mu^4])
        1 * (-x + 1)^-1 * (-x^4*y + 1)^-1
        sage: Omega(mu, 1, [1 - x*mu^3, 1 - y/mu])
        (x*y^2 + x*y + 1) * (-x + 1)^-1 * (-x*y^3 + 1)^-1
        sage: Omega(mu, 1, [1 - x*mu^4, 1 - y/mu])
        (x*y^3 + x*y^2 + x*y + 1) * (-x + 1)^-1 * (-x*y^4 + 1)^-1

        sage: Omega(mu, 1, [1 - x*mu^2, 1 - y/mu, 1 - z/mu])
        (x*y*z + x*y + x*z + 1) *
        (-x + 1)^-1 * (-x*y^2 + 1)^-1 * (-x*z^2 + 1)^-1
        sage: Omega(mu, 1, [1 - x*mu^2, 1 - y*mu, 1 - z/mu])
        (-x*y*z^2 - x*y*z + x*z + 1) *
        (-x + 1)^-1 * (-y + 1)^-1 * (-x*z^2 + 1)^-1 * (-y*z + 1)^-1

        sage: Omega(mu, 1, [1 - x*mu, 1 - y*mu, 1 - z*mu, 1 - w/mu])
        (x*y*z*w^2 + x*y*z*w - x*y*w - x*z*w - y*z*w + 1) *
        (-x + 1)^-1 * (-y + 1)^-1 * (-z + 1)^-1 *
        (-x*w + 1)^-1 * (-y*w + 1)^-1 * (-z*w + 1)^-1
        sage: Omega(mu, 1, [1 - x*mu, 1 - y*mu, 1 - z/mu, 1 - w/mu])
        (x^2*y*z*w + x*y^2*z*w - x*y*z*w - x*y*z - x*y*w + 1) *
        (-x + 1)^-1 * (-y + 1)^-1 *
        (-x*z + 1)^-1 * (-x*w + 1)^-1 * (-y*z + 1)^-1 * (-y*w + 1)^-1

        sage: Omega(mu, mu^-2, [1 - x*mu, 1 - y/mu])
        x^2 * (-x + 1)^-1 * (-x*y + 1)^-1
        sage: Omega(mu, mu^-1, [1 - x*mu, 1 - y/mu])
        x * (-x + 1)^-1 * (-x*y + 1)^-1
        sage: Omega(mu, mu, [1 - x*mu, 1 - y/mu])
        (-x*y + y + 1) * (-x + 1)^-1 * (-x*y + 1)^-1
        sage: Omega(mu, mu^2, [1 - x*mu, 1 - y/mu])
        (-x*y^2 - x*y + y^2 + y + 1) * (-x + 1)^-1 * (-x*y + 1)^-1

    TESTS::

        sage: Omega(mu, 1, [1 - x*mu])
        1 * (-x + 1)^-1
        sage: Omega(mu, 1, [1 - x/mu])
        1
        sage: Omega(mu, 0, [1 - x*mu])
        0
        sage: Omega(mu, L(1), [])
        1
        sage: Omega(mu, L(0), [])
        0
        sage: Omega(mu, 2, [])
        2
    """
    if op != operator.ge:
        raise NotImplementedError('At the moment, only Omega_ge is implemented.')

    if denominator is None:
        numerator = expression.numerator()
        denominator = expression.denominator()
    else:
        numerator = expression

    from sage.arith.misc import factor
    from sage.structure.factorization import Factorization

    if isinstance(denominator, (list, tuple)):
        factors_denominator = denominator
    else:
        if not isinstance(denominator, Factorization):
            denominator = factor(denominator)
        if not denominator.is_integral():
            raise ValueError('Factorization {} of the denominator '
                             'contains negative exponents.'.format(denominator))
        numerator *= denominator.unit()
        factors_denominator = tuple(factor
                                    for factor, exponent in denominator
                                    for _ in range(exponent))

    if not factors_denominator:
        try:
            var = numerator.parent()(var)
        except (TypeError, ValueError):
            return Factorization([(numerator, 1)])
        else:
            return Factorization([(numerator.subs({var: 1}), 1)])
    R = factors_denominator[0].parent()
    var = repr(var)
    L0 = LaurentPolynomialRing(
        R.base_ring(), tuple(v for v in R.variable_names() if v != var))
    L = LaurentPolynomialRing(L0, var)
    numerator = L(numerator)
    if numerator == 0:
        return Factorization([], unit=numerator)
    factors_denominator = tuple(L(factor) for factor in factors_denominator)

    def decode_factor(factor):
        D = factor.dict()
        if len(D) != 2 or D.get(0, 0) != 1:
            raise NotImplementedError('Cannot handle factor {}'.format(factor))
        D.pop(0)
        exponent, coefficient = next(iteritems(D))
        return -coefficient, exponent
    # Below we sort to make the caching more efficient. Doing this here
    # (in contrast to directly in Omega_higher) results in much cleaner
    # code and prevents an additional substitution or passing of a permutation.
    decoded_factors = tuple(decode_factor(factor)
                            for factor in factors_denominator)

    result_numerator = 0
    result_factors_denominator = None
    for a, c in iteritems(numerator.dict()):
        n, fd = _Omega_(a, decoded_factors)
        if result_factors_denominator is None:
            result_factors_denominator = fd
        else:
            assert result_factors_denominator == fd
        result_numerator += c * n

    return Factorization([(result_numerator, 1)] +
                         list((f, -1) for f in result_factors_denominator),
                         sort=False)


def _Omega_(a, decoded_factors):
    r"""
    Helper function for :func:`Omega` which accesses the low level functions
    and does the substituting.

    INPUT:

    - ``a`` -- an integer.

    - ``decoded_factors`` -- a tuple or list of pairs `(z, e)` representing
      a factor `1 - \lambda^e z`.

    OUTPUT:

    A pair representing a quotient as follows: Its first component is the
    numerator as a laurent polynomial, its second component a factorization
    of the denominator as a tuple of laurent polynomials.
    """
    # Below we sort to make the caching more efficient. Doing this here
    # (in contrast to directly in Omega_higher) results in much cleaner
    # code and prevents an additional substitution or passing of a permutation.
    values, exponents = zip(*sorted(decoded_factors, key=lambda k: -k[1]))
    n, fd = Omega_higher(a, exponents)
    rules = dict(zip(n.parent().gens(), values))
    return n.subs(rules), tuple(f.subs(rules) for f in fd)
