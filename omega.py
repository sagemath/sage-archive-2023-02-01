
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
from six import iteritems, itervalues

import operator
from sage.misc.cachefunc import cached_function


def partition(items, predicate=bool):
    r"""
    Split ``items`` into two parts by the given ``predicate``.

    INPUT:

    - ``item`` -- an iterator.

    OUTPUT:

    A pair of iterators; the first contains the elements not satisfying
    the ``predicate``, the second the elements satisfying the ``predicate``.

    Source of the code:
    `http://nedbatchelder.com/blog/201306/filter_a_list_into_two_parts.html`_

    EXAMPLES:

        sage: E, O = partition(srange(10), is_odd)
        sage: tuple(E), tuple(O)
        ((0, 2, 4, 6, 8), (1, 3, 5, 7, 9))
    """
    from itertools import tee
    a, b = tee((predicate(item), item) for item in items)
    return ((item for pred, item in a if not pred),
            (item for pred, item in b if pred))


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
    from sage.misc.misc_c import prod

    return sum(prod(xx**pp for xx, pp in zip(x, p))
               for p in IntegerVectors(j, length=len(x)))


def _laurent_polynomial_ring_(n, m):
    from itertools import chain
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing

    if n + m == 0:
        return ZZ, tuple()
    L = LaurentPolynomialRing(ZZ, ', '.join(chain(
        iter('x{}'.format(nn) for nn in range(n)),
        iter('y{}'.format(mm) for mm in range(m)))))
    return L, L.gens()


def Omega_numerator_P(a, x, y):
    import logging
    logger = logging.getLogger(__name__)

    n = len(x)
    if n == 1:
        x0 = x[0]
        result = x0**(-a) + \
            (prod(1 - x0*yy for yy in y) *
             sum(HomogenousSymmetricFunction(j, y) * (1-x0**(j-a))
                 for j in srange(a))
             if a > 0 else 0)
    else:
        Pprev = Omega_numerator_P(a, x[:n-1], y)
        x1 = x[n-1]
        x2 = x[n-2]
        logger.debug('Omega_numerator: P(%s): substituting...', n)
        p1 = Pprev.subs({x2: x1})
        p2 = Pprev
        logger.debug('Omega_numerator: P(%s): preparing...', n)
        dividend = x1 * (1-x2) * prod(1 - x2*yy for yy in y) * p1 - \
                x2 * (1-x1) * prod(1 - x1*yy for yy in y) * p2
        logger.debug('Omega_numerator: P(%s): dividing...', n)
        q, r = dividend.quo_rem(x1 - x2)
        assert r == 0
        result = q
    logger.debug('Omega_numerator: P(%s) has %s terms', n, result.number_of_terms())
    return result


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

        sage: L.<x0, x1, x2, y0, y1, y2> = LaurentPolynomialRing(ZZ)
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
    import logging
    logger = logging.getLogger(__name__)
    logger.info('Omega_numerator: a=%s, n=%s, m=%s', a, n, m)

    from sage.arith.srange import srange
    from sage.misc.misc_c import prod
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing

    XY, xy = _laurent_polynomial_ring_(n, m)
    x = xy[:n]
    y = xy[n:]

    if m == 0:
        return XY(1 - (prod(prod(f) for f in Omega_factors_denominator(n, m)) *
                       sum(HomogenousSymmetricFunction(j, xy)
                           for j in srange(-a))
                       if a < 0 else 0))
    elif n == 0:
        return XY(sum(HomogenousSymmetricFunction(j, xy)
                      for j in srange(a+1)))
    else:
        return XY(Omega_numerator_P(a, x, y))


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
    of the denominator as a tuple of laurent polynomials, where each
    laurent polynomial `z` represents a factor `1 - z`.

    EXAMPLES::

        sage: Omega_higher(0, (1, -2))
        (1, (z0, z0^2*z1))
        sage: Omega_higher(0, (1, -3))
        (1, (z0, z0^3*z1))
        sage: Omega_higher(0, (1, -4))
        (1, (z0, z0^4*z1))

        sage: Omega_higher(0, (2, -1))
        (z0*z1 + 1, (z0, z0*z1^2))
        sage: Omega_higher(0, (3, -1))
        (z0*z1^2 + z0*z1 + 1, (z0, z0*z1^3))
        sage: Omega_higher(0, (4, -1))
        (z0*z1^3 + z0*z1^2 + z0*z1 + 1, (z0, z0*z1^4))

        sage: Omega_higher(0, (1, 1, -2))
        (-z0^2*z1*z2 - z0*z1^2*z2 + z0*z1*z2 + 1, (z0, z1, z0^2*z2, z1^2*z2))
        sage: Omega_higher(0, (2, -1, -1))
        (z0*z1*z2 + z0*z1 + z0*z2 + 1, (z0, z0*z1^2, z0*z2^2))
        sage: Omega_higher(0, (2, 1, -1))
        (-z0*z1*z2^2 - z0*z1*z2 + z0*z2 + 1, (z0, z1, z0*z2^2, z1*z2))
    """
    import logging
    logger = logging.getLogger(__name__)
    logger.info('Omega_higher: a=%s, exponents=%s', a, exponents)

    from sage.misc.functional import cyclotomic_polynomial
    from sage.misc.misc_c import prod
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    from sage.rings.rational_field import QQ

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

    Z = L.change_ring(ZZ)

    def de_power(expression):
        expression = Z(expression)
        for e, var in zip(exponents, L.gens()):
            if abs(e) == 1:
                continue
            expression = subs_power(expression, var, abs(e))
        return expression

    xy_vars = _laurent_polynomial_ring_(n, m)[1]
    x_vars = iter(xy_vars[:n])
    y_vars = iter(xy_vars[n:])
    rules = {next(x_vars) if e > 0 else next(y_vars):
             powers[abs(e)]**j * var
             for e, var in zip(exponents, L.gens()) for j in range(abs(e))}
    logger.debug('Omega_higher: preparing denominator')
    factors_denominator = tuple(de_power(1 - prod(f.subs(rules) for f in factors))
                                for factors in Omega_factors_denominator(x, y))

    logger.debug('Omega_higher: preparing numerator')
    numerator = de_power(Omega_numerator(a, n, m).subs(rules))

    logger.info('Omega_higher: completed')
    return numerator, factors_denominator


def Omega(var, expression, denominator=None, op=operator.ge,
          Factorization_sort=False, Factorization_simplify=True):
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

        sage: L.<mu, x, y, z, w> = LaurentPolynomialRing(ZZ)

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
        sage: Omega(mu, 2*mu, [])
        2
        sage: Omega(mu, 2/mu, [])
        0

    ::

        sage: Omega(mu, Factorization([(1/mu, 1), (1 - x*mu, -1),
        ....:                          (1 - y/mu, -2)], unit=2))
        2*x * (-x + 1)^-1 * (-x*y + 1)^-2
        sage: Omega(mu, Factorization([(mu, -1), (1 - x*mu, -1),
        ....:                          (1 - y/mu, -2)], unit=2))
        2*x * (-x + 1)^-1 * (-x*y + 1)^-2
        sage: Omega(mu, Factorization([(mu, -1), (1 - x, -1)]))
        0
        sage: Omega(mu, Factorization([(2, -1)]))
        1 * 2^-1

    ::

        sage: Omega(mu, 1, [1 - x*mu, 1 - z, 1 - y/mu])
        1 * (-z + 1)^-1 * (-x + 1)^-1 * (-x*y + 1)^-1
    """
    from sage.arith.misc import factor
    from sage.misc.misc_c import prod
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing_univariate
    from sage.structure.factorization import Factorization

    if op != operator.ge:
        raise NotImplementedError('At the moment, only Omega_ge is implemented.')

    if denominator is None:
        if isinstance(expression, Factorization):
            numerator = expression.unit() * \
                        prod(f**e for f, e in expression if e > 0)
            denominator = tuple(f for f, e in expression if e < 0
                                for _ in range(-e))
        else:
            numerator = expression.numerator()
            denominator = expression.denominator()
    else:
        numerator = expression
    # at this point we have numerator/denominator

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
    # at this point we have numerator/factors_denominator

    P = var.parent()
    if isinstance(P, LaurentPolynomialRing_univariate) and P.gen() == var:
        L = P
        L0 = L.base_ring()
    elif var in P.gens():
        var = repr(var)
        L0 = LaurentPolynomialRing(
            P.base_ring(), tuple(v for v in P.variable_names() if v != var))
        L = LaurentPolynomialRing(L0, var)
        var = L.gen()
    else:
        raise ValueError('{} is not a variable.'.format(var))

    other_factors = []
    to_numerator = []
    decoded_factors = []
    for factor in factors_denominator:
        factor = L(factor)
        D = factor.dict()
        if not D:
            raise ZeroDivisionError('Denominator contains a factor 0.')
        elif len(D) == 1:
            exponent, coefficient = next(iteritems(D))
            if exponent == 0:
                other_factors.append(L0(factor))
            else:
                to_numerator.append(factor)
        elif len(D) == 2:
            if D.get(0, 0) != 1:
                raise NotImplementedError('Factor {} is not normalized.'.format(factor))
            D.pop(0)
            exponent, coefficient = next(iteritems(D))
            decoded_factors.append((-coefficient, exponent))
        else:
            raise NotImplementedError('Cannot handle factor {}.'.format(factor))
    numerator = L(numerator) / prod(to_numerator)

    result_numerator, result_factors_denominator = \
        _Omega_(numerator.dict(), decoded_factors)
    if result_numerator == 0:
        return Factorization([], unit=result_numerator)

    return Factorization([(result_numerator, 1)] +
                         list((f, -1) for f in other_factors) +
                         list((1-f, -1) for f in result_factors_denominator),
                         sort=Factorization_sort,
                         simplify=Factorization_simplify)


def _Omega_(A, decoded_factors):
    r"""
    Helper function for :func:`Omega` which accesses the low level functions
    and does the substituting.

    INPUT:

    - ``A`` -- a dictionary mapping `a` to `c` representing a summand
      `c\lambda^a` of the numerator.

    - ``decoded_factors`` -- a tuple or list of pairs `(z, e)` representing
      a factor `1 - z \lambda^e`.

    OUTPUT:

    A pair representing a quotient as follows: Its first component is the
    numerator as a laurent polynomial, its second component a factorization
    of the denominator as a tuple of laurent polynomials, where each
    laurent polynomial `z` represents a factor `1 - z`.

    TESTS:

    Extensive testing of this function is done in :func:`Omega`.

    ::

        sage: _Omega_({0: 2, 1: 40, -1: -3}, [])
        (42, ())
        sage: _Omega_({-1: 42}, [])
        (0, ())

    """
    if not decoded_factors:
        return sum(c for a, c in iteritems(A) if a >= 0), tuple()

    # Below we sort to make the caching more efficient. Doing this here
    # (in contrast to directly in Omega_higher) results in much cleaner
    # code and prevents an additional substitution or passing of a permutation.
    values, exponents = zip(*sorted(decoded_factors, key=lambda k: -k[1]))

    numerator = 0
    factors_denominator = None
    rules = None
    for a, c in iteritems(A):
        n, fd = Omega_higher(a, exponents)
        if factors_denominator is None:
            factors_denominator = fd
        else:
            assert factors_denominator == fd
        if rules is None:
            rules = dict(zip(n.parent().gens(), values))
        numerator += c * n.subs(rules)

    if numerator == 0:
        factors_denominator = tuple()
    return numerator, tuple(f.subs(rules) for f in factors_denominator)


def generating_function_of_polyhedron(polyhedron, indices=None):
    r"""
    Return the generating function of the integer points of
    the polyhedron's orthant with only nonnegative coordinates.

    EXAMPLES::

        sage: P2 = (
        ....:   Polyhedron(ieqs=[(0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, -1)]),
        ....:   Polyhedron(ieqs=[(0, -1, 0, 1), (0, 1, 0, 0), (0, 0, 1, 0)]))
        sage: generating_function_of_polyhedron(P2[0])
        1 * (-y1 + 1)^-1 * (-y0 + 1)^-1 * (-y0*y2 + 1)^-1
        sage: generating_function_of_polyhedron(P2[1])
        1 * (-y1 + 1)^-1 * (-y2 + 1)^-1 * (-y0*y2 + 1)^-1
        sage: generating_function_of_polyhedron(P2[0] & P2[1])
        1 * (-y1 + 1)^-1 * (-y0*y2 + 1)^-1

    ::

        sage: P3 = (
        ....:   Polyhedron(
        ....:     ieqs=[(0, 0, 0, 0, 1), (0, 0, 0, 1, 0),
        ....:           (0, 0, 1, 0, -1), (-1, 1, 0, -1, -1)]),
        ....:   Polyhedron(
        ....:     ieqs=[(0, 0, -1, 0, 1), (0, 1, 0, 0, -1),
        ....:           (0, 0, 0, 1, 0), (0, 0, 1, 0, 0), (-1, 1, -1, -1, 0)]),
        ....:   Polyhedron(
        ....:     ieqs=[(1, -1, 0, 1, 1), (1, -1, 1, 1, 0),
        ....:           (0, 0, 0, 0, 1), (0, 0, 0, 1, 0), (0, 0, 1, 0, 0),
        ....:           (1, 0, 1, 1, -1), (0, 1, 0, 0, 0), (1, 1, 1, 0, -1)]),
        ....:   Polyhedron(
        ....:     ieqs=[(0, 1, 0, -1, 0), (0, -1, 0, 0, 1),
        ....:           (-1, 0, -1, -1, 1), (0, 0, 1, 0, 0), (0, 0, 0, 1, 0)]),
        ....:   Polyhedron(
        ....:     ieqs=[(0, 1, 0, 0, 0), (0, 0, 1, 0, 0),
        ....:           (-1, -1, -1, 0, 1), (0, -1, 0, 1, 0)]))
        sage: def intersect(I):
        ....:     I = iter(I)
        ....:     result = next(I)
        ....:     for i in I:
        ....:         result &= i
        ....:     return result
        sage: for J in subsets(range(len(P3))):  # TODO: check more results
        ....:     if not J:
        ....:         continue
        ....:     P = intersect([P3[j] for j in J])
        ....:     print('{}: {}'.format(J, P.inequalities()))
        ....:     print(generating_function_of_polyhedron(P))
        [0]: (An inequality (0, 0, 0, 1) x + 0 >= 0,
              An inequality (0, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, -1) x + 0 >= 0,
              An inequality (1, 0, -1, -1) x - 1 >= 0)
        y0 * (-y0 + 1)^-1 * (-y0*y2 + 1)^-1 * (-y1 + 1)^-1 * (-y0*y1*y3 + 1)^-1
        [1]: (An inequality (0, -1, 0, 1) x + 0 >= 0,
              An inequality (0, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, 0) x + 0 >= 0,
              An inequality (1, -1, -1, 0) x - 1 >= 0,
              An inequality (1, 0, 0, -1) x + 0 >= 0)
        (-y0^6*y1^2*y2^2*y3^3 - y0^6*y1^2*y2*y3^3 + y0^5*y1^2*y2*y3^3 +
         y0^5*y1^2*y2*y3^2 + y0^4*y1*y2^2*y3^2 + 2*y0^4*y1*y2*y3^2 +
         y0^4*y1*y3^2 - y0^3*y1*y2*y3^2 - y0^3*y1*y2*y3 - y0^3*y1*y3^2 -
         y0^3*y1*y3 - y0^2*y2*y3 - y0^2*y3 + y0*y3 + y0) *
        (-y0 + 1)^-1 * (-y0*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1 *
        (-y0*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1 *
        (-y0^2*y1*y3 + 1)^-1 * (-y0^2*y1*y2*y3 + 1)^-1
        [0, 1]: (An inequality (1, -1, -1, 0) x - 1 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0,
                 An inequality (0, 0, 1, 0) x + 0 >= 0)
        y0 * (-y0 + 1)^-1 * (-y0*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1
        [2]: (An inequality (-1, 0, 1, 1) x + 1 >= 0,
              An inequality (-1, 1, 1, 0) x + 1 >= 0,
              An inequality (0, 0, 0, 1) x + 0 >= 0,
              An inequality (0, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, 0) x + 0 >= 0,
              An inequality (0, 1, 1, -1) x + 1 >= 0,
              An inequality (1, 0, 0, 0) x + 0 >= 0,
              An inequality (1, 1, 0, -1) x + 1 >= 0)
        (y0^7*y1^6*y2^3*y3^5 + y0^7*y1^5*y2^4*y3^4 + y0^6*y1^7*y2^2*y3^5 -
         y0^7*y1^5*y2^3*y3^4 + y0^6*y1^6*y2^3*y3^4 - y0^6*y1^6*y2^2*y3^5 -
         2*y0^6*y1^6*y2^2*y3^4 - 3*y0^6*y1^5*y2^3*y3^4 - y0^6*y1^5*y2^2*y3^5 -
         y0^6*y1^5*y2^3*y3^3 - y0^6*y1^4*y2^4*y3^3 + y0^6*y1^5*y2^2*y3^4 -
         2*y0^5*y1^6*y2^2*y3^4 - 2*y0^6*y1^4*y2^3*y3^4 - y0^5*y1^6*y2*y3^5 +
         y0^6*y1^5*y2^2*y3^3 + y0^6*y1^4*y2^3*y3^3 - y0^5*y1^5*y2^3*y3^3 -
         y0^6*y1^3*y2^4*y3^3 + y0^6*y1^4*y2^2*y3^4 - y0^5*y1^5*y2^2*y3^4 +
         y0^5*y1^5*y2*y3^5 + 3*y0^5*y1^5*y2^2*y3^3 + y0^6*y1^3*y2^3*y3^3 +
         2*y0^5*y1^5*y2*y3^4 - y0^4*y1^6*y2*y3^4 + 4*y0^5*y1^4*y2^2*y3^4 +
         y0^5*y1^4*y2^3*y3^2 + 3*y0^5*y1^4*y2^2*y3^3 + 4*y0^5*y1^3*y2^3*y3^3 -
         y0^5*y1^4*y2*y3^4 + 3*y0^4*y1^5*y2*y3^4 + y0^5*y1^3*y2^2*y3^4 -
         y0^5*y1^4*y2^2*y3^2 + y0^5*y1^3*y2^3*y3^2 + y0^5*y1^2*y2^4*y3^2 -
         y0^5*y1^4*y2*y3^3 + 2*y0^4*y1^5*y2*y3^3 - 2*y0^5*y1^3*y2^2*y3^3 +
         5*y0^4*y1^4*y2^2*y3^3 + y0^5*y1^2*y2^3*y3^3 - y0^5*y1^3*y2^2*y3^2 -
         y0^5*y1^2*y2^3*y3^2 + 2*y0^4*y1^3*y2^3*y3^2 - 4*y0^4*y1^4*y2*y3^3 +
         2*y0^3*y1^5*y2*y3^3 - y0^5*y1^2*y2^2*y3^3 - y0^4*y1^3*y2^2*y3^3 +
         y0^3*y1^5*y3^4 - y0^4*y1^3*y2*y3^4 - y0^4*y1^4*y2*y3^2 -
         5*y0^4*y1^3*y2^2*y3^2 + y0^3*y1^4*y2^2*y3^2 - y0^4*y1^2*y2^3*y3^2 -
         2*y0^4*y1^3*y2*y3^3 - y0^3*y1^4*y2*y3^3 - 3*y0^4*y1^2*y2^2*y3^3 -
         y0^3*y1^4*y3^4 - y0^4*y1^2*y2^3*y3 + y0^4*y1^3*y2*y3^2 -
         3*y0^3*y1^4*y2*y3^2 - y0^4*y1^2*y2^2*y3^2 - 2*y0^3*y1^3*y2^2*y3^2 -
         y0^4*y1*y2^3*y3^2 - 2*y0^3*y1^4*y3^3 + y0^4*y1^2*y2*y3^3 -
         5*y0^3*y1^3*y2*y3^3 + y0^4*y1^2*y2^2*y3 - y0^3*y1^3*y2^2*y3 +
         y0^4*y1^2*y2*y3^2 - y0^3*y1^3*y2*y3^2 - y0^2*y1^4*y2*y3^2 +
         y0^4*y1*y2^2*y3^2 - 4*y0^3*y1^2*y2^2*y3^2 + y0^3*y1^3*y3^3 -
         2*y0^2*y1^4*y3^3 + y0^3*y1^2*y2*y3^3 + y0^3*y1^3*y2*y3 -
         y0^3*y1*y2^3*y3 + y0^3*y1^3*y3^2 + 5*y0^3*y1^2*y2*y3^2 -
         2*y0^2*y1^3*y2*y3^2 + y0^3*y1*y2^2*y3^2 + y0^2*y1^3*y3^3 +
         y0^3*y1^2*y2*y3 + y0^2*y1^3*y2*y3 + 2*y0^3*y1*y2^2*y3 -
         y0^2*y1^2*y2^2*y3 + 3*y0^2*y1^3*y3^2 + 4*y0^2*y1^2*y2*y3^2 +
         y0^2*y1^2*y3^3 - y0^3*y1*y2*y3 + 4*y0^2*y1^2*y2*y3 +
         2*y0^2*y1*y2^2*y3 + y0^2*y1^2*y3^2 + y0*y1^3*y3^2 +
         2*y0^2*y1*y2*y3^2 + y0^2*y1*y2^2 - y0^2*y1^2*y3 - y0^2*y1*y2*y3 +
         y0*y1^2*y2*y3 + y0^2*y2^2*y3 - y0^2*y1*y3^2 + y0*y1^2*y3^2 -
         y0^2*y1*y2 - y0^2*y1*y3 - y0*y1^2*y3 - y0^2*y2*y3 - 2*y0*y1*y3^2 -
         y0*y1*y2 - 3*y0*y1*y3 - 2*y0*y2*y3 -
         y0*y2 + y0*y3 - y1*y3 + y0 + y3 + 1) *
         (-y0*y1*y3 + 1)^-1 * (-y1 + 1)^-1 * (-y0*y1*y3 + 1)^-1 *
         (-y0*y2 + 1)^-1 * (-y0*y2*y3 + 1)^-1 * (-y1*y3 + 1)^-1 *
         (-y2 + 1)^-1 * (-y0^2*y1*y2*y3 + 1)^-1 *
         (-y0*y1^2*y3 + 1)^-1 * (-y0*y1*y2 + 1)^-1
        [0, 2]: (An inequality (-1, 1, 1, 0) x + 1 >= 0,
                 An inequality (1, 0, -1, 0) x - 1 >= 0,
                 An inequality (0, 0, 1, 0) x + 0 >= 0)
        y0 * (-y0*y2 + 1)^-1 * (-y1 + 1)^-1 * (-y0*y1*y3 + 1)^-1
        [1, 2]: (An inequality (0, -1, 0, 1) x + 0 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0,
                 An inequality (1, 0, 0, -1) x + 0 >= 0,
                 An inequality (1, -1, 0, 0) x - 1 >= 0)
        (y0^4*y1*y2^2*y3^2 - y0^3*y1*y2*y3^2 - y0^3*y1*y2*y3 -
         y0^2*y2*y3 + y0*y3 + y0) *
        (-y0*y1*y3 + 1)^-1 * (-y0*y2 + 1)^-1 *
        (-y0*y2*y3 + 1)^-1 * (-y0^2*y1*y2*y3 + 1)^-1
        [0, 1, 2]: (An inequality (0, 1, 0, 0) x + 0 >= 0,
                    An inequality (1, -1, 0, 0) x - 1 >= 0)
        y0 * (-y0*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1
        [3]: (An inequality (-1, 0, 0, 1) x + 0 >= 0,
              An inequality (0, -1, -1, 1) x - 1 >= 0,
              An inequality (0, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, 0) x + 0 >= 0,
              An inequality (1, 0, -1, 0) x + 0 >= 0)
        (-y0*y1*y3^2 - y0*y3^2 + y0*y3 + y3) *
        (-y0*y2*y3 + 1)^-1 * (-y3 + 1)^-1 * (-y1*y3 + 1)^-1 *
        (-y0*y3 + 1)^-1 * (-y0*y1*y3 + 1)^-1
        [0, 3]: ()
        0
        [1, 3]: (An inequality (1, -1, -1, 0) x - 1 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0,
                 An inequality (0, 0, 1, 0) x + 0 >= 0)
        y0*y3 * (-y0*y3 + 1)^-1 * (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 1, 3]: ()
        0
        [2, 3]: (An inequality (1, 0, -1, 0) x + 0 >= 0,
                 An inequality (-1, 1, 1, 0) x + 1 >= 0,
                 An inequality (0, 0, 1, 0) x + 0 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0)
        (y0^2*y1^2*y2*y3^4 - y0^2*y1*y2*y3^3 - y0*y1*y2*y3^3 -
         y0*y1*y3^2 + y0*y3 + y3) *
        (-y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1 *
        (-y0*y1*y3 + 1)^-1 * (-y0*y1*y2*y3^2 + 1)^-1
        [0, 2, 3]: ()
        0
        [1, 2, 3]: (An inequality (0, 1, 0, 0) x + 0 >= 0,
                    An inequality (1, -1, 0, 0) x - 1 >= 0)
        y0*y3 * (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 1, 2, 3]: ()
        0
        [4]: (An inequality (-1, -1, 0, 1) x - 1 >= 0,
              An inequality (-1, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, 0) x + 0 >= 0,
              An inequality (1, 0, 0, 0) x + 0 >= 0)
        y3 * (-y2 + 1)^-1 * (-y3 + 1)^-1 * (-y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 4]: ()
        0
        [1, 4]: ()
        0
        [0, 1, 4]: ()
        0
        [2, 4]: (An inequality (-1, 0, 1, 0) x + 0 >= 0,
                 An inequality (1, 0, 0, 0) x + 0 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0)
        y3 * (-y2 + 1)^-1 * (-y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 2, 4]: ()
        0
        [1, 2, 4]: ()
        0
        [0, 1, 2, 4]: ()
        0
        [3, 4]: (An inequality (0, 1, 0, 0) x + 0 >= 0,
                 An inequality (-1, -1, 0, 1) x - 1 >= 0,
                 An inequality (1, 0, 0, 0) x + 0 >= 0)
        y3 * (-y3 + 1)^-1 * (-y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 3, 4]: ()
        0
        [1, 3, 4]: ()
        0
        [0, 1, 3, 4]: ()
        0
        [2, 3, 4]: (An inequality (0, 1, 0, 0) x + 0 >= 0,
                    An inequality (1, 0, 0, 0) x + 0 >= 0)
        y3 * (-y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 2, 3, 4]: ()
        0
        [1, 2, 3, 4]: ()
        0
        [0, 1, 2, 3, 4]: ()
        0
    """
    from sage.geometry.polyhedron.representation import Hrepresentation
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    from sage.structure.factorization import Factorization

    try:
        if polyhedron.is_empty():
            return Factorization([], unit=0)
    except AttributeError:
        pass

    def inequalities_coeffs(inequalities):
        for entry in inequalities:
            if isinstance(entry, (tuple, list)):
                yield tuple(entry)
            elif isinstance(entry, Hrepresentation):
                if entry.is_inequality():
                    yield tuple(entry.vector())
                elif entry.is_equation():
                    e = tuple(entry.vector())
                    yield e
                    yield tuple(-ee for ee in e)
                else:
                    raise ValueError(
                        'Cannot handle Hrepresentation {}.'.format(entry))
            else:
                raise ValueError('Cannot handle {}.'.format(entry))

    try:
        inequalities = polyhedron.Hrepresentation()
    except AttributeError:
        inequalities = polyhedron
    inequalities = tuple(inequalities_coeffs(inequalities))
    if not inequalities:
        raise ValueError('no inequality given')

    if indices is None:
        indices = range(len(inequalities[0]) - 1)
    B = LaurentPolynomialRing(
        ZZ,
        ', '.join('y{}'.format(k) for k in indices),
        sparse=True)

    n = len(B.gens()) + 1
    if any(len(ineq) != n for ineq in inequalities):
        raise ValueError('Not all coefficient vectors of the inequalities '
                         'have the same length.')

    def is_unit_vector(it):
        found = 0
        for e in it:
            if e != 0:
                if e != 1:
                    return False
                else:
                    found += 1
                    if found >= 2:
                        return False
        return True

    numerator = B(1)
    terms = B.gens()
    L = B
    for i, coeffs in enumerate(inequalities):
        if is_unit_vector(coeffs):
            continue
        L = LaurentPolynomialRing(L, 'lambda{}'.format(i), sparse=True)
        l = L.gen()
        it_coeffs = iter(coeffs)
        numerator *= l**next(it_coeffs)
        assert numerator.parent() == L
        terms = tuple(l**c * t for c, t in zip(it_coeffs, terms))

    def decode_factor(factor):
        D = factor.dict()
        assert len(D) == 1
        exponent, coefficient = next(iteritems(D))
        return coefficient, exponent

    while repr(numerator.parent().gen()).startswith('lambda'):
        decoded_factors, other_factors = \
            partition((decode_factor(factor) for factor in terms),
                      lambda factor: factor[1] == 0)
        other_factors = tuple(factor[0] for factor in other_factors)
        numerator, factors_denominator = \
            _Omega_(numerator.dict(), tuple(decoded_factors))
        terms = other_factors + factors_denominator

    return Factorization([(numerator, 1)] +
                         list((1-t, -1) for t in terms),
                         sort=False, simplify=False)
