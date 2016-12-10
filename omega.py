
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


def Omega_numerator_P(a, x, y, t):
    import logging
    logger = logging.getLogger(__name__)

    from sage.arith.srange import srange
    from sage.misc.misc_c import prod

    n = len(x)
    if n == 1:
        x0 = t
        result = x0**(-a) + \
            (prod(1 - x0*yy for yy in y) *
             sum(HomogenousSymmetricFunction(j, y) * (1-x0**(j-a))
                 for j in srange(a))
             if a > 0 else 0)
    else:
        Pprev = Omega_numerator_P(a, x[:n-1], y, t)
        x2 = x[n-2]
        logger.debug('Omega_numerator: P(%s): substituting...', n)
        x1 = t
        p1 = Pprev
        p2 = Pprev.subs({t: x2})
        logger.debug('Omega_numerator: P(%s): preparing...', n)
        dividend = x1 * (1-x2) * prod(1 - x2*yy for yy in y) * p1 - \
                x2 * (1-x1) * prod(1 - x1*yy for yy in y) * p2
        logger.debug('Omega_numerator: P(%s): dividing...', n)
        q, r = dividend.quo_rem(x1 - x2)
        assert r == 0
        result = q
    logger.debug('Omega_numerator: P(%s) has %s terms', n, result.number_of_terms())
    return result


def Omega_numerator(a, x, y, t):
    r"""
    Return the numerator of `\Omega_{\ge}` of the expression
    specified by the input.

    To be more precise, this calculates

    .. MATH::

        \Omega_{\ge} \frac{\mu^a}{
        (1 - x_1 \mu) \dots (1 - x_n \mu)
        (1 - y_1 / \mu) \dots (1 - y_m / \mu)

    and returns its numerator.

    INPUT:

    - ``a`` -- an integer.

    - ``x`` and ``y`` -- a tuple of tuples of laurent polynomials. The
      flattened ``x`` contains `x_1,...,x_n`, the flattened ``y`` the
      `y_1,...,y_m`.

    - ``t`` -- a temporary laurent polynomial variable used for substituting.

    OUTPUT:

    A laurent polynomial.

    EXAMPLES::

        sage: L.<x0, x1, x2, x3, y0, y1, t> = LaurentPolynomialRing(ZZ)
        sage: Omega_numerator(0, ((x0,),), ((y0,),), t)
        1
        sage: Omega_numerator(0, ((x0,), (x1,)), ((y0,),), t)
        -x0*x1*y0 + 1
        sage: Omega_numerator(0, ((x0,),), ((y0,), (y1,)), t)
        1
        sage: Omega_numerator(0, ((x0,), (x1,), (x2,)), ((y0,),), t)
        x0*x1*x2*y0^2 + x0*x1*x2*y0 - x0*x1*y0 - x0*x2*y0 - x1*x2*y0 + 1
        sage: Omega_numerator(0, ((x0,), (x1,)), ((y0,), (y1,)), t)
        x0^2*x1*y0*y1 + x0*x1^2*y0*y1 - x0*x1*y0*y1 - x0*x1*y0 - x0*x1*y1 + 1

        sage: Omega_numerator(-2, ((x0,),), ((y0,),), t)
        x0^2
        sage: Omega_numerator(-1, ((x0,),), ((y0,),), t)
        x0
        sage: Omega_numerator(1, ((x0,),), ((y0,),), t)
        -x0*y0 + y0 + 1
        sage: Omega_numerator(2, ((x0,),), ((y0,),), t)
        -x0*y0^2 - x0*y0 + y0^2 + y0 + 1

    TESTS::

        sage: Omega_factors_denominator((), ())
        ()
        sage: Omega_numerator(0, (), (), t)
        1
        sage: Omega_numerator(+2, (), (), t)
        1
        sage: Omega_numerator(-2, (), (), t)
        0

        sage: Omega_factors_denominator(((x0,),), ())
        (-x0 + 1,)
        sage: Omega_numerator(0, ((x0,),), (), t)
        1
        sage: Omega_numerator(+2, ((x0,),), (), t)
        1
        sage: Omega_numerator(-2, ((x0,),), (), t)
        x0^2

        sage: Omega_factors_denominator((), ((y0,),))
        ()
        sage: Omega_numerator(0, (), ((y0,),), t)
        1
        sage: Omega_numerator(+2, (), ((y0,),), t)
        y0^2 + y0 + 1
        sage: Omega_numerator(-2, (), ((y0,),), t)
        0

    ::

        sage: L.<X, Y, t> = LaurentPolynomialRing(ZZ)
        sage: Omega_numerator(2, ((X,),), ((Y,),), t)
        -X*Y^2 - X*Y + Y^2 + Y + 1
    """
    from sage.arith.srange import srange
    from sage.misc.misc_c import prod

    x_flat = sum(x, tuple())
    y_flat = sum(y, tuple())
    n = len(x_flat)
    m = len(y_flat)
    xy = x_flat + y_flat

    import logging
    logger = logging.getLogger(__name__)
    logger.info('Omega_numerator: a=%s, n=%s, m=%s', a, n, m)

    if m == 0:
        result = 1 - (prod(Omega_factors_denominator(x, y)) *
                      sum(HomogenousSymmetricFunction(j, xy)
                          for j in srange(-a))
                      if a < 0 else 0)
    elif n == 0:
        result = sum(HomogenousSymmetricFunction(j, xy)
                     for j in srange(a+1))
    else:
        result = Omega_numerator_P(a, x_flat, y_flat, t).subs({t: x_flat[-1]})
    L = t.parent()
    result = L(result)

    logger.info('Omega_numerator: %s terms', result.number_of_terms())
    return result


@cached_function
def Omega_factors_denominator(x, y):
    r"""
    Return the denominator of `\Omega_{\ge}` of the expression
    specified by the input.

    To be more precise, this calculates

    .. MATH::

        \Omega_{\ge} \frac{1}{
        (1 - x_1 \mu) \dots (1 - x_n \mu)
        (1 - y_1 / \mu) \dots (1 - y_m / \mu)

    and returns a factorization of its denominator.

    INPUT:

    - ``x`` and ``y`` -- a tuple of tuples of laurent polynomials. The
      flattened ``x`` contains `x_1,...,x_n`, the flattened ``y`` the
      `y_1,...,y_m`.

    OUTPUT:

    A factorization of the denominator as
    a tuple of laurent polynomials.

    EXAMPLES::

        sage: L.<x0, x1, x2, x3, y0, y1> = LaurentPolynomialRing(ZZ)
        sage: Omega_factors_denominator(((x0,),), ((y0,),))
        (-x0 + 1, -x0*y0 + 1)
        sage: Omega_factors_denominator(((x0,),), ((y0,), (y1,)))
        (-x0 + 1, -x0*y0 + 1, -x0*y1 + 1)
        sage: Omega_factors_denominator(((x0,), (x1,)), ((y0,),))
        (-x0 + 1, -x1 + 1, -x0*y0 + 1, -x1*y0 + 1)
        sage: Omega_factors_denominator(((x0,), (x1,), (x2,)), ((y0,),))
        (-x0 + 1, -x1 + 1, -x2 + 1, -x0*y0 + 1, -x1*y0 + 1, -x2*y0 + 1)
        sage: Omega_factors_denominator(((x0,), (x1,)), ((y0,), (y1,)))
        (-x0 + 1, -x1 + 1, -x0*y0 + 1, -x0*y1 + 1, -x1*y0 + 1, -x1*y1 + 1)

    ::

        sage: B.<zeta> = ZZ.extension(cyclotomic_polynomial(3))
        sage: L.<x, y> = LaurentPolynomialRing(B)
        sage: Omega_factors_denominator(((x, -x),), ((y,),))
        (-x^2 + 1, -x^2*y^2 + 1)
        sage: Omega_factors_denominator(((x, -x),), ((y, zeta*y, zeta^2*y),))
        (-x^2 + 1, -x^6*y^6 + 1)
        sage: Omega_factors_denominator(((x, -x),), ((y, -y),))
        (-x^2 + 1, -x^2*y^2 + 1, -x^2*y^2 + 1)

    TESTS::

        sage: L.<x0, y0> = LaurentPolynomialRing(ZZ)
        sage: Omega_factors_denominator((), ())
        ()
        sage: Omega_factors_denominator(((x0,),), ())
        (-x0 + 1,)
        sage: Omega_factors_denominator((), ((y0,),))
        ()
    """
    import logging
    logger = logging.getLogger(__name__)

    from sage.misc.misc_c import prod

    result = tuple(prod(1 - xx for xx in gx) for gx in x) + \
             sum(((prod(1 - xx*yy for xx in gx for yy in gy),)
                  if len(gx) != len(gy)
                  else tuple(prod(1 - xx*yy for xx in gx) for yy in gy)
                  for gx in x for gy in y),
                 tuple())

    logger.info('Omega_denominator: %s factors', len(result))
    return result


@cached_function
def Omega_higher(a, exponents):
    r"""
    Return `\Omega_{\ge}` of the expression specified by the input.

    To be more precise, this calculates

    .. MATH::

        \Omega_{\ge} \frac{\mu^a}{
        (1 - z_1 \mu^{e_1}) \dots (1 - z_n \mu^{e_n})

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

    ::

        sage: Omega_higher(0, (2, -2))
        (-z0*z1 + 1, (z0, z0*z1, z0*z1))
        sage: Omega_higher(0, (2, -3))
        (z0^2*z1 + 1, (z0, z0^3*z1^2))
        sage: Omega_higher(0, (3, 1, -3))
        (-z0^3*z1^3*z2^3 + 2*z0^2*z1^3*z2^2 - z0*z1^3*z2
         + z0^2*z2^2 - 2*z0*z2 + 1,
         (z0, z1, z0*z2, z0*z2, z0*z2, z1^3*z2))

    ::

        sage: Omega_higher(0, (3, 6, -1))
        (-z0*z1*z2^8 - z0*z1*z2^7 - z0*z1*z2^6 - z0*z1*z2^5 - z0*z1*z2^4 +
         z1*z2^5 - z0*z1*z2^3 + z1*z2^4 - z0*z1*z2^2 + z1*z2^3 -
         z0*z1*z2 + z0*z2^2 + z1*z2^2 + z0*z2 + z1*z2 + 1,
         (z0, z1, z0*z2^3, z1*z2^6))

    ::

        sage: Omega_higher(0, (2, 2, 1, 1, 1, 1, 1, -1, -1))[0].number_of_terms()  # long time
        27837
    """
    import logging
    logger = logging.getLogger(__name__)
    logger.info('Omega_higher: a=%s, exponents=%s', a, exponents)

    from sage.arith.misc import lcm
    from sage.arith.srange import srange
    from sage.misc.functional import cyclotomic_polynomial
    from sage.misc.misc_c import prod
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    from sage.rings.rational_field import QQ

    if not exponents or any(e == 0 for e in exponents):
        raise NotImplementedError

    rou = sorted(set(abs(e) for e in exponents) - set([1]))
    ellcm = lcm(rou)
    B = QQ.extension(cyclotomic_polynomial(ellcm), 'zeta')
    zeta = B.gen()
    z_names = tuple('z{}'.format(i) for i in range(len(exponents)))
    L = LaurentPolynomialRing(B, ('t',) + z_names)
    t = L.gens()[0]
    Z = LaurentPolynomialRing(ZZ, z_names)
    powers = {i: L(zeta**(ellcm//i)) for i in rou}
    powers[2] = L(-1)
    powers[1] = L(1)
    exponents_and_values = tuple(
        (e, tuple(powers[abs(e)]**j * z for j in srange(abs(e))))
        for z, e in zip(L.gens()[1:], exponents))
    x = tuple(v for e, v in exponents_and_values if e > 0)
    y = tuple(v for e, v in exponents_and_values if e < 0)

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

    def de_power(expression):
        expression = Z(expression)
        for e, var in zip(exponents, Z.gens()):
            if abs(e) == 1:
                continue
            expression = subs_power(expression, var, abs(e))
        return expression

    logger.debug('Omega_higher: preparing denominator')
    factors_denominator = tuple(de_power(1 - factor)
                                for factor in Omega_factors_denominator(x, y))

    logger.debug('Omega_higher: preparing numerator')
    numerator = de_power(Omega_numerator(a, x, y, t))

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
      `c\mu^a` of the numerator.

    - ``decoded_factors`` -- a tuple or list of pairs `(z, e)` representing
      a factor `1 - z \mu^e`.

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


def pretty_inequality(ineq, indices=None):
    from sage.symbolic.ring import SR
    from sage.modules.free_module_element import vector

    if indices is None:
        indices = range(len(ineq)-1)
    vars = vector([1] + list(SR("b{}".format(i)) for i in indices))
    v = vector(ineq)
    positive_part = vector([max(c, 0) for c in v])
    negative_part = - (v - positive_part)
    assert v == positive_part - negative_part
    return '{} >= {}'.format(positive_part*vars, negative_part*vars)


def prepare_inequalities(inequalities, B):
    r"""
    EXAMPLES::

        sage: B = LaurentPolynomialRing(ZZ, 'y', 3)
        sage: prepare_inequalities([(0, -1, 1, 0), (2, -1, -1, 1)], B)
        ([(2, -2, -1, 1)], 1, {y2: y2, y1: y1, y0: y0*y1})
        sage: prepare_inequalities([(-1, -1, 1, 0), (2, -1, -1, 1)], B)
        ([(1, -2, -1, 1)], y1, {y2: y2, y1: y1, y0: y0*y1})
        sage: prepare_inequalities([(2, -1, 1, 0), (2, -1, -1, 1)], B)
        ([(4, -2, -1, 1)], y1^-2, {y2: y2, y1: y1, y0: y0*y1})

    TESTS::

        sage: B = LaurentPolynomialRing(ZZ, 'y', 3)
        sage: prepare_inequalities([(1, 1, -1, 0), (1, -1, 0, 1),
        ....:                       (2, -1, -1, 3)], B)
        ([(-3, 2, 1, 3)], y0^-1*y2^-2, {y2: y2, y1: y0*y1*y2, y0: y0*y2})

        sage: B = LaurentPolynomialRing(ZZ, 'y', 4)
        sage: prepare_inequalities([(-1, 1, -1, 0, 0)], B)
        ([], y0, {y3: y3, y2: y2, y1: y0*y1, y0: y0})
        sage: prepare_inequalities([(0, 0, -1, 0, 1), (0, 0, 1, 0, 0),
        ....:                       (0, 1, 0, 0, -1), (-1, 1, -1, 0, 0)], B)
        ([(1, 1, 0, 0, -1)], y0, {y3: y3, y2: y2, y1: y0*y1*y3, y0: y0})
        sage: prepare_inequalities([(-2, 1, -1, 0, 0)], B)
        ([], y0^2, {y3: y3, y2: y2, y1: y0*y1, y0: y0})

        sage: prepare_inequalities([(0, -1, 1, 0, 0), (-2, 0, -1, 0, 1),
        ....:                       (0, -1, 0, 1, 0), (-3, 0, 0, -1, 1)], B)
        ([(1, 0, -1, 1, 1)],
         y3^3,
         {y3: y3, y2: y2*y3, y1: y1, y0: y0*y1*y2*y3})
        sage: prepare_inequalities([(0, -1, 1, 0, 0), (-3, 0, -1, 0, 1),
        ....:                       (0, -1, 0, 1, 0), (-2, 0, 0, -1, 1)], B)
        ([(1, 0, 1, -1, 1)],
         y3^3,
         {y3: y3, y2: y2, y1: y1*y3, y0: y0*y1*y2*y3})
        sage: prepare_inequalities([(0, -1, 1, 0, 0), (-2, 0, -1, 0, 1),
        ....:                       (-3, -1, 0, 1, 0), (0, 0, 0, -1, 1)], B)
        ([(1, 0, -1, 1, 1)],
         y2^3*y3^3,
         {y3: y3, y2: y2*y3, y1: y1, y0: y0*y1*y2*y3})
        sage: prepare_inequalities([(0, -1, 1, 0, 0), (-3, 0, -1, 0, 1),
        ....:                       (-2, -1, 0, 1, 0), (0, 0, 0, -1, 1)], B)
        ([(1, 0, 1, -1, 1)],
         y2^2*y3^3,
         {y3: y3, y2: y2, y1: y1*y3, y0: y0*y1*y2*y3})
        sage: prepare_inequalities([(-2, -1, 1, 0, 0), (0, 0, -1, 0, 1),
        ....:                       (0, -1, 0, 1, 0), (-3, 0, 0, -1, 1)], B)
        ([(1, 0, -1, 1, 1)],
         y1^2*y3^3,
         {y3: y3, y2: y2*y3, y1: y1, y0: y0*y1*y2*y3})
        sage: prepare_inequalities([(-3, -1, 1, 0, 0), (0, 0, -1, 0, 1),
        ....:                       (0, -1, 0, 1, 0), (-2, 0, 0, -1, 1)], B)
        ([(1, 0, 1, -1, 1)],
         y1^3*y3^3,
         {y3: y3, y2: y2, y1: y1*y3, y0: y0*y1*y2*y3})
        sage: prepare_inequalities([(-2, -1, 1, 0, 0), (0, 0, -1, 0, 1),
        ....:                       (-3, -1, 0, 1, 0), (0, 0, 0, -1, 1)], B)
        ([(1, 0, -1, 1, 1)],
         y1^2*y2^3*y3^3,
         {y3: y3, y2: y2*y3, y1: y1, y0: y0*y1*y2*y3})
        sage: prepare_inequalities([(-3, -1, 1, 0, 0), (0, 0, -1, 0, 1),
        ....:                       (-2, -1, 0, 1, 0), (0, 0, 0, -1, 1)], B)
        ([(1, 0, 1, -1, 1)],
         y1^3*y2^2*y3^3,
         {y3: y3, y2: y2, y1: y1*y3, y0: y0*y1*y2*y3})
    """
    import logging
    logger = logging.getLogger(__name__)

    from itertools import takewhile
    from sage.graphs.digraph import DiGraph
    from sage.matrix.constructor import matrix
    from sage.modules.free_module_element import vector
    from sage.rings.integer_ring import ZZ

    inequalities_filtered = []
    chain_links = {}
    for coeffs in inequalities:
        dim = len(coeffs)
        constant = coeffs[0]
        ones = tuple(i+1 for i, c in enumerate(coeffs[1:]) if c == 1)
        mones = tuple(i+1 for i, c in enumerate(coeffs[1:]) if c == -1)
        absgetwo = tuple(i+1 for i, c in enumerate(coeffs[1:]) if abs(c) >= 2)
        if len(ones) == 1 and not mones and not absgetwo:
            if constant >= 0:
                logger.debug('generating_function_of_polyhedron: '
                             'skipping %s', pretty_inequality(coeffs))
            else:
                # TODO: could be skipped...
                inequalities_filtered.append(coeffs)
        elif len(ones) == 1 and len(mones) == 1 and not absgetwo:
            chain_links[(mones[0], ones[0])] = constant
        else:
            inequalities_filtered.append(coeffs)

    G = DiGraph(chain_links, format='list_of_edges')
    potential = {}
    paths = {}
    D = {}
    inequalities_extra = []
    for i in range(dim):
        D[(i, i)] = 1
    for v in G.topological_sort():
        NP = iter(sorted(((n, potential[n] + chain_links[(n, v)])
                          for n in G.neighbor_in_iterator(v)),
                         key=lambda k: k[1]))
        n, p = next(NP, (None, 0))
        potential[v] = p
        D[(0, v)] = -p
        paths[v] = paths.get(n, ()) + (v,)
        for u in paths[v]:
            D[(u, v)] = 1

        for n, p in NP:
            ell = len(tuple(takewhile(lambda u: u[0] == u[1],
                                      zip(paths[n], paths[v]))))
            coeffs = dim*[0]
            for u in paths[v][ell:]:
                coeffs[u] = 1
            for u in paths[n][ell:]:
                coeffs[u] = -1
            coeffs[0] = p - potential[v]
            inequalities_extra.append(tuple(coeffs))
    T = matrix(ZZ, dim, dim, D)

    inequalities = list(tuple(T*vector(ieq))
                        for ieq in inequalities_filtered) + \
                   inequalities_extra

    rules_pre = iter((y, B({tuple(row[1:]): 1}))
                     for y, row in zip((1,) + B.gens(), T.rows()))
    factor = next(rules_pre)[1]
    rules = dict(rules_pre)

    return inequalities, factor, rules


def prepare_equations(equations, B):
    r"""
    EXAMPLES::

        sage: B = LaurentPolynomialRing(ZZ, 'y', 4)
        sage: prepare_equations([(1, 1, 1, -1, 0)], B)
        (y2, {y1: y1*y2, y0: y0*y2}, (2,))
        sage: prepare_equations([(0, 1, 0, -1, 0)], B)
        (1, {y0: y0*y2}, (2,))
        sage: prepare_equations([(-1, 0, 1, -1, -1), (1, 1, 0, 1, 2)], B)
        (y2^-1, {y1: y1*y2^2*y3^-1, y0: y0*y2*y3^-1}, (2, 3))

    TESTS::

        sage: B = LaurentPolynomialRing(ZZ, 'y', 4)
        sage: prepare_equations([(0, 0, 1, 0, -1), (-1, 1, -1, -1, 0)], B)
        (y2^-1, {y1: y1*y2^-1*y3, y0: y0*y2}, (2, 3))

        sage: B = LaurentPolynomialRing(ZZ, 'y', 5)
        sage: prepare_equations([(0, 0, 0, 1, 0, -1), (0, 1, 0, 0, -1, 0),
        ....:                    (0, 1, -1, 0, 0, 0)], B)
        (1, {y2: y2*y4, y0: y0*y1*y3}, (1, 3, 4))
    """
    import logging
    logger = logging.getLogger(__name__)

    from sage.matrix.constructor import matrix
    from sage.misc.misc_c import prod

    E = matrix(equations)
    if not E:
        return 1, {}, ()

    indices_nonzero = tuple(i for i, col in enumerate(E.columns())
                            if i > 0 and not col.is_zero())
    indices = []
    r = 0
    for i in reversed(indices_nonzero):
        indices.append(i)
        r1 = E.matrix_from_columns(indices).rank()
        if r1 > r:
            r = r1
            if len(indices) >= E.nrows():
                break
        else:
            indices = indices[:-1]
    assert len(indices) == E.nrows()
    indices = tuple(reversed(indices))
    indicesn = (0,) + tuple(i for i in indices_nonzero if i not in indices)
    T = E.matrix_from_columns(indices).inverse()
    TE = T*E

    gens = (1,) + B.gens()
    z = tuple(gens[i] for i in indices)
    gens_cols = zip(gens, TE.columns())
    rules_pre = iter((y, y * prod(zz**(-c) for zz, c in zip(z, col)))
                     for y, col in (gens_cols[i] for i in indicesn))
    factor = next(rules_pre)[1]
    rules = dict(rules_pre)

    return factor, rules, tuple(i-1 for i in indices)



def generating_function_of_polyhedron(polyhedron, split=False,
                                      result_as_tuple=None, **kwds):
    r"""
    Return the generating function of the integer points of
    the polyhedron's orthant with only nonnegative coordinates.

    EXAMPLES::

        sage: P2 = (
        ....:   Polyhedron(ieqs=[(0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, -1)]),
        ....:   Polyhedron(ieqs=[(0, -1, 0, 1), (0, 1, 0, 0), (0, 0, 1, 0)]))
        sage: generating_function_of_polyhedron(P2[0], sort_factors=True)
        1 * (-y0 + 1)^-1 * (-y1 + 1)^-1 * (-y0*y2 + 1)^-1
        sage: generating_function_of_polyhedron(P2[1], sort_factors=True)
        1 * (-y1 + 1)^-1 * (-y2 + 1)^-1 * (-y0*y2 + 1)^-1
        sage: (P2[0] & P2[1]).Hrepresentation()
        (An equation (1, 0, -1) x + 0 == 0,
         An inequality (1, 0, 0) x + 0 >= 0,
         An inequality (0, 1, 0) x + 0 >= 0)
        sage: generating_function_of_polyhedron(P2[0] & P2[1], sort_factors=True)
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
        ....:     print('{}: {}'.format(J, P.Hrepresentation()))
        ....:     print(generating_function_of_polyhedron(P, sort_factors=True))
        [0]: (An inequality (0, 0, 0, 1) x + 0 >= 0,
              An inequality (0, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, -1) x + 0 >= 0,
              An inequality (1, 0, -1, -1) x - 1 >= 0)
        y0 * (-y0 + 1)^-1 * (-y1 + 1)^-1 * (-y0*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1
        [1]: (An inequality (0, -1, 0, 1) x + 0 >= 0,
              An inequality (0, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, 0) x + 0 >= 0,
              An inequality (1, -1, -1, 0) x - 1 >= 0,
              An inequality (1, 0, 0, -1) x + 0 >= 0)
        (-y0^2*y2*y3 - y0^2*y3 + y0*y3 + y0) *
        (-y0 + 1)^-1 * (-y0*y2 + 1)^-1 * (-y0*y3 + 1)^-1 *
        (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 1]: (An equation (0, 1, 0, -1) x + 0 == 0,
                 An inequality (1, -1, -1, 0) x - 1 >= 0,
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
         (-y1 + 1)^-1 * (-y2 + 1)^-1 * (-y0*y2 + 1)^-1 * (-y1*y3 + 1)^-1 *
         (-y0*y1*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1 * (-y0*y1*y3 + 1)^-1 *
         (-y0*y2*y3 + 1)^-1 * (-y0*y1^2*y3 + 1)^-1 * (-y0^2*y1*y2*y3 + 1)^-1
        [0, 2]: (An equation (1, 0, -1, -1) x - 1 == 0,
                 An inequality (-1, 1, 1, 0) x + 1 >= 0,
                 An inequality (1, 0, -1, 0) x - 1 >= 0,
                 An inequality (0, 0, 1, 0) x + 0 >= 0)
        y0 * (-y1 + 1)^-1 * (-y0*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1
        [1, 2]: (An equation (1, -1, -1, 0) x - 1 == 0,
                 An inequality (0, -1, 0, 1) x + 0 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0,
                 An inequality (1, 0, 0, -1) x + 0 >= 0,
                 An inequality (1, -1, 0, 0) x - 1 >= 0)
        (-y0^2*y2*y3 + y0*y3 + y0) *
        (-y0*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 1, 2]: (An equation (0, 1, 0, -1) x + 0 == 0,
                    An equation (1, -1, -1, 0) x - 1 == 0,
                    An inequality (0, 1, 0, 0) x + 0 >= 0,
                    An inequality (1, -1, 0, 0) x - 1 >= 0)
        y0 * (-y0*y2 + 1)^-1 * (-y0*y1*y3 + 1)^-1
        [3]: (An inequality (-1, 0, 0, 1) x + 0 >= 0,
              An inequality (0, -1, -1, 1) x - 1 >= 0,
              An inequality (0, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, 0) x + 0 >= 0,
              An inequality (1, 0, -1, 0) x + 0 >= 0)
        (-y0*y1*y3^2 - y0*y3^2 + y0*y3 + y3) *
        (-y3 + 1)^-1 * (-y0*y3 + 1)^-1 *
        (-y1*y3 + 1)^-1 * (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 3]: (An equation -1 == 0,)
        0
        [1, 3]: (An equation (1, 0, 0, -1) x + 0 == 0,
                 An inequality (1, -1, -1, 0) x - 1 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0,
                 An inequality (0, 0, 1, 0) x + 0 >= 0)
        y0*y3 * (-y0*y3 + 1)^-1 * (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 1, 3]: (An equation -1 == 0,)
        0
        [2, 3]: (An equation (0, 1, 1, -1) x + 1 == 0,
                 An inequality (1, 0, -1, 0) x + 0 >= 0,
                 An inequality (-1, 1, 1, 0) x + 1 >= 0,
                 An inequality (0, 0, 1, 0) x + 0 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0)
        (-y0*y1*y3^2 + y0*y3 + y3) *
        (-y1*y3 + 1)^-1 * (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 2, 3]: (An equation -1 == 0,)
        0
        [1, 2, 3]: (An equation (1, 0, 0, -1) x + 0 == 0,
                    An equation (1, -1, -1, 0) x - 1 == 0,
                    An inequality (0, 1, 0, 0) x + 0 >= 0,
                    An inequality (1, -1, 0, 0) x - 1 >= 0)
        y0*y3 * (-y0*y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 1, 2, 3]: (An equation -1 == 0,)
        0
        [4]: (An inequality (-1, -1, 0, 1) x - 1 >= 0,
              An inequality (-1, 0, 1, 0) x + 0 >= 0,
              An inequality (0, 1, 0, 0) x + 0 >= 0,
              An inequality (1, 0, 0, 0) x + 0 >= 0)
        y3 * (-y2 + 1)^-1 * (-y3 + 1)^-1 * (-y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 4]: (An equation -1 == 0,)
        0
        [1, 4]: (An equation -1 == 0,)
        0
        [0, 1, 4]: (An equation -1 == 0,)
        0
        [2, 4]: (An equation (1, 1, 0, -1) x + 1 == 0,
                 An inequality (-1, 0, 1, 0) x + 0 >= 0,
                 An inequality (1, 0, 0, 0) x + 0 >= 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0)
        y3 * (-y2 + 1)^-1 * (-y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 2, 4]: (An equation -1 == 0,)
        0
        [1, 2, 4]: (An equation -1 == 0,)
        0
        [0, 1, 2, 4]: (An equation -1 == 0,)
        0
        [3, 4]: (An equation (1, 0, -1, 0) x + 0 == 0,
                 An inequality (0, 1, 0, 0) x + 0 >= 0,
                 An inequality (-1, -1, 0, 1) x - 1 >= 0,
                 An inequality (1, 0, 0, 0) x + 0 >= 0)
        y3 * (-y3 + 1)^-1 * (-y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 3, 4]: (An equation -1 == 0,)
        0
        [1, 3, 4]: (An equation -1 == 0,)
        0
        [0, 1, 3, 4]: (An equation -1 == 0,)
        0
        [2, 3, 4]: (An equation (1, 1, 0, -1) x + 1 == 0,
                    An equation (1, 0, -1, 0) x + 0 == 0,
                    An inequality (0, 1, 0, 0) x + 0 >= 0,
                    An inequality (1, 0, 0, 0) x + 0 >= 0)
        y3 * (-y1*y3 + 1)^-1 * (-y0*y2*y3 + 1)^-1
        [0, 2, 3, 4]: (An equation -1 == 0,)
        0
        [1, 2, 3, 4]: (An equation -1 == 0,)
        0
        [0, 1, 2, 3, 4]: (An equation -1 == 0,)
        0

    TESTS::

        sage: generating_function_of_polyhedron(
        ....:     Polyhedron(ieqs=[(0, 0, 1, 0, 0), (-1, 1, -1, 0, 0)]),
        ....:     sort_factors=True)
        y0 * (-y0 + 1)^-1 * (-y2 + 1)^-1 * (-y3 + 1)^-1 * (-y0*y1 + 1)^-1
        sage: generating_function_of_polyhedron(
        ....:     Polyhedron(ieqs=[(0, 0, -1, 0, 1), (0, 0, 1, 0, 0),
        ....:                      (0, 1, 0, 0, -1), (-1, 1, -1, 0, 0)]),
        ....:     sort_factors=True)
        (-y0^2*y3 + y0*y3 + y0) *
        (-y0 + 1)^-1 * (-y2 + 1)^-1 * (-y0*y3 + 1)^-1 * (-y0*y1*y3 + 1)^-1

        sage: generating_function_of_polyhedron(
        ....:     Polyhedron(ieqs=[(0, 1, 0, -1, 0, 0), (0, 0, 0, 1, 0, 0)],
        ....:                eqns=[(0, 0, 0, 1, 0, -1), (0, 1, 0, 0, -1, 0),
        ....:                      (0, 1, -1, 0, 0, 0)]),
        ....:     sort_factors=True)
        1 * (-y0*y1*y3 + 1)^-1 * (-y0*y1*y2*y3*y4 + 1)^-1

    ::

        sage: G = generating_function_of_polyhedron(P2[0], sort_factors=True)
        sage: S = generating_function_of_polyhedron(P2[0], sort_factors=True,
        ....:                                       split=True)
        sage: sum(S) == G.value()
        True

        sage: G = generating_function_of_polyhedron(P2[1], sort_factors=True)
        sage: S = generating_function_of_polyhedron(P2[1], sort_factors=True,
        ....:                                       split=True)
        sage: sum(S) == G.value()
        True

        sage: G = generating_function_of_polyhedron(P3[0], sort_factors=True)
        sage: S = generating_function_of_polyhedron(P3[0], sort_factors=True,
        ....:                                       split=True)
        sage: sum(S) == G.value()
        True
    """
    import logging
    logger = logging.getLogger(__name__)

    from sage.combinat.permutation import Permutations
    from sage.geometry.polyhedron.constructor import Polyhedron

    if polyhedron.is_empty():
        from sage.structure.factorization import Factorization
        return Factorization([], unit=0)

    logger.info('generating_function_of_polyhedron: %s', polyhedron)

    if split is False:
        result = _generating_function_of_polyhedron_(polyhedron, **kwds)
        if result_as_tuple is True:
            return (result,)
        else:
            return result

    d = polyhedron.ambient_dim()
    if d <= 1:
        raise ValueError('Cannot do splitting with only '
                         'dimension {}.'.format(d))

    if split is True:
        split = iter(
            (Polyhedron(
                ieqs=[tuple(1 if i==b else (-1 if i==a or i==0 and a > b else 0)
                            for i in range(d+1))
                      for a, b in zip(pi[:-1], pi[1:])]),
             'b{}'.format(pi[0]-1) +
             ''.join((' <= ' if a < b else ' < ') +
                     'b{}'.format(b-1)
                     for a, b in zip(pi[:-1], pi[1:])))
            for pi in Permutations(d))
    else:
        split = iter((ph, ph.repr_pretty_Hrepresentation(prefix='b'))
                     for ph in split)

    result = []
    for split_polyhedron, pi_log in split:
        logger.info('split polyhedron by %s', pi_log)
        result.append(generating_function_of_polyhedron(
            polyhedron & split_polyhedron, split=False, **kwds))
    if result_as_tuple is False:
        raise ValueError('Cannot unpack result.')
    return result


def _generating_function_of_polyhedron_(
        polyhedron, indices=None, **kwds):
    r"""
    Helper function for :func:`generating_function_of_polyhedron` which
    does the actual computation of the generating function.
    """
    import logging
    logger = logging.getLogger(__name__)

    logger.info('polyhedron: %s',
                polyhedron.repr_pretty_Hrepresentation(prefix='b'))

    if polyhedron.is_empty():
        return Factorization([], unit=0)

    Hrepr = polyhedron.Hrepresentation()

    inequalities = tuple(tuple(entry)
                         for entry in Hrepr if entry.is_inequality())
    equations = tuple(tuple(entry)
                      for entry in Hrepr if entry.is_equation())
    if len(inequalities) + len(equations) != len(Hrepr):
        raise ValueError('Cannot handle {}.'.format(polyhedron))

    if not inequalities:
        raise NotImplementedError('no inequality given')

    if indices is None:
        indices = range(len(inequalities[0]) - 1)

    n = len(indices) + 1
    if any(len(ieq) != n for ieq in inequalities):
        raise ValueError('Not all coefficient vectors of the inequalities '
                         'have the same length.')
    if any(len(eq) != n for eq in equations):
        raise ValueError('Not all coefficient vectors of the equations '
                         'have the same length.')

    logger.info('generating_function_of_polyhedron: '
                '%s inequalities', len(inequalities))

    return __generating_function_of_polyhedron__(
        indices, inequalities, equations, **kwds)


def __generating_function_of_polyhedron__(
        indices, inequalities, equations,
        Factorization_sort=False, Factorization_simplify=False,
        sort_factors=False):
    r"""
    Helper function for :func:`generating_function_of_polyhedron` which
    does the actual computation of the generating function.
    """
    import logging
    logger = logging.getLogger(__name__)

    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    from sage.structure.factorization import Factorization

    B = LaurentPolynomialRing(
        ZZ,
        tuple('y{}'.format(k) for k in indices),
        sparse=True)

    gf_extra_factor_equations, rules_equations, indices_equations = \
            prepare_equations(equations, B)

    inequalities, gf_extra_factor_inequalities, rules_inequalities = \
        prepare_inequalities(inequalities, B)

    numerator = B(1)
    terms = B.gens()
    L = B
    for i, coeffs in enumerate(inequalities):
        L = LaurentPolynomialRing(L, 'mu{}'.format(i), sparse=True)
        l = L.gen()
        logger.debug('generating_function_of_polyhedron: '
                     '%s --> %s', l, pretty_inequality(coeffs))
        it_coeffs = iter(coeffs)
        numerator *= l**next(it_coeffs)
        assert numerator.parent() == L
        terms = tuple(l**c * t for c, t in zip(it_coeffs, terms))
    assert all(y == t for y, t in
               (zip(B.gens(), terms)[i] for i in indices_equations))
    terms = tuple(t for i, t in enumerate(terms)
                  if i not in indices_equations)

    logger.info('generating_function_of_polyhedron: '
                'terms denominator %s', terms)

    def decode_factor(factor):
        D = factor.dict()
        assert len(D) == 1
        exponent, coefficient = next(iteritems(D))
        return coefficient, exponent

    while repr(numerator.parent().gen()).startswith('mu'):
        logger.info('generating_function_of_polyhedron: '
                    'applying Omega[%s]...', numerator.parent().gen())
        logger.info('...on terms denominator %s', terms)
        logger.info('...(numerator has %s)', numerator.number_of_terms())

        decoded_factors, other_factors = \
            partition((decode_factor(factor) for factor in terms),
                      lambda factor: factor[1] == 0)
        other_factors = tuple(factor[0] for factor in other_factors)
        numerator, factors_denominator = \
            _Omega_(numerator.dict(), tuple(decoded_factors))
        terms = other_factors + factors_denominator

    numerator = (numerator.subs(rules_inequalities) *
                 gf_extra_factor_inequalities).subs(rules_equations) * \
                 gf_extra_factor_equations
    terms = tuple(t.subs(rules_inequalities).subs(rules_equations)
                  for t in terms)

    if sort_factors:
        def key(t):
            D = t.dict().popitem()[0]
            return (-sum(abs(d) for d in D), D)
        terms = sorted(terms, key=key, reverse=True)
    return Factorization([(numerator, 1)] +
                         list((1-t, -1) for t in terms),
                         sort=Factorization_sort,
                         simplify=Factorization_simplify)
