r"""
MacMahon's Omega Operator
"""
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

    - ``Factorization_sort`` (default: ``False``) and
      ``Factorization_simplify`` (default: ``True``) -- are passed on to
      :class:`sage.structure.factorization.Factorization` when creating
      the result.

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

        sage: from sage.rings.polynomial.omega import _Omega_
        sage: _Omega_({0: 2, 1: 40, -1: -3}, [])
        (42, ())
        sage: _Omega_({-1: 42}, [])
        (0, ())
    """
    if not decoded_factors:
        return sum(c for a, c in iteritems(A) if a >= 0), tuple()

    # Below we sort to make the caching more efficient. Doing this here
    # (in contrast to directly in Omega_ge) results in much cleaner
    # code and prevents an additional substitution or passing of a permutation.
    values, exponents = zip(*sorted(decoded_factors, key=lambda k: -k[1]))

    numerator = 0
    factors_denominator = None
    rules = None
    for a, c in iteritems(A):
        n, fd = Omega_ge(a, exponents)
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


@cached_function
def Omega_ge(a, exponents):
    r"""
    Return `\Omega_{\ge}` of the expression specified by the input.

    To be more precise, this calculates

    .. MATH::

        \Omega_{\ge} \frac{\mu^a}{
        (1 - z_1 \mu^{e_1}) \dots (1 - z_n \mu^{e_n})}

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

        sage: from sage.rings.polynomial.omega import Omega_ge
        sage: Omega_ge(0, (1, -2))
        (1, (z0, z0^2*z1))
        sage: Omega_ge(0, (1, -3))
        (1, (z0, z0^3*z1))
        sage: Omega_ge(0, (1, -4))
        (1, (z0, z0^4*z1))

        sage: Omega_ge(0, (2, -1))
        (z0*z1 + 1, (z0, z0*z1^2))
        sage: Omega_ge(0, (3, -1))
        (z0*z1^2 + z0*z1 + 1, (z0, z0*z1^3))
        sage: Omega_ge(0, (4, -1))
        (z0*z1^3 + z0*z1^2 + z0*z1 + 1, (z0, z0*z1^4))

        sage: Omega_ge(0, (1, 1, -2))
        (-z0^2*z1*z2 - z0*z1^2*z2 + z0*z1*z2 + 1, (z0, z1, z0^2*z2, z1^2*z2))
        sage: Omega_ge(0, (2, -1, -1))
        (z0*z1*z2 + z0*z1 + z0*z2 + 1, (z0, z0*z1^2, z0*z2^2))
        sage: Omega_ge(0, (2, 1, -1))
        (-z0*z1*z2^2 - z0*z1*z2 + z0*z2 + 1, (z0, z1, z0*z2^2, z1*z2))

    ::

        sage: Omega_ge(0, (2, -2))
        (-z0*z1 + 1, (z0, z0*z1, z0*z1))
        sage: Omega_ge(0, (2, -3))
        (z0^2*z1 + 1, (z0, z0^3*z1^2))
        sage: Omega_ge(0, (3, 1, -3))
        (-z0^3*z1^3*z2^3 + 2*z0^2*z1^3*z2^2 - z0*z1^3*z2
         + z0^2*z2^2 - 2*z0*z2 + 1,
         (z0, z1, z0*z2, z0*z2, z0*z2, z1^3*z2))

    ::

        sage: Omega_ge(0, (3, 6, -1))
        (-z0*z1*z2^8 - z0*z1*z2^7 - z0*z1*z2^6 - z0*z1*z2^5 - z0*z1*z2^4 +
         z1*z2^5 - z0*z1*z2^3 + z1*z2^4 - z0*z1*z2^2 + z1*z2^3 -
         z0*z1*z2 + z0*z2^2 + z1*z2^2 + z0*z2 + z1*z2 + 1,
         (z0, z1, z0*z2^3, z1*z2^6))

    TESTS::

        sage: Omega_ge(0, (2, 2, 1, 1, 1, 1, 1, -1, -1))[0].number_of_terms()  # long time
        27837

    ::

        sage: Omega_ge(1, (2,))
        (1, (z0,))
    """
    import logging
    logger = logging.getLogger(__name__)
    logger.info('Omega_ge: a=%s, exponents=%s', a, exponents)

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
    L = LaurentPolynomialRing(B, ('t',) + z_names, len(z_names) + 1)
    t = L.gens()[0]
    Z = LaurentPolynomialRing(ZZ, z_names, len(z_names))
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

    logger.debug('Omega_ge: preparing denominator')
    factors_denominator = tuple(de_power(1 - factor)
                                for factor in Omega_factors_denominator(x, y))

    logger.debug('Omega_ge: preparing numerator')
    numerator = de_power(Omega_numerator(a, x, y, t))

    logger.info('Omega_ge: completed')
    return numerator, factors_denominator


def Omega_numerator(a, x, y, t):
    r"""
    Return the numerator of `\Omega_{\ge}` of the expression
    specified by the input.

    To be more precise, this calculates

    .. MATH::

        \Omega_{\ge} \frac{\mu^a}{
        (1 - x_1 \mu) \dots (1 - x_n \mu)
        (1 - y_1 / \mu) \dots (1 - y_m / \mu)}

    and returns its numerator.

    This function is meant to be a helper function of :func:`Omega`.

    INPUT:

    - ``a`` -- an integer.

    - ``x`` and ``y`` -- a tuple of tuples of laurent polynomials. The
      flattened ``x`` contains `x_1,...,x_n`, the flattened ``y`` the
      `y_1,...,y_m`.

    - ``t`` -- a temporary laurent polynomial variable used for substituting.

    OUTPUT:

    A laurent polynomial.

    EXAMPLES::

        sage: from sage.rings.polynomial.omega import Omega_numerator, Omega_factors_denominator

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
                      sum(homogenous_symmetric_function(j, xy)
                          for j in srange(-a))
                      if a < 0 else 0)
    elif n == 0:
        result = sum(homogenous_symmetric_function(j, xy)
                     for j in srange(a+1))
    else:
        result = _Omega_numerator_P_(a, x_flat, y_flat, t).subs({t: x_flat[-1]})
    L = t.parent()
    result = L(result)

    logger.info('Omega_numerator: %s terms', result.number_of_terms())
    return result


def _Omega_numerator_P_(a, x, y, t):
    r"""
    Helper function for :func:`Omega_numerator`.

    INPUT:

    - ``a`` -- an integer.

    - ``x`` and ``y`` -- a tuple of laurent polynomials.

    - ``t`` -- a temporary laurent polynomial variable used for substituting.

    OUTPUT:

    A laurent polynomial.

    TESTS::

        sage: from sage.rings.polynomial.omega import _Omega_numerator_P_
        sage: L.<x0, x1, y0, y1, t> = LaurentPolynomialRing(ZZ)
        sage: _Omega_numerator_P_(0, (x0, x1), (y0,), t).subs({t: x1})
        -x0*x1*y0 + 1
    """
    import logging
    logger = logging.getLogger(__name__)

    from sage.arith.srange import srange
    from sage.misc.misc_c import prod

    n = len(x)
    if n == 1:
        x0 = t
        result = x0**(-a) + \
            (prod(1 - x0*yy for yy in y) *
             sum(homogenous_symmetric_function(j, y) * (1-x0**(j-a))
                 for j in srange(a))
             if a > 0 else 0)
    else:
        Pprev = _Omega_numerator_P_(a, x[:n-1], y, t)
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


@cached_function
def Omega_factors_denominator(x, y):
    r"""
    Return the denominator of `\Omega_{\ge}` of the expression
    specified by the input.

    To be more precise, this calculates

    .. MATH::

        \Omega_{\ge} \frac{1}{
        (1 - x_1 \mu) \dots (1 - x_n \mu)
        (1 - y_1 / \mu) \dots (1 - y_m / \mu)}

    and returns a factorization of its denominator.

    This function is meant to be a helper function of :func:`Omega`.

    INPUT:

    - ``x`` and ``y`` -- a tuple of tuples of laurent polynomials. The
      flattened ``x`` contains `x_1,...,x_n`, the flattened ``y`` the
      `y_1,...,y_m`.

    OUTPUT:

    A factorization of the denominator as
    a tuple of laurent polynomials.

    EXAMPLES::

        sage: from sage.rings.polynomial.omega import Omega_factors_denominator

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


def partition(items, predicate=bool):
    r"""
    Split ``items`` into two parts by the given ``predicate``.

    INPUT:

    - ``item`` -- an iterator.

    - ``predicate`` -- a function.

    OUTPUT:

    A pair of iterators; the first contains the elements not satisfying
    the ``predicate``, the second the elements satisfying the ``predicate``.

    Source of the code:
    `http://nedbatchelder.com/blog/201306/filter_a_list_into_two_parts.html
    <http://nedbatchelder.com/blog/201306/filter_a_list_into_two_parts.html>`_

    EXAMPLES::

        sage: from sage.rings.polynomial.omega import partition
        sage: E, O = partition(srange(10), is_odd)
        sage: tuple(E), tuple(O)
        ((0, 2, 4, 6, 8), (1, 3, 5, 7, 9))
    """
    from itertools import tee
    a, b = tee((predicate(item), item) for item in items)
    return ((item for pred, item in a if not pred),
            (item for pred, item in b if pred))


def homogenous_symmetric_function(j, x):
    r"""
    Return a complete homogeneous symmetric polynomial.

    INPUT:

    - ``j`` -- the degree as a nonnegative integer.

    - ``x`` -- an iterable of variables.

    OUTPUT:

    A polynomial whose type is determined by the input ``x``.

    EXAMPLES::

        sage: from sage.rings.polynomial.omega import homogenous_symmetric_function
        sage: P = PolynomialRing(ZZ, 'X', 3)
        sage: homogenous_symmetric_function(0, P.gens())
        1
        sage: homogenous_symmetric_function(1, P.gens())
        X0 + X1 + X2
        sage: homogenous_symmetric_function(2, P.gens())
        X0^2 + X0*X1 + X1^2 + X0*X2 + X1*X2 + X2^2
        sage: homogenous_symmetric_function(3, P.gens())
        X0^3 + X0^2*X1 + X0*X1^2 + X1^3 + X0^2*X2 +
        X0*X1*X2 + X1^2*X2 + X0*X2^2 + X1*X2^2 + X2^3
    """
    from sage.combinat.integer_vector import IntegerVectors
    from sage.misc.misc_c import prod

    return sum(prod(xx**pp for xx, pp in zip(x, p))
               for p in IntegerVectors(j, length=len(x)))
