
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


def Omega_P(a, x, y):
    if len(x) == 1:
        return x[0]**(-a) + \
            (prod(1 - x[0]*yy for yy in y) *
             sum(HomogenousSymmetricFunction(j, y) * (1-x[0]**(j-a))
                 for j in srange(a))
             if a > 0 else 0)

    return (x[-1] * (1-x[-2]) *
            prod(1 - x[-2]*yy for yy in y) *
            Omega_P(a, x[:-2] + x[-1:], y)
            -
            x[-2] * (1-x[-1]) *
            prod(1 - x[-1]*yy for yy in y) *
            Omega_P(a, x[:-1], y))  /  (x[-1] - x[-2])


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

    if not flat_y:
        numerator = 1 - (prod(factors_denominator) *
                         sum(HomogenousSymmetricFunction(j, flat_x)
                             for j in srange(-a))
                         if a < 0 else 0)
        factors_denominator = \
            tuple(tuple(1 - xx for xx in gx) for gx in x)

    elif not flat_x:
        numerator = sum(HomogenousSymmetricFunction(j, flat_x)
                        for j in srange(a+1))
        factors_denominator = (tuple(),)

    else:
        numerator = Omega_P(a, flat_x, flat_y)
        factors_denominator = \
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
        sage: Omega_higher(0, [(x, 2), (y, -1)])
        (x*y + 1, (-x + 1, -x*y^2 + 1))
        sage: Omega_higher(0, [(x, 1), (y, 1), (z, -2)])
        (-x^2*y*z - x*y^2*z + x*y*z + 1,
         (-x + 1, -y + 1, -x^2*z + 1, -y^2*z + 1))
        sage: Omega_higher(0, [(x, 1), (y, -3)])
        (1, (-x + 1, -x^3*y + 1))
        sage: Omega_higher(0, [(x, 1), (y, -4)])
        (1, (-x + 1, -x^4*y + 1))
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
    def subs_Omega(factor, v):
        f = Omega_map[v]
        value = f.value
        exponent = abs(f.exponent)
        p = tuple(v.dict().popitem()[0]).index(1)
        def subs_e(e):
            e = list(e)
            assert e[p] % exponent == 0
            e[p] = e[p] // exponent
            return tuple(e)
        P = factor.parent()
        result = P({subs_e(e): c for e, c in iteritems(factor.dict())})
        return result.subs({v: value})

    vars_Omega = L_high.gens()[:nv]
    def subs_all_Omega(factor):
        for v in vars_Omega:
            factor = subs_Omega(factor, v)
        return factor

    factors_denominator = tuple(subs_all_Omega(factor)
                                for factor in factors_denominator)

    from sage.rings.fraction_field import FractionField_generic
    if isinstance(numerator.parent(), FractionField_generic):
        numerator = subs_all_Omega(L_high(numerator.numerator())) / \
                    subs_all_Omega(L_high(numerator.denominator()))
    else:
        numerator = subs_all_Omega(numerator)

    return numerator, factors_denominator


def Omega(var, numerator, factors_denominator, operator=operator.ge):
    pass


class OmegaGroupElement(IndexedFreeAbelianGroup.Element):

    def __init__(self, parent, x, normalize=True):
        if normalize:
            from sage.misc.misc_c import prod
            L = parent.factor_ring()
            constant = prod(
                (f.constant_coefficient()**e for f, e in iteritems(x)), L.one())
            x = {f/f.constant_coefficient(): e for f, e in iteritems(x)}
        super(OmegaGroup.Element, self).__init__(parent, x)


    def Omega(var, operator=operator.ge):
        r"""
        """
        pass


class OmegaGroup(IndexedFreeAbelianGroup):

    Element = OmegaGroupElement

    @staticmethod
    def __classcall__(cls, base_ring, names):
        from sage.structure.unique_representation import UniqueRepresentation
        names = LaurentPolynomialRing(QQ, names).variable_names()
        return UniqueRepresentation.__classcall__(cls, base_ring, names)


    def __init__(self, base_ring, names):
        r"""
        EXAMPLES::

            sage: G = OmegaGroup(QQ, 'x, y'); G
            OmegaGroup in x, y over Rational Field
        """
        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
        L = LaurentPolynomialRing(base_ring, names)
        self._factor_ring_ = L
        super(OmegaGroup, self).__init__(indices=L, prefix='',
                                         names=names,
                                         bracket='(', latex_bracket='(')


    def factor_ring(self):
        return self._factor_ring_


    def _repr_(self):
        return 'OmegaGroup in {} over {}'.format(
            ', '.join(self.variable_names()),
            self.factor_ring().base_ring())


    def _element_constructor_(self, x=None):
        r"""
        TESTS::

            sage: G = OmegaGroup(QQ, 'mu, x, y')
            sage: var('mu, x, y')
            (mu, x, y)
            sage: G(1 / (1-mu*x) / (1-y/mu))
            (1 - mu^-1*y)^-1*(-mu*x + 1)^-1
       """
        if isinstance(x, Expression):
            import operator
            from sage.symbolic.operators import mul_vararg
            factors = []
            if x.operator() == mul_vararg:
                x = list(
                    (f.operands() if f.operator() == operator.pow else (f, 1))
                    for f in x.operands())
            elif x.operator() == operator.pow:
                x = [x.operands()]
            else:
                x = [(x, 1)]
            L = self.factor_ring()
            x = list((L(f), e) for f, e in x)
        return super(OmegaGroup, self)._element_constructor_(x)
