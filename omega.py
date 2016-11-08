
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
from sage.groups.indexed_free_group import IndexedFreeAbelianGroup
from six import iteritems


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


def Omega_Fundamental(a, x, y):
    r"""
    EXAMPLES::

        sage: L.<x, y, z, w> = LaurentPolynomialRing(QQ)
        sage: Omega_Fundamental(0, [x], [y])
        (1, (-x + 1, -x*y + 1))
        sage: Omega_Fundamental(0, [x], [y, z])
        (1, (-x + 1, -x*y + 1, -x*z + 1))
        sage: Omega_Fundamental(0, [x, y], [z])
        (-x*y*z + 1, (-x + 1, -y + 1, -x*z + 1, -y*z + 1))
        sage: Omega_Fundamental(0, [x, y, z], [w])
        (x*y*z*w^2 + x*y*z*w - x*y*w - x*z*w - y*z*w + 1,
         (-x + 1, -y + 1, -z + 1, -x*w + 1, -y*w + 1, -z*w + 1))
        sage: Omega_Fundamental(0, [x, y], [z, w])
        (x^2*y*z*w + x*y^2*z*w - x*y*z*w - x*y*z - x*y*w + 1,
         (-x + 1, -y + 1, -x*z + 1, -x*w + 1, -y*z + 1, -y*w + 1))
    """
    if not y:
        factors_denominator = tuple(1 - xx for xx in x)
        return (1 - (prod(factors_denominator) *
                     sum(HomogenousSymmetricFunction(j, x)
                         for j in srange(-a))
                     if a < 0 else 0),
                factors_denominator)

    if not x:
        return (sum(HomogenousSymmetricFunction(j, x)
                    for j in srange(a+1)),
                tuple())

    return (Omega_P(a, x, y),
            tuple(1 - xx for xx in x) +
            tuple(1 - xx*yy for xx in x for yy in y))


class OmegaGroupElement(IndexedFreeAbelianGroup.Element):

    def __init__(self, parent, x, normalize=True):
        if normalize:
            from sage.misc.misc_c import prod
            L = parent.factor_ring()
            constant = prod(
                (f.constant_coefficient()**e for f, e in iteritems(x)), L.one())
            x = {f/f.constant_coefficient(): e for f, e in iteritems(x)}
        super(OmegaGroup.Element, self).__init__(parent, x)


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
