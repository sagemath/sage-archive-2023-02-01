"""
Differentials

This module provides differentials on function fields.

EXAMPLES::

    sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
    sage: x.differential()
    d(x)
    sage: w = y.differential(); w
    (x*y^2 + 1/x*y) d(x)
    sage: w.parent()
    Space of differentials of Function field in y defined by y^3 + x^3*y + x

AUTHORS:

- Kwankyu Lee (2016): initial version

"""
#*****************************************************************************
#       Copyright (C) 2016 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import ModuleElement
from sage.structure.sage_object import richcmp

from sage.categories.modules import Modules

from .function_field import is_RationalFunctionField

def differential(field, f, t=None):
    """
    Return the differential `fdt`.

    This is a helper function to construct differentials regardless of
    the type of function fields the differential belongs to.

    INPUT:

    - ``field`` -- function field to which the differential belongs

    - ``f``, ``t`` -- elements of the function field

    If ``t`` is ``None`` (default), then the differential is `fdx` where
    `x` is the generator of the base rational function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
        sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
        sage: x.differential()
        d(x)
        sage: y.differential()
        (x*y^2 + 1/x*y) d(x)
    """
    f = field(f)
    if t is not None:
        t = field(t)

    return FunctionFieldDifferential(field, f, t)

class FunctionFieldDifferential(ModuleElement):
    """
    Base class for differentials in function fields.
    """
    def __init__(self, field, f, t=None):
        """
        Initialize differential `fdt`.

        INPUT:

        - ``field`` -- function field to which this differential belongs

        - ``f``, ``t`` -- elements of the function field

        If ``t`` is ``None`` (default), then the differential is `fdx` where
        x is the generator of the base rational function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: y.differential()
            (x*y^2 + 1/x*y) d(x)

            sage: F.<x>=FunctionField(GF(7))
            sage: f = x/(x^2+x+1)
            sage: f.differential()
            ((6*x^2 + 1)/(x^4 + 2*x^3 + 3*x^2 + 2*x + 1)) d(x)
        """
        ModuleElement.__init__(self, field.space_of_differentials())

        if t is not None:
            der = field.derivation()
            f *= der(t)

        self._field = field
        self._f = f

    def _repr_(self):
        """
        Return the string representation of this differential.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: y.differential()
            (x*y^2 + 1/x*y) d(x)

            sage: F.<x>=FunctionField(QQ)
            sage: f = 1/x
            sage: f.differential()
            (-1/x^2) d(x)
        """
        r =  'd({})'.format(self._field.base_field().gen())
        if self._f != 1:
            r = '({})'.format(self._f) + ' ' + r
        return r

    def _richcmp_(self, other, op):
        """
        Compare this differential and ``other`` with respect to ``op``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: w1 = y.differential()
            sage: w2 = L(x).differential()
            sage: w3 = (x*y).differential()
            sage: w1 < w2
            False
            sage: w2 < w1
            True
            sage: w3 == x * w1 + y * w2
            True

            sage: F.<x>=FunctionField(QQ)
            sage: w1 = ((x^2+x+1)^10).differential()
            sage: w2 = (x^2+x+1).differential()
            sage: w1 < w2
            False
            sage: w1 > w2
            True
            sage: w1 == 10*(x^2+x+1)^9 * w2
            True
        """
        return richcmp(self._f, other._f, op)

    def _add_(self, other):
        """
        Return the sum of this differential and `other`.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: w1 = y.differential()
            sage: w2 = (1/y).differential()
            sage: w1 + w2
            (((x^3 + 1)/x^2)*y^2 + 1/x*y) d(x)

            sage: F.<x> = FunctionField(QQ)
            sage: w1 = (1/x).differential()
            sage: w2 = (x^2+1).differential()
            sage: w1 + w2
            ((2*x^3 - 1)/x^2) d(x)
        """
        return differential(self._field, self._f + other._f)

    def _neg_(self):
        """
        Return the negation of this differential.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: w1 = y.differential()
            sage: w2 = (-y).differential()
            sage: -w1 == w2
            True

            sage: F.<x> = FunctionField(QQ)
            sage: w1 = (1/x).differential()
            sage: w2 = (-1/x).differential()
            sage: -w1 == w2
            True
        """
        return differential(self._field, -self._f)

    def _rmul_(self, f):
        """
        Return this differential multiplied by function ``f``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: w1 = (1/y).differential()
            sage: w2 = (-1/y^2) * y.differential()
            sage: w1 == w2
            True

            sage: F.<x>=FunctionField(QQ)
            sage: w1 = (x^2*(x^2+x+1)).differential()
            sage: w2 = (x^2).differential()
            sage: w3 = (x^2+x+1).differential()
            sage: w1 == (x^2) * w3 + (x^2+x+1) * w2
            True
        """
        return differential(self._field, f * self._f)

    def divisor(self):
        """
        Return the divisor of the differential.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: w = (1/y) * y.differential()
            sage: w.divisor()
            -1*Place (1/x, 1/x^3*y^2 + 1/x)
             - Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1)
             - Place (x, y)
             + Place (x + 2, y + 3)
             + Place (x^6 + 3*x^5 + 4*x^4 + 2*x^3 + x^2 + 3*x + 4, y + x^5)

            sage: F.<x> = FunctionField(QQ)
            sage: w = (1/x).differential()
            sage: w.divisor()
            -2*Place (x)
        """
        F = self._field
        x = F.base_field().gen()
        return self._f.divisor() + (-2) * F(x).divisor_of_poles() + F.different()

    def residue(self, place):
        """
        Return the residue of the differential at the place.

        The residue is an element of the residue field of the place.

        EXAMPLES:

        We verify the residue theorem in a rational function field::

            sage: F.<x> = FunctionField(GF(4))
            sage: f = 0
            sage: while f == 0:
            ....:     f = F.random_element()
            sage: w = 1/f * f.differential()
            sage: d = f.divisor()
            sage: s = d.support()
            sage: sum([w.residue(p).trace() for p in s])
            0

        and also in an extension field::

            sage: K.<x> = FunctionField(GF(7)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: f = 0
            sage: while f == 0:
            ....:     f = L.random_element()
            sage: w = 1/f * f.differential()
            sage: d = f.divisor()
            sage: s = d.support()
            sage: sum([w.residue(p).trace() for p in s])
            0
        """
        R,fr_R,to_R = place._residue_field()

        # Step 1: compute f such that fds equals this differential.
        s = place.local_uniformizer()
        dxds = ~(s.derivative())
        g = self._f * dxds

        # Step 2: compute c that is the coefficient of s^-1 in
        # the power series expansion of f
        r = g._valuation(place)
        if r >= 0:
            return R(0)
        else:
            g_shifted = g * s**(-r)
            c = g_shifted.hasse_derivative(-r-1, s)
            return to_R(c)

    def cartier(self):
        """
        Return the image of the differential under Cartier operation.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: f = x/y
            sage: w = 1/f * f.differential()
            sage: w.cartier() == w
            True

            sage: F.<x> = FunctionField(GF(4))
            sage: f = x/(x^2+x+1)
            sage: w = 1/f * f.differential()
            sage: w.cartier() == w
            True
        """
        hasse = self._field.hasse_derivation()
        power_repr = hasse._prime_power_representation(self._f)
        return differential(self._field, power_repr[-1])

class DifferentialsSpace(Parent):
    """
    Space of differentials of a function field.

    EXAMPLES::

        sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
        sage: L.<y>=K.extension(Y^3+x+x^3*Y)
        sage: L.space_of_differentials()
        Space of differentials of Function field in y defined by y^3 + x^3*y + x
    """
    def __init__(self, field):
        """
        Initialize the space of differentials of the function ``field``.

        INPUT:

        - ``field`` -- function field to which this space belongs
        """
        Parent.__init__(self, base=field, category=Modules(field))
        self._field = field

    def _repr_(self):
        """
        Return the string representation of this space of differentials.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: w = y.differential()
            sage: w.parent()
            Space of differentials of Function field in y defined by y^3 + x^3*y + x
        """
        return "Space of differentials of %s"%(self._field,)

    def _element_constructor_(self, f):
        """
        Construct differential ``df`` in this space from ``f``.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: S = L.space_of_differentials()
            sage: S(y)
            (x*y^2 + 1/x*y) d(x)
            sage: S(y) in S
            True
        """
        return differential(self._field, 1, f)

    def function_field(self):
        """
        Return the function field to which this space of differentials
        is attached.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: S = L.space_of_differentials()
            sage: S.function_field()
            Function field in y defined by y^3 + x^3*y + x
        """
        return self._field


