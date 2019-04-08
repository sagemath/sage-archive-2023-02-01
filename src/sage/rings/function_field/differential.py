"""
Differentials of function fields

Sage provides basic arithmetic and advanced computations with differentials on
global function fields.

EXAMPLES:

The module of differentials on a function field forms an one-dimensional vector space over
the function field::

    sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
    sage: f = x + y
    sage: g = 1 / y
    sage: df = f.differential()
    sage: dg = g.differential()
    sage: dfdg = f.derivative() / g.derivative()
    sage: df == dfdg * dg
    True
    sage: df
    (x*y^2 + 1/x*y + 1) d(x)
    sage: df.parent()
    Space of differentials of Function field in y defined by y^3 + x^3*y + x

We can compute a canonical divisor::

    sage: k = df.divisor()
    sage: k.degree()
    4
    sage: k.degree() == 2 * L.genus() - 2
    True

Exact differentials vanish and logarithmic differentials are stable under the
Cartier operation::

    sage: df.cartier()
    0
    sage: w = 1/f * df
    sage: w.cartier() == w
    True

AUTHORS:

- Kwankyu Lee (2017-04-30): initial version

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

import operator

from sage.misc.latex import latex

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import ModuleElement
from sage.structure.richcmp import richcmp
from sage.sets.family import Family

from sage.categories.modules import Modules
from sage.categories.morphism import Morphism
from sage.categories.action import Action

from .function_field import FunctionField

class FunctionFieldDifferential(ModuleElement):
    """
    Base class for differentials on function fields.
    """
    pass

class FunctionFieldDifferential_global(FunctionFieldDifferential):
    """
    Differentials on global function fields.

    INPUT:

    - ``f`` -- element of the function field

    - ``t`` -- element of the function field; if `t` is not specified, `t`
      is the generator of the base rational function field

    EXAMPLES::

        sage: F.<x>=FunctionField(GF(7))
        sage: f = x/(x^2 + x + 1)
        sage: f.differential()
        ((6*x^2 + 1)/(x^4 + 2*x^3 + 3*x^2 + 2*x + 1)) d(x)

        sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
        sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
        sage: y.differential()
        (x*y^2 + 1/x*y) d(x)
    """
    def __init__(self, parent, f, t=None):
        """
        Initialize the differential `fdt`.

        TESTS::

            sage: F.<x>=FunctionField(GF(7))
            sage: f = x/(x^2 + x + 1)
            sage: w = f.differential()
            sage: TestSuite(w).run()
        """
        ModuleElement.__init__(self, parent)

        if t is not None:
            der = parent.function_field().derivation()
            f *= der(t)

        self._f = f

    def _repr_(self):
        """
        Return the string representation of the differential.

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
        F = self.parent().function_field()

        if self._f.is_zero(): # zero differential
            return '0'

        r =  'd({})'.format(F.base_field().gen())

        if self._f.is_one():
            return r

        return '({})'.format(self._f) + ' ' + r

    def _latex_(self):
        r"""
        Return a latex representation of the differential.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: w = y.differential()
            sage: latex(w)
            \left( x y^{2} + \frac{1}{x} y \right)\, dx
        """
        F = self.parent().function_field()

        if self._f.is_zero(): # zero differential
            return '0'

        r =  'd{}'.format(F.base_field().gen())

        if self._f.is_one():
            return r

        return '\\left(' + latex(self._f) + '\\right)\\,' + r

    def _richcmp_(self, other, op):
        """
        Compare the differential and the other differential with respect to the
        comparison operator.

        INPUT:

        - ``other`` -- differential

        - ``op`` -- comparison operator

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
        Return the sum of the differential and the other differential.

        INPUT:

        - ``other`` -- differential

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
        W = self.parent()
        return W.element_class(W, self._f + other._f)

    def _neg_(self):
        """
        Return the negation of the differential.

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
        W = self.parent()
        return W.element_class(W, -self._f)

    def _rmul_(self, f):
        """
        Return the differential multiplied by the element of the function
        field.

        INPUT:

        - ``f`` -- element of the function field

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
        W = self.parent()
        return W.element_class(W, f * self._f)

    def divisor(self):
        """
        Return the divisor of the differential.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: w = (1/y) * y.differential()
            sage: w.divisor()
            - Place (1/x, 1/x^3*y^2 + 1/x)
             - Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1)
             - Place (x, y)
             + Place (x + 2, y + 3)
             + Place (x^6 + 3*x^5 + 4*x^4 + 2*x^3 + x^2 + 3*x + 4, y + x^5)

            sage: F.<x> = FunctionField(QQ)
            sage: w = (1/x).differential()
            sage: w.divisor()
            -2*Place (x)
        """
        F = self.parent().function_field()
        x = F.base_field().gen()
        return self._f.divisor() + (-2) * F(x).divisor_of_poles() + F.different()

    def residue(self, place):
        """
        Return the residue of the differential at the place.

        INPUT:

        - ``place`` -- place of the function field

        OUTPUT:

        - an element of the residue field of the place

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
        r = g.valuation(place)
        if r >= 0:
            return R.zero()
        else:
            g_shifted = g * s**(-r)
            c = g_shifted.higher_derivative(-r-1, s)
            return to_R(c)

    def cartier(self):
        """
        Return the image of the differential by the Cartier operator.

        The Cartier operator operates on differentials. Let `x` be a separating
        element of the function field.  If a differential `\omega` is written
        `\omega=(f_0^p+f_1^px+\dots+f_{p-1}^px^{p-1})dx` (prime-power
        representation), then the Cartier operator maps `\omega` to
        `f_{p-1}dx`. It is known that this definition does not depend on the
        choice of `x`.

        The Cartier operator has interesting properties. Notably, the set of
        exact differentials is precisely the kernel of the Cartier operator and
        logarithmic differentials are stable under the Cartier operation.

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
        W = self.parent()
        F = W.function_field()
        der = F.higher_derivation()
        power_repr = der._prime_power_representation(self._f)
        return W.element_class(W, power_repr[-1])

    def monomial_coefficients(self, copy=True):
        """
        Return a dictionary whose keys are indices of basis
        elements in the support of ``self`` and whose values
        are the corresponding coefficients.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)
            sage: d = y.differential()
            sage: d
            ((4*x/(x^7 + 3))*y^2 + ((4*x^7 + 1)/(x^8 + 3*x))*y + x^4/(x^7 + 3)) d(x)
            sage: d.monomial_coefficients()
            {0: (4*x/(x^7 + 3))*y^2 + ((4*x^7 + 1)/(x^8 + 3*x))*y + x^4/(x^7 + 3)}
        """
        return {0: self._f}

class DifferentialsSpaceMorphism(Morphism):
    """
    Base class for morphisms between differential spaces of function fields.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: mor = L.space_of_differentials().coerce_map_from(K.space_of_differentials())
        sage: isinstance(mor, sage.rings.function_field.differential.DifferentialsSpaceMorphism)
        True

    TESTS::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 5*x^3)
        sage: mor2 = L.space_of_differentials().coerce_map_from(K.space_of_differentials())

        sage: mor == mor2
        False
        sage: mor == mor
        True
    """

    def _repr_(self):
        """
        Return the string representation of this morphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: L.space_of_differentials().coerce_map_from(K.space_of_differentials())
            Morphism of Differential Spaces:
              From: Space of differentials of Rational function field in x over Rational Field
              To:   Space of differentials of Function field in y defined by y^2 - x*y + 4*x^3

        """
        s = "Morphism of Differential Spaces:"
        s += "\n  From: {}".format(self.domain())
        s += "\n  To:   {}".format(self.codomain())
        return s

    def is_injective(self):
        """
        Return ``True``, since the morphism is injective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: L.space_of_differentials().coerce_map_from(K.space_of_differentials()).is_injective()
            True
        """
        return True

    def is_surjective(self):
        """
        Return ``False``, since the morphism is not surjective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: L.space_of_differentials().coerce_map_from(K.space_of_differentials()).is_surjective()
            False
        """
        return False

    def __hash__(self):
        r"""
        Return a hash value of this morphism.

        .. TODO::

            Make FunctionFieldDifferential hashable, so we can drop this
            method and use Morphism's default __hash__() method instead.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQbar); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: mor = L.space_of_differentials().coerce_map_from(K.space_of_differentials())
            sage: hash(mor)  # random

        """
        return hash((self.domain(), self.codomain()))

    def _call_(self, v):
        """
        Map ``v`` to the codomain differential space.

        INPUT:

        - ``v`` -- a differential in the domain differential space

        EXAMPLES::

            sage: K.<x> = FunctionField(QQbar); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: mor = L.space_of_differentials().coerce_map_from(K.space_of_differentials())
            sage: mor(x.differential()).parent()
            Space of differentials of Function field in y defined by y^2 - x*y + 4*x^3

        """
        codomain = self.codomain()
        return FunctionFieldDifferential_global(codomain, codomain.function_field()(v._f))

class DifferentialMultiplicationAction(Action):
    """
    Implement multiplication of a differential by a function.

    INPUT:

    - ``G`` -- must be a function field

    - ``S`` -- must be a space of differentials for a function field

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: mor = K.space_of_differentials().get_action(L)
        sage: isinstance(mor, sage.rings.function_field.differential.DifferentialMultiplicationAction)
        True

    """
    def __init__(self, G, S):
        """
        EXAMPLES::

            sage: from sage.rings.function_field.differential import DifferentialMultiplicationAction
            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: R.<x>=FunctionField(QQ);

            sage: DifferentialMultiplicationAction(K, L)
            Traceback (most recent call last):
            ...
            TypeError: (Function field in y defined by y^3 + x^3*y + x) must be a DifferentialsSpace


            sage: DifferentialMultiplicationAction(K.space_of_differentials(), L.space_of_differentials())
            Traceback (most recent call last):
            ...
            TypeError: (Space of differentials of Rational function field in x over Finite Field in z2 of size 2^2) must be a FunctionField

            sage: DifferentialMultiplicationAction(R, L.space_of_differentials())
            Traceback (most recent call last):
            ...
            TypeError: 'Rational function field in x over Rational Field'
               and 'Space of differentials of Function field in y defined by y^3 + x^3*y + x' are not compatible

            sage: DifferentialMultiplicationAction(K, L.space_of_differentials())
            Left action by Rational function field in x over Finite Field in z2 of size 2^2
               on Space of differentials of Function field in y defined by y^3 + x^3*y + x
        """
        if not isinstance(G, FunctionField):
            raise TypeError("(%s) must be a FunctionField" % repr(G))
        if not isinstance(S, DifferentialsSpace):
            raise TypeError("(%s) must be a DifferentialsSpace" % repr(S))
        if not G.space_of_differentials().has_coerce_map_from(S) \
           and not S.has_coerce_map_from(G.space_of_differentials()):
            raise TypeError("'%s' and '%s' are not compatible" % (repr(G), repr(S)))
        Action.__init__(self, G, S, True, operator.pow)

    def _act_(self, f, d):
        r"""
        Return the product ``f d``.

        INPUT:

        - ``f`` -- an element in a function field

        - ``d`` -- a differential in a compatible function field
        """
        if f.parent().space_of_differentials().has_coerce_map_from(d.parent()):
            return f * f.parent().space_of_differentials()(d)
        if d.parent().function_field().has_coerce_map_from(f.parent()):
            return d.parent().function_field()(f) * d
        raise TypeError("unsupported operands for DifferentialMultiplicationAction")

class DifferentialsSpace(UniqueRepresentation, Parent):
    """
    Space of differentials of a function field.

    INPUT:

    - ``field`` -- function field

    EXAMPLES::

        sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
        sage: L.<y>=K.extension(Y^3+x+x^3*Y)
        sage: L.space_of_differentials()
        Space of differentials of Function field in y defined by y^3 + x^3*y + x
    """
    Element = FunctionFieldDifferential_global

    def __init__(self, field, category=None):
        """
        Initialize the space of differentials of the function field.

        TESTS::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: W = L.space_of_differentials()
            sage: TestSuite(W).run()
        """
        Parent.__init__(self, base=field, category=Modules(field).FiniteDimensional().WithBasis().or_subcategory(category))

    def _repr_(self):
        """
        Return the string representation of the space of differentials.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: w = y.differential()
            sage: w.parent()
            Space of differentials of Function field in y defined by y^3 + x^3*y + x
        """
        return "Space of differentials of {}".format(self.base())

    def _element_constructor_(self, f):
        """
        Construct differential `df` in the space from `f`.

        INPUT:

        - ``f`` -- element of the function field

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: S = L.space_of_differentials()
            sage: S(y)
            (x*y^2 + 1/x*y) d(x)
            sage: S(y) in S
            True
            sage: S(1)
            0
        """
        if f in self.base():
            return self.element_class(self, self.base().one(), f)

        raise ValueError

    def _get_action_(self, G, op, self_on_left):
        """
        EXAMPLES::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: K.space_of_differentials().get_action(L)
            Left action by Function field in y defined by y^3 + x^3*y + x
               on Space of differentials of Rational function field in x over Finite Field in z2 of size 2^2
            sage: y * x.differential()
            (y) d(x)
        """
        if op is operator.mul and G is not self.function_field() and isinstance(G, FunctionField):
            if G.space_of_differentials().has_coerce_map_from(self):
                return DifferentialMultiplicationAction(G, self)

    def _coerce_map_from_(self, S):
        """
        Define coercions.

        We can coerce from any DifferentialsSpace whose underlying field
        can be coerced into our underlying field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: L.space_of_differentials().coerce_map_from(K.space_of_differentials())
            Morphism of Differential Spaces:
              From: Space of differentials of Rational function field in x over Rational Field
              To:   Space of differentials of Function field in y defined by y^2 - x*y + 4*x^3

        """
        if isinstance(S, DifferentialsSpace):
            if self.function_field().has_coerce_map_from(S.function_field()):
                return DifferentialsSpaceMorphism(S, self)

    def function_field(self):
        """
        Return the function field to which the space of differentials
        is attached.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: S = L.space_of_differentials()
            sage: S.function_field()
            Function field in y defined by y^3 + x^3*y + x
        """
        return self.base()

    def _an_element_(self):
        """
        Return a differential.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: S = L.space_of_differentials()
            sage: S.an_element()  # random
            (x*y^2 + 1/x*y) d(x)
        """
        F = self.base()
        return self.element_class(self, F.one(), F.an_element())

    def basis(self):
        """
        Return a basis.

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(4)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^3+x+x^3*Y)
            sage: S = L.space_of_differentials()
            sage: S.basis()
            Family (d(x),)
        """
        return Family([self.element_class(self, self.base().one())])

