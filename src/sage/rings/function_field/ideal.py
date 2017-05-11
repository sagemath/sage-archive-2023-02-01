r"""
Ideals

Ideals of an order of a function field include all fractional ideals of the order.
Sage provides basic arithmetic with fractional ideals.

The fractional ideals of the maximal order of a global function field forms a multiplicative
monoid. Sage allows advanced arithmetic with the fractional ideals. For example, an ideal
of the maximal order can be factored into a product of prime ideals.

EXAMPLES:

Ideals in the maximal order of a rational function field::

    sage: K.<x> = FunctionField(QQ)
    sage: O = K.maximal_order()
    sage: I = O.ideal(x^3 + 1); I
    Ideal (x^3 + 1) of Maximal order of Rational function field in x over Rational Field
    sage: I^2
    Ideal (x^6 + 2*x^3 + 1) of Maximal order of Rational function field in x over Rational Field
    sage: ~I
    Ideal (1/(x^3 + 1)) of Maximal order of Rational function field in x over Rational Field
    sage: ~I * I
    Ideal (1) of Maximal order of Rational function field in x over Rational Field

Ideals in the equation order of an extension of a rational function field::

    sage: K.<x> = FunctionField(QQ); R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - x^3 - 1)
    sage: O = L.equation_order()
    sage: I = O.ideal(y); I
    Ideal (x^3 + 1, -y) of Order in Function field in y defined by y^2 - x^3 - 1
    sage: I^2
    Ideal (x^3 + 1, (-x^3 - 1)*y) of Order in Function field in y defined by y^2 - x^3 - 1

Ideals in the maximal order of a global function field::

    sage: K.<x> = FunctionField(GF(2)); R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - x^3*y - x)
    sage: O = L.maximal_order()
    sage: I = O.ideal(y)
    sage: I^2
    Ideal (x) of Maximal order of Function field in y defined by y^2 + x^3*y + x
    sage: ~I
    Ideal (1, 1/x*y) of Maximal order of Function field in y defined by y^2 + x^3*y + x
    sage: ~I * I
    Ideal (1) of Maximal order of Function field in y defined by y^2 + x^3*y + x

    sage: J = O.ideal(x+y) * I
    sage: J.factor()
    (Ideal (x, y) of Maximal order of Function field in y defined by y^2 + x^3*y + x)^2 *
    (Ideal (x^3 + x + 1, y + x) of Maximal order of Function field in y defined by y^2 + x^3*y + x)

Ideals in the maximal infinite order of a global function field::

    sage: K.<x> = FunctionField(GF(3^2)); R.<t> = K[]
    sage: F.<y> = K.extension(t^3 + t^2 - x^4)
    sage: Oinf = F.maximal_order_infinite()
    sage: I = Oinf.ideal(1/y)
    sage: I + I == I
    True
    sage: I^2
    Ideal (1/x^3,1/x^4*y) of Maximal infinite order of Function field in y defined by y^3 + y^2 + 2*x^4
    sage: ~I
    Ideal (x,y) of Maximal infinite order of Function field in y defined by y^3 + y^2 + 2*x^4
    sage: ~I * I
    Ideal (1) of Maximal infinite order of Function field in y defined by y^3 + y^2 + 2*x^4
    sage: I.factor()
    (Ideal (1/x,1/x^3*y^2) of Maximal infinite order of Function field in y defined by y^3 + y^2 + 2*x^4)^4

AUTHORS:

- William Stein (2010): initial version

- Maarten Derickx (2011-09-14): fixed ideal_with_gens_over_base()

- Kwankyu Lee (2017-04-30): added ideals for global function fields

"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#       Copyright (C) 2011 Maarten Derickx <m.derickx.student@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator
import itertools

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod

from sage.structure.sage_object import SageObject, richcmp
from sage.structure.parent import Parent
from sage.structure.factorization import Factorization
from sage.structure.element import AlgebraElement
from sage.modules.free_module_element import vector

from sage.categories.monoids import Monoids

import sage.rings.integer_ring
from sage.rings.infinity import infinity

from .place import FunctionFieldPlace_rational
from .place import FunctionFieldPlace_global

lazy_import('sage.matrix.constructor', 'matrix')
lazy_import('sage.rings.function_field.divisor', 'FunctionFieldDivisor')

def is_Ideal(x):
    """
    Return True if ``x`` is an ideal in a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(7))
        sage: O = K.equation_order()
        sage: I = O.ideal(x^3+1)
        sage: from sage.rings.function_field.ideal import is_Ideal
        sage: is_Ideal(I)
        True
    """
    return isinstance(x, FunctionFieldIdeal)

class FunctionFieldIdeal(AlgebraElement):
    """
    Fractional ideals of function fields.
    """
    def __init__(self, ring):
        """
        Initialize.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: O = K.equation_order()
            sage: I = O.ideal(x^3+1)
            sage: isinstance(I, sage.rings.function_field.ideal.FunctionFieldIdeal)
            True
        """
        AlgebraElement.__init__(self, ring.ideal_monoid())
        self._ring = ring

    def _div_(self, other):
        """
        Return the ideal divided by the other ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: O = K.equation_order()
            sage: I = O.ideal(x^3+1)
            sage: I / I
            Ideal (1) of Maximal order of Rational function field in x
            over Finite Field of size 7
        """
        return self * ~other

    def gens_reduced(self):
        r"""
        Return reduced generators. This just returns the generators for now.

        This method is provided so that ideals in funtion fields have the method
        :meth:`gens_reduced()`, just like ideals of number fields. Sage linear algebra
        machinery sometimes requires this.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: O = K.equation_order()
            sage: I = O.ideal(x,x^2,x^2+x)
            sage: I.gens_reduced()
            (x,)
        """
        return self.gens()

    def ring(self):
        """
        Return the ring to which the ideal belongs.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: O = K.equation_order()
            sage: I = O.ideal(x,x^2,x^2+x)
            sage: I.ring()
            Maximal order of Rational function field in x over Finite Field of size 7
        """
        return self._ring

class FunctionFieldIdeal_rational(FunctionFieldIdeal):
    """
    Fractional ideals of the maximal order of a rational function field.
    """
    def __init__(self, ring, gen):
        """
        Initialize.

        INPUT:

        - ``ring`` -- the maximal order of the rational function field.

        - ``gen`` -- generator of the ideal, an element of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(1/(x^2+x)); I
            Ideal (1/(x^2 + x)) of Maximal order of Rational function field in x over Rational Field
        """
        FunctionFieldIdeal.__init__(self, ring)
        self._gen = gen

    def __hash__(self):
        """
        Return hash computed from the data.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(1/(x^2+x))
            sage: d = { I: 1, I^2: 2 }
        """
        return hash( (self._ring, self._gen) )

    def _repr_(self):
        """
        Return string representation.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x,1/(x+1)); I
            Ideal (1/(x + 1)) of Maximal order of Rational function field in x over Rational Field
        """
        return "Ideal (%s) of %s"%(self._gen, self._ring)

    def __contains__(self, element):
        """
        Test if ``element`` is in the ideal.

        INPUT:

        - ``element`` -- element of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(1/(x+1))
            sage: x in I
            True
        """
        return (element / self._gen) in self._ring

    def _richcmp_(self, other, op):
        """
        Compare the element with the other element with respect
        to the comparison operator.

        INPUT:

        - ``other`` -- element

        - ``op`` -- comparison operator

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x,x^2+1)
            sage: J = O.ideal(x^2+x+1,x)
            sage: I == J
            True
            sage: I = O.ideal(x)
            sage: J = O.ideal(x+1)
            sage: I < J
            True
        """
        return richcmp(self._gen, other._gen, op)

    def _add_(self, other):
        """
        Add the ideal with the other ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x,x^2+1)
            sage: J = O.ideal(x^2+x+1,x)
            sage: I + J == J + I
            True
        """
        return self._ring.ideal([self._gen, other._gen])

    def _mul_(self, other):
        """
        Multiply the ideal with the other ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x,x^2+x)
            sage: J = O.ideal(x^2,x)
            sage: I * J == J * I
            True
        """
        return self._ring.ideal([self._gen * other._gen])

    def _rmul_(self, other):
        """
        Multiply ``other`` with the ideal on the right.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^3+x^2)
            sage: 2 * I
            Ideal (x^3 + x^2) of Maximal order of Rational function field in x over Rational Field
            sage: x * I
            Ideal (x^4 + x^3) of Maximal order of Rational function field in x over Rational Field
        """
        return self._ring.ideal([other * self._gen])

    def __invert__(self):
        """
        Return the ideal inverse of the fractional ideal.

        EXAMPLES::

            sage: F.<x> = FunctionField(QQ)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x/(x^2+1))
            sage: ~I
            Ideal ((x^2 + 1)/x) of Maximal order of Rational function field
            in x over Rational Field
        """
        return self._ring.ideal([~(self._gen)])

    def denominator(self):
        """
        Return the denominator of the fractional ideal.

        EXAMPLES::

            sage: F.<x> = FunctionField(QQ)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x/(x^2+1))
            sage: I.denominator()
            x^2 + 1
        """
        return self._gen.denominator()

    def is_prime(self):
        """
        Return True if this is a prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^3+x^2)
            sage: [f.is_prime() for f,m in I.factor()]
            [True, True]
        """
        return self._gen.denominator() == 1 and self._gen.numerator().is_prime()

    @cached_method
    def module(self):
        """
        Return the module, that is the ideal viewed as a module over the ring.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^3+x^2)
            sage: I.module()
            Free module of degree 1 and rank 1 over Maximal order of Rational
            function field in x over Rational Field
            Echelon basis matrix:
            [x^3 + x^2]
            sage: J = 0*I
            sage: J.module()
            Free module of degree 1 and rank 0 over Maximal order of Rational
            function field in x over Rational Field
            Echelon basis matrix:
            []
        """
        V, fr, to = self.ring().fraction_field().vector_space()
        return V.span([to(g) for g in self.gens()], base_ring=self.ring())

    def gen(self):
        """
        Return the unique generator of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2+x)
            sage: I.gen()
            x^2 + x
        """
        return self._gen

    def gens(self):
        """
        Return the tuple of the unique generator of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2+x)
            sage: I.gens()
            (x^2 + x,)
        """
        return (self._gen,)

    def gens_over_base(self):
        """
        Return the generator of the ideal as a rank one module over the maximal
        order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2+x)
            sage: I.gens_over_base()
            (x^2 + x,)
        """
        return (self._gen,)

    def valuation(self, ideal):
        """
        Return valuation value of the ideal at the prime ideal.

        INPUT:

        - ``ideal`` -- fractional ideal

        EXAMPLES::

            sage: F.<x> = FunctionField(QQ)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x^2*(x^2+x+1)^3)
            sage: [f.valuation(I) for f,_ in I.factor()]
            [2, 3]
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        O = self.ring()
        d = ideal.denominator()
        return self._valuation(d*ideal) - self._valuation(O.ideal(d))

    def _valuation(self, ideal):
        """
        Return valuation of the integral ideal at the prime ideal.

        INPUT:

        - ``ideal`` -- ideal

        EXAMPLES::

            sage: F.<x> = FunctionField(QQ)
            sage: O = F.maximal_order()
            sage: p = O.ideal(x)
            sage: p.valuation(O.ideal(x+1))  # indirect doctest
            0
            sage: p.valuation(O.ideal(x^2))  # indirect doctest
            2
            sage: p.valuation(O.ideal(1/x^3))  # indirect doctest
            -3
            sage: p.valuation(O.ideal(0))  # indirect doctest
            +Infinity
        """
        return ideal.gen().valuation(self.gen())

    def factor(self):
        """
        Return the factorization of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^3*(x+1)^2)
            sage: I.factor()
            (Ideal (x) of Maximal order of Rational function field in x
            over Finite Field in z2 of size 2^2)^3 *
            (Ideal (x + 1) of Maximal order of Rational function field in x
            over Finite Field in z2 of size 2^2)^2
        """
        factors = self._factor()
        return Factorization(factors, cr=True)

    def _factor(self):
        """
        Return the list of prime and multiplicity pairs of the
        factorization of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^3*(x+1)^2)
            sage: I._factor()
            [(Ideal (x) of Maximal order of Rational function field in x
            over Finite Field in z2 of size 2^2, 3),
            (Ideal (x + 1) of Maximal order of Rational function field in x
            over Finite Field in z2 of size 2^2, 2)]
        """
        factors = []
        for f,m in self._gen.factor():
             factors.append( (self.ring().ideal(f), m) )
        return factors

    def place(self):
        """
        Return the place associated with the prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2+x+1)
            sage: I.place()
            Traceback (most recent call last):
            ...
            TypeError: not a prime ideal
            sage: I = O.ideal(x^3+x+1)
            sage: I.place()
            Place (x^3 + x + 1)
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")
        return FunctionFieldPlace_rational(self.ring().fraction_field(), self)

    def divisor(self):
        """
        Return divisor corresponding to the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x*(x+1)^2/(x^2+x+1))
            sage: I.divisor()
            Place (x) + 2*Place (x + 1) - Place (x + z2) - Place (x + z2 + 1)
        """
        if self._gen == 0:
            raise ValueError("not defined for zero element")

        data = {prime.place(): multiplicity for prime, multiplicity in self._factor()}
        return FunctionFieldDivisor(self.ring().fraction_field(), data)

    def divisor_of_zeros(self):
        """
        Return divisor of zeros corresponding to the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x*(x+1)^2/(x^2+x+1))
            sage: I.divisor_of_zeros()
            Place (x) + 2*Place (x + 1)
        """
        if self._gen == 0:
            raise ValueError("not defined for zero element")

        data = {prime.place(): multiplicity for prime, multiplicity in self._factor() if multiplicity > 0}
        return FunctionFieldDivisor(self.ring().fraction_field(), data)

    def divisor_of_poles(self):
        """
        Return divisor corresponding to the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x*(x+1)^2/(x^2+x+1))
            sage: I.divisor_of_poles()
            Place (x + z2) + Place (x + z2 + 1)
        """
        if self._gen == 0:
            raise ValueError("not defined for zero element")

        data = {prime.place(): - multiplicity for prime, multiplicity in self._factor() if multiplicity < 0}
        return FunctionFieldDivisor(self.ring().fraction_field(), data)

class FunctionFieldIdeal_module(FunctionFieldIdeal):
    """
    A fractional ideal specified by a finitely generated module over
    the integers of the base field.

    EXAMPLES:

    An ideal in an extension of a rational function field::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x^3 - 1)
        sage: O = L.equation_order()
        sage: I = O.ideal(y)
        sage: I
        Ideal (x^3 + 1, -y) of Order in Function field in y defined by y^2 - x^3 - 1
        sage: I^2
        Ideal (x^3 + 1, (-x^3 - 1)*y) of Order in Function field in y defined by y^2 - x^3 - 1
    """
    def __init__(self, ring, module):
        """
        Initialize.

        INPUT:

        - ``ring`` -- order in a function field

        - ``module`` -- module

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(y)
            sage: type(I)
            <class 'sage.rings.function_field.ideal.FunctionFieldIdeal_module'>
        """
        FunctionFieldIdeal.__init__(self, ring)

        self._module = module
        self._structure = ring.fraction_field().vector_space()

        V, from_V, to_V = self._structure
        self._gens = tuple([from_V(a) for a in module.basis()])

    def __contains__(self, x):
        """
        Return True if ``x`` is in the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal_with_gens_over_base([1, y]);  I
            Ideal (1, y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: y in I
            True
            sage: y/x in I
            False
            sage: y^2 - 2 in I
            True
        """
        return self._structure[2](x) in self._module

    def __hash__(self):
        """
        Return the hash of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal_with_gens_over_base([1, y])
            sage: d = {I: 1}  # indirect doctest
        """
        return hash((self._ring,self._module))

    def __eq__(self, other):
        """
        Test equality of the ideal with the other ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal_with_gens_over_base([1, y])
            sage: I + I == I  # indirect doctest
            True
        """
        if not isinstance(other, FunctionFieldIdeal_module):
            other = self.ring().ideal(other)
        if self.ring() != other.ring():
            raise ValueError("rings must be the same")

        if (self.module().is_submodule(other.module()) and
            other.module().is_submodule(self.module())):
            return True
        else:
            return False

    def __cmp__(self, other):
        """
        Compare the ideal and the other ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y*(y+1)); J = O.ideal((y^2-2)*(y+1))
            sage: I+J == J+I            # indirect test
            True
            sage: I == J
            False
            sage: I < J
            True
            sage: J < I
            False
        """
        if not isinstance(other, FunctionFieldIdeal_module):
            other = self.ring().ideal(other)
        if self.ring() != other.ring():
            raise ValueError("rings must be the same")

        return cmp(self.module(), other.module())

    def module(self):
        """
        Return module over the maximal order of the base field that
        underlies the ideal.

        The formation of the module is compatible with the vector
        space corresponding to the function field.

        OUTPUT:

        - a module over the maximal order of the base field of the ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: O = K.maximal_order(); O
            Maximal order of Rational function field in x over Finite Field of size 7
            sage: K.polynomial_ring()
            Univariate Polynomial Ring in x over Rational function field in x over Finite Field of size 7
            sage: I = O.ideal(x^2 + 1)
            sage: I.gens()
            (x^2 + 1,)
            sage: I.module()
            Free module of degree 1 and rank 1 over Maximal order of Rational function field in x over Finite Field of size 7
            Echelon basis matrix:
            [x^2 + 1]
            sage: V, from_V, to_V = K.vector_space(); V
            Vector space of dimension 1 over Rational function field in x over Finite Field of size 7
            sage: I.module().is_submodule(V)
            True
        """
        return self._module

    def __repr__(self):
        """
        Return a string representation of the ideal.

        EXAMPLES::

            sage: P.<a,b,c> = QQ[]
            sage: P*[a^2,a*b+c,c^3] # indirect doctest
            Ideal (a^2, a*b + c, c^3) of Multivariate Polynomial Ring in a, b, c over Rational Field
        """
        return "Ideal (%s) of %s"%(', '.join([repr(g) for g in self.gens()]), self.ring())

    def __nonzero__(self):
        r"""
        Return True if the ideal is not a zero ideal.

        TESTS::

            sage: I = ZZ.ideal(5)
            sage: bool(I)
            True

        ::

            sage: I = ZZ['x'].ideal(0)
            sage: bool(I)
            False

        ::

            sage: I = ZZ['x'].ideal(ZZ['x'].gen()^2)
            sage: bool(I)
            True

        ::

            sage: I = QQ['x', 'y'].ideal(0)
            sage: bool(I)
            False
        """
        for g in self.gens():
            if not g.is_zero():
                return True
        return False

    def base_ring(self):
        r"""
        Returns the base ring of the ideal.

        EXAMPLES::

            sage: R = ZZ
            sage: I = 3*R; I
            Principal ideal (3) of Integer Ring
            sage: J = 2*I; J
            Principal ideal (6) of Integer Ring
            sage: I.base_ring(); J.base_ring()
            Integer Ring
            Integer Ring

        We construct an example of an ideal of a quotient ring::

            sage: R = PolynomialRing(QQ, 'x'); x = R.gen()
            sage: I = R.ideal(x^2 - 2)
            sage: I.base_ring()
            Rational Field

        And `p`-adic numbers::

            sage: R = Zp(7, prec=10); R
            7-adic Ring with capped relative precision 10
            sage: I = 7*R; I
            Principal ideal (7 + O(7^11)) of 7-adic Ring with capped relative precision 10
            sage: I.base_ring()
            7-adic Ring with capped relative precision 10
        """
        return self.ring().base_ring()

    def apply_morphism(self, phi):
        r"""
        Apply the morphism ``phi`` to every element of the ideal and
        return an ideal in the domain of ``phi``.

        INPUT:

        - ``phi`` -- morphism

        EXAMPLES::

            sage: psi = CC['x'].hom([-CC['x'].0])
            sage: J = ideal([CC['x'].0 + 1]); J
            Principal ideal (x + 1.00000000000000) of Univariate Polynomial Ring in x over Complex Field with 53 bits of precision
            sage: psi(J)
            Principal ideal (x - 1.00000000000000) of Univariate Polynomial Ring in x over Complex Field with 53 bits of precision
            sage: J.apply_morphism(psi)
            Principal ideal (x - 1.00000000000000) of Univariate Polynomial Ring in x over Complex Field with 53 bits of precision

        ::

            sage: psi = ZZ['x'].hom([-ZZ['x'].0])
            sage: J = ideal([ZZ['x'].0, 2]); J
            Ideal (x, 2) of Univariate Polynomial Ring in x over Integer Ring
            sage: psi(J)
            Ideal (-x, 2) of Univariate Polynomial Ring in x over Integer Ring
            sage: J.apply_morphism(psi)
            Ideal (-x, 2) of Univariate Polynomial Ring in x over Integer Ring

        TESTS::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: A = K.ideal(a)
            sage: taus = K.embeddings(K)
            sage: A.apply_morphism(taus[0]) # identity
            Fractional ideal (a)
            sage: A.apply_morphism(taus[1]) # complex conjugation
            Fractional ideal (-a)
            sage: A.apply_morphism(taus[0]) == A.apply_morphism(taus[1])
            True

        ::

            sage: K.<a> = NumberField(x^2 + 5)
            sage: B = K.ideal([2, a + 1]); B
            Fractional ideal (2, a + 1)
            sage: taus = K.embeddings(K)
            sage: B.apply_morphism(taus[0]) # identity
            Fractional ideal (2, a + 1)

        Since 2 is totally ramified, complex conjugation fixes it::

            sage: B.apply_morphism(taus[1]) # complex conjugation
            Fractional ideal (2, a + 1)
            sage: taus[1](B)
            Fractional ideal (2, a + 1)
        """
        from sage.categories.morphism import is_Morphism
        if not is_Morphism(phi):
            raise TypeError("phi must be a morphism")
        # delegate: morphisms know how to apply themselves to ideals
        return phi(self)

    def gens(self):
        """
        Return a set of generators of the ideal.

        This is the set of generators provided during the creation of the ideal.

        EXAMPLE::

            sage: P.<x,y> = PolynomialRing(QQ,2)
            sage: I = Ideal([x,y+1]); I
            Ideal (x, y + 1) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: I.gens()
            [x, y + 1]

        ::

            sage: ZZ.ideal(5,10).gens()
            (5,)
        """
        return self._gens

    def gen(self, i):
        """
        Return the ``i``-th generator in the current basis of the ideal.

        EXAMPLE::

            sage: P.<x,y> = PolynomialRing(QQ,2)
            sage: I = Ideal([x,y+1]); I
            Ideal (x, y + 1) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: I.gen(1)
            y + 1

            sage: ZZ.ideal(5,10).gen()
            5
        """
        return self._gens[i]

    def ngens(self):
        """
        Return the number of generators in the basis.

        EXAMPLE::

            sage: P.<x,y> = PolynomialRing(QQ,2)
            sage: I = Ideal([x,y+1]); I
            Ideal (x, y + 1) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: I.ngens()
            2

            sage: ZZ.ideal(5,10).ngens()
            1
        """
        return len(self.__gens)

    def __add__(self, other):
        """
        This method just makes sure that the ideal and ``other`` are ideals in the
        same ring and then calls :meth:`_add_`. If you want to change the
        behaviour of ideal addition in a subclass of
        :class:`Ideal_generic` please overwrite :meth:`_add_` and not
        :meth:`__add__`.

        EXAMPLE::

            sage: P.<x,y> = PolynomialRing(QQ,2)
            sage: I = Ideal([x,y+1])
            sage: I + I == I  # indirect doctest
            True
        """
        if not is_Ideal(other):
            other = self.ring().ideal(other)
        return self._add_(other)

    def _add_(self, other):
        """
        Add the ideal on the left to the other ideal.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: I = [x + y]*P
            sage: I + [y + z]
            Ideal (x + y, y + z) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        return self.ring().ideal(self.gens() + other.gens())

    def __radd__(self, other):
        """
        Add the ideal on the right to the other ideal.

        This makes sure that ``other`` is in the same ring of the ideal.

        EXAMPLES::

            sage: P.<x,y,z> = QQ[]
            sage: I = [x + y]*P
            sage: [y + z] + I
            Ideal (x + y, y + z) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        if not is_Ideal(other):
            other = self.ring().ideal(other)
        return self.ring().ideal(self.gens() + other.gens())

    def __mul__(self, other):
        """
        This method just makes sure that ``other`` is an ideal in the same ring
        with the ideal and then calls :meth:`_mul_`. If you want to change the
        behaviour of ideal multiplication in a subclass of
        :class:`Ideal_generic` please overwrite :meth:`_mul_` and not
        :meth:`__mul__`.

        EXAMPLE::

            sage: P.<x,y,z> = QQ[]
            sage: I = [x*y + y*z, x^2 + x*y - y*x - y^2] * P
            sage: I * 2    # indirect doctest
            Ideal (2*x*y + 2*y*z, 2*x^2 - 2*y^2) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        if not is_Ideal(other):
            try:
                if self.ring().has_coerce_map_from(other):
                    return self
            except (TypeError,ArithmeticError,ValueError):
                pass
            other = self.ring().ideal(other)
        return self._mul_(other)

    def _mul_(self, other):
        """
        This is a very general implementation of Ideal multiplication.

        The number of generators of ```self * other` will be
        ``self.ngens() * other.ngens()``. So if used repeatedly this method
        will create an ideal with a uselessly large amount of generators.
        Therefore it is advisable to overwrite this method with a method that
        takes advantage of the structure of the ring your working in.

        Example::

            sage: P.<x,y,z> = QQ[]
            sage: I=P.ideal([x*y, x*z, x^2])
            sage: J=P.ideal([x^2, x*y])
            sage: I._mul_(J)
            Ideal (x^3*y, x^2*y^2, x^3*z, x^2*y*z, x^4, x^3*y) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        return self.ring().ideal([x*y for x in self.gens() for y in other.gens()])

    def __rmul__(self, other):
        """
        Multiply the ideal on the right with ``other``.

        EXAMPLE::

            sage: P.<x,y,z> = QQ[]
            sage: I = [x*y+y*z,x^2+x*y-y*x-y^2]*P
            sage: [2]*I    # indirect doctest
            Ideal (2*x*y + 2*y*z, 2*x^2 - 2*y^2) of Multivariate Polynomial Ring in x, y, z over Rational Field

        """
        if not is_Ideal(other):
            try:
                if self.ring().has_coerce_map_from(other):
                    return self
            except (TypeError,ArithmeticError,ValueError):
                pass
            other = self.ring().ideal(other)
        return self.ring().ideal([z for z in [y*x for x in self.gens() for y in other.gens()] if z])

class FunctionFieldIdeal_global(FunctionFieldIdeal):
    """
    Fractional ideals of canonical function fields
    """
    def __init__(self, ring, hnf, denominator=1):
        """
        Initialize.

        INPUT:

        - ``ring`` -- order in a function field

        - ``hnf`` -- matrix in hermite normal form

        - ``denominator`` -- denominator

        The rows of ``hnf`` is a basis of the ideal, which itself is
        ``denominator`` times the fractional ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3*y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: type(I)
            <class 'sage.rings.function_field.ideal.FunctionFieldIdeal_global'>
        """
        FunctionFieldIdeal.__init__(self, ring)

        # the denominator and the entries of the hnf are
        # univariate polynomials.
        self._hnf = hnf
        self._denominator = denominator

        self._two_gens = None
        self._is_prime = None

        # for prime ideals
        self._relative_degree = None
        self._ramification_index = None
        self._prime_below = None
        self._beta = None

        # beta in matrix form for fast multiplicaton
        self._beta_matrix = None

        # (p, q) with irreducibl polynomial p and q an element of O in vector
        # form, together generating the prime ideal. This data is obtained by
        # Kummer's theorem when this prime ideal is constructed. This is used
        # for fast multiplication with other ideal.
        self._kummer_form = None

    def __nonzero__(self):
        """
        Test if the ideal is zero.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3*y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y); I
            Ideal (x, y) of Maximal order of Function field in y defined by y^2 + x^3*y + x
            sage: I.is_zero()
            False
            sage: J = 0*I; J
            Zero ideal of Maximal order of Function field in y defined by y^2 + x^3*y + x
            sage: J.is_zero()
            True

            sage: K.<x>=FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y>=K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y); I
            Ideal (x^2 + 1, y) of Maximal order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
            sage: I.is_zero()
            False
            sage: J = 0*I; J
            Zero ideal of Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: J.is_zero()
            True
        """
        return self._hnf.nrows() != 0

    def __hash__(self):
        """
        Return hash of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: { I: 2 }[I] == 2
            True

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: { I: 2 }[I] == 2
            True
        """
        return hash((self._ring, self._hnf, self._denominator))

    def __contains__(self, x):
        """
        Return True if ``x`` is in the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal([y]); I
            Ideal (x^3 + 1, y) of Maximal order of Function field in y
            defined by y^2 + 6*x^3 + 6
            sage: x * y in I
            True
            sage: y / x in I
            False
            sage: y^2 - 2 in I
            False

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal([y]); I
            Ideal (x^2 + 1, y) of Maximal order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
            sage: x * y in I
            True
            sage: y / x in I
            False
            sage: y^2 - 2 in I
            False
        """
        vec = self.ring().coordinate_vector(self._denominator * x)
        v = []
        for e in vec:
            if e.denominator() != 1:
                return False
            v.append(e.numerator())
        vec = vector(v)
        return vec in self._hnf.image()

    def __invert__(self):
        """
        Return the inverse fractional ideal of the ideal.

        EXAMPLES::
            sage: K.<x> = FunctionField(GF(7)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: ~I
            Ideal (1, (1/(x^3 + 1))*y) of Maximal order of Function field in y defined by y^2 + 6*x^3 + 6
            sage: I^(-1)
            Ideal (1, (1/(x^3 + 1))*y) of Maximal order of Function field in y defined by y^2 + 6*x^3 + 6
            sage: ~I * I
            Ideal (1) of Maximal order of Function field in y defined by y^2 + 6*x^3 + 6

            sage: K.<x> = FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: ~I
            Ideal (x, (x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: I^(-1)
            Ideal (x, (x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: ~I * I
            Ideal (1) of Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        R = self.ring()
        T = R._codifferent_matrix()
        I = self * R.codifferent()
        J = I._denominator * (I._hnf * T).inverse()
        return R._ideal_from_vectors(J.columns())

    def _richcmp_(self, other, op):
        """
        Compare the ideal with the other ideal with respect to ``op``.

        INPUT:

        - ``other`` -- ideal

        - ``op`` -- comparison operator

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: I == I + I
            True
            sage: I == I * I
            False

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: I == I + I
            True
            sage: I == I * I
            False
            sage: I < I * I
            True
            sage: I > I * I
            False
        """
        return richcmp((self._denominator, self._hnf), (other._denominator, other._hnf), op)

    def _repr_(self):
        """
        Return string representation.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y); I
            Ideal (x, y) of Maximal order of Function field in y defined by y^2 + x^3*y + x

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y); I
            Ideal (x^2 + 1, y) of Maximal order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
        """
        if self.is_zero():
            return "Zero ideal of %s" % (self._ring,)

        s = ', '.join([repr(g) for g in self.gens_two()])
        return "Ideal (%s) of %s"%(s, self._ring)

    def _add_(self, other):
        """
        Add with other ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I + J
            Ideal (x, y) of Maximal order of Function field in y defined by y^2 + x^3*y + x

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I + J
            Ideal (1, y) of Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        ds = self._denominator
        do = other._denominator
        vecs1 = [do * r for r in self._hnf]
        vecs2 = [ds * r for r in other._hnf]
        return self._ring._ideal_from_vectors_and_denominator(vecs1 + vecs2, ds * do)

    def _mul_(self, other):
        """
        Multiply with other ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I * J
            Ideal (x^4 + x^2 + x, x*y + x^4 + x) of Maximal order
            of Function field in y defined by y^2 + x^3*y + x

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I * J
            Ideal ((x^4 + x^3 + x + 1)/x, (x + 1)*y + (x^2 + 1)/x) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        O = self._ring
        mul = O._mul_vecs

        if self._kummer_form is not None:
            p, q = self._kummer_form
            vecs = list(p * other._hnf) + [mul(q, v) for v in other._hnf]
        elif other._kummer_form is not None:
            p, q = other._kummer_form
            vecs = list(p * self._hnf) + [mul(q, v) for v in self._hnf]
        else:
            vecs = [mul(r1,r2) for r1 in self._hnf for r2 in other._hnf]

        return O._ideal_from_vectors_and_denominator(vecs, self._denominator * other._denominator)

    def _rmul_(self, other):
        """
        Multiply other and the ideal on the right.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: J = O.ideal(x)
            sage: x * I ==  I * J
            True

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: J = O.ideal(x)
            sage: x * I ==  I * J
            True
        """
        O = self._ring
        mul = O._mul_vecs

        # compute the vector form of other element
        v = O.coordinate_vector(other)
        d = v.denominator()
        vec = vector([(d * c).numerator() for c in v])

        # multiply with the ideal
        vecs = [mul(vec, r) for r in self._hnf]

        return O._ideal_from_vectors_and_denominator(vecs, d * self._denominator)

    def intersect(self, other):
        """
        Intersect the ideal with the other ideal as fractional ideals.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: J = O.ideal(x)
            sage: I.intersect(J) == I * J * (I + J)^-1
            True

        """
        from sage.matrix.special import block_matrix

        A = self._hnf
        B = other._hnf

        ds = self.denominator()
        do = other.denominator()
        d = ds.lcm(do)
        if not d.is_one():
            A = (d // ds) * A
            B = (d // do) * B

        MS = A.matrix_space()
        I = MS.identity_matrix()
        O = MS.zero()
        n = A.ncols()

        # intersect the row spaces of A and B
        M = block_matrix([[I,I],[A,O],[O,B]])
        _,U = M.hermite_form_reversed(transformation=True)
        vecs = [U[i][:n] for i in range(n)]

        return self._ring._ideal_from_vectors_and_denominator(vecs, d)

    def hnf(self):
        """
        Return the matrix in hermite normal form representing the ideal.

        See also :meth:`denominator`

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y*(y+1)); I.hnf()
            [x^6 + x^3         0]
            [  x^3 + 1         1]
        """
        return self._hnf

    def denominator(self):
        """
        Return the denominator of the fractional ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y/(y+1))
            sage: d = I.denominator(); d
            x^3
            sage: d in O
            True
        """
        return self._denominator

    @cached_method
    def module(self):
        """
        Return the module, that is the ideal viewed as a module
        over the base maximal order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: F.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x,1/y)
            sage: I.module()
            Free module of degree 2 and rank 2 over Maximal order
            of Rational function field in x over Finite Field of size 7
            Echelon basis matrix:
            [          1           0]
            [          0 1/(x^3 + 1)]
        """
        F = self.ring().fraction_field()
        V, fr, to = F.vector_space()
        O = F.base_field().maximal_order()
        return V.span([to(g) for g in self.gens_over_base()], base_ring=O)

    @cached_method
    def gens_over_base(self):
        """
        Return the generators of the ideal as a module over the
        maximal order of the base rational function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: I.gens_over_base()
            (x^4 + x^2 + x, y + x)

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: I.gens_over_base()
            (x^3 + 1, y + x)
        """
        gens, d  = self._gens_over_base()
        return tuple([~d * b for b in gens])

    @cached_method
    def _gens_over_base(self):
        """
        Return the generators of the integral ideal, which is
        the denominator times of the fractional ideal, together
        with the denominator.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3*y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: I._gens_over_base()
            ([x, y], x)
        """
        gens = []
        for row in self._hnf:
            gens.append(sum([c1*c2 for c1,c2 in zip(row, self._ring.basis())]))
        return gens, self._denominator

    def gens(self):
        """
        Return a set of generators of the ideal.

        This provides whatever set of generators as quickly
        as possible.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: I.gens()
            (x^4 + x^2 + x, y + x)

            sage: L.<y> = K.extension(Y^2 +Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x+y)
            sage: I.gens()
            (x^3 + 1, y + x)
        """
        if self._two_gens is None:
            return self.gens_over_base()
        return self._two_gens

    def gens_two(self):
        """
        Return two generators of the fractional ideal.

        If the ideal is principal, one generator may be returned.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I  # indirect doctest
            Ideal (x^3 + x^2 + x, y) of Maximal order of Function field
            in y defined by y^3 + x^6 + x^4 + x^2
            sage: ~I  # indirect doctest
            Ideal (1, (1/(x^6 + x^4 + x^2))*y^2) of Maximal order of Function field
            in y defined by y^3 + x^6 + x^4 + x^2

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I  # indirect doctest
            Ideal (x^2 + 1, y) of Maximal order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
            sage: ~I  # indirect doctest
            Ideal (x, (x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        if self._two_gens is None:
            if self._kummer_form is not None: # prime ideal
                O = self.ring()
                _g1,_g2 = self._kummer_form
                g1 = O.fraction_field()(_g1)
                g2 = sum([c1*c2 for c1,c2 in zip(_g2, O.basis())])
                self._two_gens = (g1,g2)
            else:
                self._two_gens = self._gens_two()
        d = self.denominator()
        return tuple(e/d for e in self._two_gens)

    def _gens_two(self):
        """
        Return a set of two generators of the integral ideal, that is
        the denominator times the fractional ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x^2,x*y,x+y)
            sage: I._gens_two()
            (x, y)
        """
        O = self.ring()
        F = O.fraction_field()
        hnf = self._hnf

        norm = reduce(operator.mul, hnf.diagonal())

        if norm.is_constant(): # unit ideal
            return (F(1),)

        # one generator; see .ideal_below()
        l = hnf[0][0]
        p = l.degree()
        l = F(l)

        if self._hnf == O.ideal(l)._hnf: # principal ideal
            return (l,)

        R = hnf.base_ring()

        basis = []
        for row in hnf:
            basis.append(sum([c1*c2 for c1,c2 in zip(row, O.basis())]))

        n = len(basis)
        alpha = None; found = False

        def check(alpha):
            alpha_norm = alpha.norm().numerator() # denominator is 1
            found = norm.gcd(alpha_norm // norm) == 1
            return found

        # Trial 1: search for alpha among generators
        for alpha in basis:
            if check(alpha):
                return (l, alpha)

        # Trial 2: exhaustive search for alpha using only polynomials
        # with coefficients 0 or 1
        for d in range(p):
            G = itertools.product(itertools.product([0,1],repeat=d+1), repeat=n)
            for g in G:
                alpha = sum([R(c1)*c2 for c1,c2 in zip(g, basis)])
                if check(alpha):
                    return (l, alpha)

        # Trial 3: exhaustive search for alpha using all polynomials
        for d in range(p):
            G = itertools.product(R.polynomials(max_degree=d), repeat=n)
            for g in G:
                # discard duplicate cases
                if max(c.degree() for c in g) != d: continue
                for j in range(n):
                    if g[j] != 0: break
                if g[j].leading_coefficient() != 1: continue

                alpha = sum([c1*c2 for c1,c2 in zip(g, basis)])
                if check(alpha):
                    return (l, alpha)

    @cached_method
    def basis_matrix(self):
        """
        Return the matrix of basis vectors of the ideal as a module.

        The basis matrix is by definition the hermite norm form of the ideal
        divided by the denominator.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x,1/y)
            sage: I.denominator() * I.basis_matrix() == I.hnf()
            True
        """
        O = self.ring()
        d = self.denominator()
        m = (d * self).hnf()
        if d != 1:
            m = ~d * m
        return m

    def is_integral(self):
        """
        Return True if this is an integral ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x,1/y)
            sage: I.is_integral()
            False
            sage: J = I.denominator() * I
            sage: J.is_integral()
            True

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x,1/y)
            sage: I.is_integral()
            False
            sage: J = I.denominator() * I
            sage: J.is_integral()
            True
        """
        return self.denominator() == 1

    def ideal_below(self):
        """
        Return the ideal below the ideal.

        This is defined only for integral ideals.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x,1/y)
            sage: I.ideal_below()
            Traceback (most recent call last):
            ...
            TypeError: not an integral ideal
            sage: J = I.denominator() * I
            sage: J.ideal_below()
            Ideal (x^3 + x^2 + x) of Maximal order of Rational function field
            in x over Finite Field of size 2

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x,1/y)
            sage: I.ideal_below()
            Traceback (most recent call last):
            ...
            TypeError: not an integral ideal
            sage: J = I.denominator() * I
            sage: J.ideal_below()
            Ideal (x^3 + x) of Maximal order of Rational function field
            in x over Finite Field of size 2
        """
        if not self.is_integral():
            raise TypeError("not an integral ideal")

        K = self.ring().fraction_field().base_field().maximal_order()

        # self._hnf is in reversed hermite normal form, that is, lower
        # triangular form. Thus the generator of the ideal below is
        # just the (0,0) entry of the normal form. When self._hnf was in
        # hermite normal form, that is, upper triangular form, then the
        # generator can be computed in the following way:
        #
        #   m = matrix([hnf[0].parent().gen(0)] + list(hnf))
        #   _,T = m.hermite_form(transformation=True)
        #   return T[-1][0]
        #
        # This is certainly less efficient! This is an argument for using
        # reversed hermite normal form for ideal representation.
        l = self._hnf[0][0]

        return K.ideal(l)

    def factor(self):
        """
        Return the factorization of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I == I.factor().prod()
            True

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I == I.factor().prod()
            True
        """
        return Factorization(self._factor(), cr=True)

    def _factor(self):
        """
        Return the factorization of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I == I.factor().prod()  # indirect doctest
            True
        """
        O = self.ring()
        F = O.fraction_field()
        o = F.base_field().maximal_order()

        # First we collect primes below self
        d = self._denominator
        i = d * self

        factors = []
        primes = set([o.ideal(p) for p,_ in d.factor()] + [p for p,_ in i.ideal_below().factor()])
        for prime in primes:
            qs = [q[0] for q in O.decomposition(prime)]
            for q in qs:
                exp = q.valuation(self)
                if exp != 0:
                    factors.append((q,exp))
        return factors

    def norm(self):
        """
        Return the norm of the fractional ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: i1 = O.ideal(x)
            sage: i2 = O.ideal(y)
            sage: i3 = i1 * i2
            sage: i3.norm() == i1.norm() * i2.norm()
            True
            sage: i1.norm()
            x^3
            sage: i1.norm() == x ** F.degree()
            True
            sage: i2.norm()
            x^6 + x^4 + x^2
            sage: i2.norm() == y.norm()
            True

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: i1 = O.ideal(x)
            sage: i2 = O.ideal(y)
            sage: i3 = i1 * i2
            sage: i3.norm() == i1.norm() * i2.norm()
            True
            sage: i1.norm()
            x^2
            sage: i1.norm() == x ** L.degree()
            True
            sage: i2.norm()
            (x^2 + 1)/x
            sage: i2.norm() == y.norm()
            True
        """
        return reduce(operator.mul, self.basis_matrix().diagonal())

    def divisor(self):
        """
        Return the divisor corresponding to the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I.divisor()
            2*Place (x, (1/(x^3 + x^2 + x))*y^2)
             + 2*Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2+Y+x+1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I.divisor()
            -1*Place (x, x*y)
             + 2*Place (x + 1, x*y)
        """
        data = {prime.place(): multiplicity for prime, multiplicity in self.factor()}
        return FunctionFieldDivisor(self.ring().fraction_field(), data)

    def divisor_of_zeros(self):
        """
        Return the divisor of zeros corresponding to the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2+Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I.divisor_of_zeros()
            2*Place (x + 1, x*y)
        """
        data = {prime.place(): multiplicity for prime, multiplicity in self.factor() if multiplicity > 0}
        return FunctionFieldDivisor(self.ring().fraction_field(), data)

    def divisor_of_poles(self):
        """
        Return the divisor of poles corresponding to the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I.divisor_of_poles()
            Place (x, x*y)
        """
        data = {prime.place(): - multiplicity for prime, multiplicity in self.factor() if multiplicity < 0}
        return FunctionFieldDivisor(self.ring().fraction_field(), data)

    def is_prime(self):
        """
        Return True if the ideal is a prime ideal.

        If checked to be a prime ideal, then the ideal can be used
        as a prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.is_prime() for f,_ in I.factor()]
            [True, True]

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.is_prime() for f,_ in I.factor()]
            [True, True]
        """
        if self._is_prime is None:
            factors = self.factor()
            self._is_prime = len(factors) == 1 and factors[0][1] == 1
            if self._is_prime:
                prime = factors[0][0]
                assert self == prime
                self._relative_degree = prime._relative_degree
                self._ramification_index = prime._ramification_index
                self._prime_below = prime._prime_below
                self._beta = prime._beta
                self._beta_matrix = prime._beta_matrix
        return self._is_prime

    ###################################################
    # The following methods are only for prime ideals #
    ###################################################

    def valuation(self, ideal):
        """
        Return the valuation of the ideal at the prime ideal.

        INPUT:

        - ``ideal`` -- fractional ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x, (1/(x^3 + x^2 + x))*y^2)
            sage: I.is_prime()
            True
            sage: J = O.ideal(y)
            sage: I.valuation(J)
            2

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.valuation(I) for f,_ in I.factor()]
            [-1, 2]

        The method closely follows Algorithm 4.8.17 of [Coh1993]_.
        """
        if ideal.is_zero(): return infinity

        O = self.ring()
        F = O.fraction_field()
        n = F.degree()

        # beta_matrix is for fast multiplication with beta. For performance,
        # this is computed here rather than when the prime is constructed.
        if self._beta_matrix is None:
            beta = self._beta
            m = []
            for i in range(n):
                mtable_row = O._mtable[i]
                c = sum(beta[j] * mtable_row[j] for j in range(n))
                m.append(c)
            self._beta_matrix = matrix(m)

        B = self._beta_matrix

        # Step 1: compute the valuation of the denominator times the ideal
        #
        # This part is highly optimized as it is critical for
        # overall performance of the function field machinery.
        p = self.prime_below().gen().numerator()
        h = ideal._hnf.list()
        val = min([c.valuation(p) for c in h])
        i = self._ramification_index * val
        while True:
            ppow = p ** val
            h = (matrix(n, [c // ppow for c in h]) * B).list()
            val = min([c.valuation(p) for c in h])
            if val.is_zero():
                break
            i += self._ramification_index * (val - 1) + 1

        # Step 2: compute the valuation of the denominator
        j = self._ramification_index * ideal.denominator().valuation(p)

        # Step 3: return the valuation
        return i - j

    def prime_below(self):
        """
        Return the prime below of the prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.prime_below() for f,_ in I.factor()]
            [Ideal (x) of Maximal order of Rational function field in x
            over Finite Field of size 2, Ideal (x^2 + x + 1) of Maximal order
            of Rational function field in x over Finite Field of size 2]

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.prime_below() for f,_ in I.factor()]
            [Ideal (x) of Maximal order of Rational function field in x over Finite Field of size 2,
             Ideal (x + 1) of Maximal order of Rational function field in x over Finite Field of size 2]
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        return self._prime_below

    def place(self):
        """
        Return the place corresponding to the prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.place() for f,_ in I.factor()]
            [Place (x, (1/(x^3 + x^2 + x))*y^2),
             Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)]

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.place() for f,_ in I.factor()]
            [Place (x, x*y), Place (x + 1, x*y)]
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        return FunctionFieldPlace_global(self.ring().fraction_field(), self)

class FunctionFieldIdealInfinite(FunctionFieldIdeal):
    """
    Base class of ideals of maximal infinite orders
    """
    pass

class FunctionFieldIdealInfinite_rational(FunctionFieldIdealInfinite):
    """
    Fractional ideal of the maximal order of rational function field.
    """
    def __init__(self, ring, gen):
        """
        Initialize.

        INPUT:

        - ``ring`` -- infinite maximal order

        - ``gen``-- generator

        Note that the infinite maximal order is a principal ideal domain.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: K.maximal_order_infinite()
            Maximal infinite order of Rational function field in x
            over Finite Field of size 2
        """
        FunctionFieldIdealInfinite.__init__(self, ring)
        self._gen = gen

    def __hash__(self):
        """
        Return hash of the fractional ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x)
            sage: J = Oinf.ideal(1/x)
            sage: d = { I: 1, J: 2 }
        """
        return hash( (self.ring(), self._gen) )

    def _richcmp_(self, other, op):
        """
        Compare the differential and ``other`` with respect to ``op``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x+1)
            sage: J = Oinf.ideal(x^2+x)
            sage: I + J == J
            True
        """
        return richcmp(self._gen, other._gen, op)

    def _repr_(self):
        """
        Return the string representation of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x/(x^2+1))
            sage: I
            Ideal (1/x) of Maximal infinite order of Rational function field
            in x over Finite Field of size 2
        """
        return "Ideal (%s) of %s"%(self._gen, self._ring)

    def _add_(self, other):
        """
        Add the ideal with the other ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x/(x^2+1))
            sage: J = Oinf.ideal(1/(x+1))
            sage: I + J
            Ideal (1/x) of Maximal infinite order of Rational function field
            in x over Finite Field of size 2
        """
        return self._ring.ideal([self._gen, other._gen])

    def _mul_(self, other):
        """
        Multiply the ideal with the other ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x/(x^2+1))
            sage: J = Oinf.ideal(1/(x+1))
            sage: I * J
            Ideal (1/x^2) of Maximal infinite order of Rational function field
            in x over Finite Field of size 2
        """
        return self._ring.ideal([self._gen * other._gen])

    def _rmul_(self, other):
        """
        Multiply the ideal with the other ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x/(x^2+1))
            sage: x * I
            Ideal (1) of Maximal infinite order of Rational function field
            in x over Finite Field of size 2
        """
        return self._ring.ideal([other * self._gen])

    def __invert__(self):
        """
        Return the multiplicative inverse of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x/(x^2 + 1))
            sage: ~I  # indirect doctest
            Ideal (x) of Maximal infinite order of Rational function field in x
            over Finite Field of size 2
        """
        return self._ring.ideal([~self._gen])

    def is_prime(self):
        """
        Return True if the ideal is a prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x/(x^2 + 1))
            sage: I.is_prime()
            True
        """
        x = self._ring.fraction_field().gen()
        return self._gen == 1/x

    def gen(self):
        """
        Return the generator of the principal ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x+1)/(x^3+x),(x^2+1)/x^4)
            sage: I.gen()
            1/x^2
        """
        return self._gen

    def gens(self):
        """
        Return the generator of the principal ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x+1)/(x^3+x),(x^2+1)/x^4)
            sage: I.gens()
            (1/x^2,)
        """
        return (self._gen,)

    def gens_over_base(self):
        """
        Return the generator of the ideal as a rank one module
        over the infinite maximal order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x+1)/(x^3+x),(x^2+1)/x^4)
            sage: I.gens_over_base()
            (1/x^2,)
        """
        return (self._gen,)

    def factor(self):
        """
        Return the factorization of the ideal into prime ideals.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x+1)/(x^3+1))
            sage: I.factor()
            (Ideal (1/x) of Maximal infinite order of Rational function field
            in x over Finite Field of size 2)^2
        """
        g = ~(self.ring().fraction_field().gen())
        m = self._gen.denominator().degree() - self._gen.numerator().degree()
        if m == 0:
            factors = []
        else:
            factors = [(self.ring().ideal(g), m)]
        return Factorization(factors, cr=True)

    def divisor(self):
        """
        Return divisor corresponding to the fractional ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x+1)/(x^3+1))
            sage: I.divisor()
            2*Place (1/x)
        """
        if self._gen == 0:
            raise ValueError("not defined for zero element")

        data = {prime.place(): multiplicity for prime, multiplicity in self.factor()}
        return FunctionFieldDivisor(self.ring().fraction_field(), data)

    def divisor_of_zeros(self):
        """
        Return divisor corresponding to the fractional ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x+1)/(x^3+1))
            sage: I.divisor_of_zeros()
            2*Place (1/x)
        """
        if self._gen == 0:
            raise ValueError("not defined for zero element")

        data = {prime.place(): multiplicity for prime, multiplicity in self.factor() if multiplicity > 0}
        return FunctionFieldDivisor(self.ring().fraction_field(), data)

    def divisor_of_poles(self):
        """
        Return divisor corresponding to the fractional ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x + 1)/(x^3 + 1))
            sage: I.divisor_of_poles()
            0
        """
        if self._gen == 0:
            raise ValueError("not defined for zero element")

        data = {prime.place(): - multiplicity for prime, multiplicity in self.factor() if multiplicity < 0}
        return FunctionFieldDivisor(self.ring().fraction_field(), data)

    def place(self):
        """
        Return the place corresponding to the prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x+1)/(x^3+1))
            sage: p = I.factor()[0][0]
            sage: p.place()
            Place (1/x)
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        return FunctionFieldPlace_rational(self.ring().fraction_field(), self)

    def valuation(self, ideal):
        """
        Return valuation value of the ideal at the prime ideal.

        INPUT:

        - ``ideal`` -- fractional ideal

        EXAMPLES::

            sage: F.<x> = FunctionField(QQ)
            sage: O = F.maximal_order_infinite()
            sage: p = O.ideal(1/x)
            sage: p.valuation(O.ideal(x/(x+1)))
            0
            sage: p.valuation(O.ideal(0))
            +Infinity
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        f = ideal.gen()
        if f == 0:
            return infinity
        else:
            return f.denominator().degree() - f.numerator().degree()

class FunctionFieldIdealInfinite_module(FunctionFieldIdealInfinite):
    """
    A fractional ideal specified by a finitely generated module over
    the integers of the base field.
    """
    def __init__(self, ring, module):
        """
        Initialize.

        INPUT:

        - ``ring`` -- order in a function field

        - ``module`` -- module

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(y)
            sage: type(I)
            <class 'sage.rings.function_field.ideal.FunctionFieldIdeal_module'>
        """
        FunctionFieldIdealInfinite.__init__(self, ring)

        self._module = module
        self._structure = ring.fraction_field().vector_space()

        V, from_V, to_V = self._structure
        gens = tuple([from_V(a) for a in module.basis()])
        self._gens = gens

    def __contains__(self, x):
        """
        Return True if ``x`` is in the ideal.

        INPUT:

        - ``x`` -- element of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal_with_gens_over_base([1, y]);  I
            Ideal (1, y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: y in I
            True
            sage: y/x in I
            False
            sage: y^2 - 2 in I
            True
        """
        return self._structure[2](x) in self._module

    def __hash__(self):
        """
        Return the hash of the ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal_with_gens_over_base([1, y])
            sage: d = {I: 2}  # indirect doctest
        """
        return hash((self._ring,self._module))

    def __eq__(self, other):
        """
        Test equality of the ideal with the other ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal_with_gens_over_base([1, y])
            sage: I == I + I  # indirect doctest
            True
        """
        if not isinstance(other, FunctionFieldIdeal_module):
            other = self.ring().ideal(other)
        if self.ring() != other.ring():
            raise ValueError("rings must be the same")

        if (self.module().is_submodule(other.module()) and
            other.module().is_submodule(self.module())):
            return True
        else:
            return False

    def _repr_short(self):
        """
        Represent the list of generators.

        EXAMPLE::

            sage: P.<a,b,c> = QQ[]
            sage: P*[a^2,a*b+c,c^3]
            Ideal (a^2, a*b + c, c^3) of Multivariate Polynomial Ring in a, b, c over Rational Field
            sage: (P*[a^2,a*b+c,c^3])._repr_short()
            '(a^2, a*b + c, c^3)'

        If the string representation of a generator contains a line break,
        the generators are not represented from left to right but from
        top to bottom. This is the case, e.g., for matrices::

            sage: MS = MatrixSpace(QQ,2,2)
            sage: MS*[MS.1,2]
            Left Ideal
            (
              [0 1]
              [0 0],
            <BLANKLINE>
              [2 0]
              [0 2]
            )
             of Full MatrixSpace of 2 by 2 dense matrices over Rational Field


        """
        L = []
        has_return = False
        for x in self.gens():
            s = repr(x)
            if '\n' in s:
                has_return = True
                s = s.replace('\n','\n  ')
            L.append(s)
        if has_return:
            return '\n(\n  %s\n)\n'%(',\n\n  '.join(L))
        return '(%s)'%(', '.join(L))

    def _repr_(self):
        """
        Return a string representation of the ideal.

        EXAMPLES::

            sage: P.<a,b,c> = QQ[]
            sage: P*[a^2,a*b+c,c^3] # indirect doctest
            Ideal (a^2, a*b + c, c^3) of Multivariate Polynomial Ring in a, b, c over Rational Field
        """
        return "Ideal %s of %s"%(self._repr_short(), self.ring())

    def module(self):
        """
        Return module over the maximal order of the base field that
        underlies self.

        The formation of the module is compatible with the vector
        space corresponding to the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: O = K.maximal_order(); O
            Maximal order of Rational function field in x over Finite Field of size 7
            sage: K.polynomial_ring()
            Univariate Polynomial Ring in x over Rational function field in x over Finite Field of size 7
            sage: I = O.ideal([x^2 + 1, x*(x^2+1)])
            sage: I.gens()
            (x^2 + 1,)
            sage: I.module()
            Free module of degree 1 and rank 1 over Maximal order of Rational function field in x over Finite Field of size 7
            Echelon basis matrix:
            [x^2 + 1]
            sage: V, from_V, to_V = K.vector_space(); V
            Vector space of dimension 1 over Rational function field in x over Finite Field of size 7
            sage: I.module().is_submodule(V)
            True
        """
        return self._module

class FunctionFieldIdealInfinite_global(FunctionFieldIdealInfinite):
    """
    Ideals of the infinite maximal order.
    """
    def __init__(self, ring, ideal):
        """
        Initialize the ideal.

        INPUT:

        - ``ring`` -- infinite maximal order of the function field

        - ``ideal`` -- ideal in the inverted function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3+t^2-x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: Oinf.ideal(1/y)
            Ideal (1/x^2,1/x^4*y^2) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4
        """
        FunctionFieldIdealInfinite.__init__(self, ring)
        self._ideal = ideal

    def __hash__(self):
        """
        Return hash of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/y)
            sage: d = { I: 1 }

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/y)
            sage: d = { I: 1 }
        """
        return hash((self.ring(), self._ideal))

    def _repr_(self):
        """
        Return a string representation of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: Oinf.ideal(1/y)
            Ideal (1/x^2,1/x^4*y^2) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.ideal(1/y)
            Ideal (1/x,1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
        """
        return self._repr_short() + " of {}".format(self.ring())

    def _repr_short(self):
        """
        Return a string representation of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = K[]
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/y)
            sage: I._repr_short()
            'Ideal (1/x^2,1/x^4*y^2)'
        """
        gens_str = ','.join([repr(g) for g in self.gens_two()])
        return "Ideal ({})".format(gens_str)

    def _add_(self, other):
        """
        Add the ideal with the other ideal.

        INPUT:

        - ``ideal`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3+t^2-x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x^2*1/y)
            sage: J = Oinf.ideal(1/x)
            sage: I + J
            Ideal (1/x) of Maximal infinite order of Function field in y
            defined by y^3 + y^2 + 2*x^4

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x^2*1/y)
            sage: J = Oinf.ideal(1/x)
            sage: I + J
            Ideal (1/x) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
        """
        return FunctionFieldIdealInfinite_global(self._ring, self._ideal + other._ideal)

    def _mul_(self, other):
        """
        Multiply the ideal with the other ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3+t^2-x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x^2*1/y)
            sage: J = Oinf.ideal(1/x)
            sage: I * J
            Ideal (1/x^5,1/x^7*y^2) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x^2*1/y)
            sage: J = Oinf.ideal(1/x)
            sage: I * J
            Ideal (1/x^4,1/x^4*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
        """
        return FunctionFieldIdealInfinite_global(self._ring, self._ideal * other._ideal)

    def __pow__(self, n):
        """
        Raise the ideal to ``n``-th power.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3+t^2-x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: J = Oinf.ideal(1/x)
            sage: J^3
            Ideal (1/x^3) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4
        """
        return FunctionFieldIdealInfinite_global(self._ring, self._ideal ** n)

    def __invert__(self):
        """
        Return the inverted ideal of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: J = Oinf.ideal(y)
            sage: ~J
            Ideal (1/x^2,1/x^4*y^2) of Maximal infinite order
            of Function field in y defined by y^3 + y^2 + 2*x^4
            sage: J * ~J
            Ideal (1) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: J = Oinf.ideal(y)
            sage: ~J
            Ideal (1/x,1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
            sage: J * ~J
            Ideal (1) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
        """
        return FunctionFieldIdealInfinite_global(self._ring, ~ self._ideal)

    def _richcmp_(self, other, op):
        """
        Compare the ideal with the other ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3+t^2-x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x^2*1/y)
            sage: J = Oinf.ideal(1/x)
            sage: I * J == J * I
            True
            sage: I + J == J
            True
            sage: I + J == I
            False
            sage: (I < J) == (not J < I)
            True

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x^2*1/y)
            sage: J = Oinf.ideal(1/x)
            sage: I * J == J * I
            True
            sage: I + J == J
            True
            sage: I + J == I
            False
            sage: (I < J) == (not J < I)
            True
        """
        return richcmp(self._ideal, other._ideal, op)

    @property
    def _relative_degree(self):
        """
        Return the relative degree of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: [J._relative_degree for J,_ in I.factor()]
            [1]
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        return self._ideal._relative_degree

    def gens(self):
        """
        Return a set of generators of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3+t^2-x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(x+y)
            sage: I.gens()
            (x, y, 1/x^2*y^2)

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(x+y)
            sage: I.gens()
            (x, y)
        """
        F = self.ring().fraction_field()
        iF,from_iF,to_iF = F._inversion_isomorphism()
        return tuple(from_iF(b) for b in self._ideal.gens())

    def gens_two(self):
        """
        Return a set of at most two generators of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3+t^2-x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(x+y)
            sage: I.gens_two()
            (x, y)

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2+Y+x+1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(x+y)
            sage: I.gens_two()
            (x,)
        """
        F = self.ring().fraction_field()
        iF,from_iF,to_iF = F._inversion_isomorphism()
        return tuple(from_iF(b) for b in self._ideal.gens_two())

    def gens_over_base(self):
        """
        Return a set of generators of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(x + y)
            sage: I.gens_over_base()
            (x, y, 1/x^2*y^2)
        """
        F = self.ring().fraction_field()
        iF,from_iF,to_iF = F._inversion_isomorphism()
        return tuple(from_iF(g) for g in self._ideal.gens_over_base())

    def ideal_below(self):
        """
        Return a set of generators of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/y^2)
            sage: I.ideal_below()
            Ideal (x^3) of Maximal order of Rational function field
            in x over Finite Field in z2 of size 3^2
        """
        return self._ideal.ideal_below()

    def is_prime(self):
        """
        Return True if the ideal is a prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: I.factor()
            (Ideal (1/x,1/x^3*y^2) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4)^3
            sage: I.is_prime()
            False
            sage: J = I.factor()[0][0]
            sage: J.is_prime()
            True

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: I.factor()
            (Ideal (1/x,1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x)^2
            sage: I.is_prime()
            False
            sage: J = I.factor()[0][0]
            sage: J.is_prime()
            True
        """
        return self._ideal.is_prime()

    @cached_method
    def prime_below(self):
        """
        Return the prime of the base order that underlies the prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3+t^2-x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: I.factor()
            (Ideal (1/x,1/x^3*y^2) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4)^3
            sage: J = I.factor()[0][0]
            sage: J.is_prime()
            True
            sage: J.prime_below()
            Ideal (1/x) of Maximal infinite order of Rational function field
            in x over Finite Field in z2 of size 3^2

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: I.factor()
            (Ideal (1/x,1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x)^2
            sage: J = I.factor()[0][0]
            sage: J.is_prime()
            True
            sage: J.prime_below()
            Ideal (1/x) of Maximal infinite order of Rational function field in x
            over Finite Field of size 2
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        F = self.ring().fraction_field()
        K = F.base_field()
        return K.maximal_order_infinite().prime_ideal()

    def factor(self):
        """
        Return factorization of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: Oinf = F.maximal_order_infinite()
            sage: f= 1/x
            sage: I = Oinf.ideal(f)
            sage: I.factor()
            (Ideal (1/x,1/x^4*y^2 + 1/x^2*y + 1) of Maximal infinite order
            of Function field in y defined by y^3 + x^6 + x^4 + x^2) *
            (Ideal (1/x,1/x^2*y + 1) of Maximal infinite order
            of Function field in y defined by y^3 + x^6 + x^4 + x^2)

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2+Y+x+1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: f= 1/x
            sage: I = Oinf.ideal(f)
            sage: I.factor()
            (Ideal (1/x,1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x)^2
        """
        return Factorization(self._factor(), cr=True)

    def _factor(self):
        """
        Return factorization of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)
            sage: Oinf = F.maximal_order_infinite()
            sage: f= 1/x
            sage: I = Oinf.ideal(f)
            sage: I.factor()  # indirect doctest
            (Ideal (1/x,1/x^4*y^2 + 1/x^2*y + 1) of Maximal infinite order
            of Function field in y defined by y^3 + x^6 + x^4 + x^2) *
            (Ideal (1/x,1/x^2*y + 1) of Maximal infinite order
            of Function field in y defined by y^3 + x^6 + x^4 + x^2)
        """
        O = self.ring()

        factors = []
        for iprime, exp in O._to_iF(self).factor():
            prime = FunctionFieldIdealInfinite_global(O, iprime)
            factors.append((prime, exp))

        return factors

    def place(self):
        """
        Return the place corresponding to the prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: I.factor()
            (Ideal (1/x,1/x^3*y^2) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4)^3
            sage: J = I.factor()[0][0]
            sage: J.is_prime()
            True
            sage: J.place()
            Place (1/x, 1/x^3*y^2)

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: I.factor()
            (Ideal (1/x,1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x)^2
            sage: J = I.factor()[0][0]
            sage: J.is_prime()
            True
            sage: J.place()
            Place (1/x, 1/x*y)
        """
        if not self.is_prime():
            raise ValueError("not a prime ideal")

        return FunctionFieldPlace_global(self.ring().fraction_field(), self)

    def divisor(self):
        """
        Return the divisor corresponding to the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(y)
            sage: I.divisor()
            -2*Place (1/x, 1/x^4*y^2 + 1/x^2*y + 1)
             - 2*Place (1/x, 1/x^2*y + 1)

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(y)
            sage: I.divisor()
            -1*Place (1/x, 1/x*y)
        """
        F = self.ring().fraction_field()

        if self._ideal._is_prime:
            from .divisor import prime_divisor

            place = FunctionFieldPlace_global(F, self)
            div = prime_divisor(F, place)
            return div

        data = {prime.place(): multiplicity for prime, multiplicity in self.factor()}
        return FunctionFieldDivisor(self.ring().fraction_field(), data)

    def divisor_of_zeros(self):
        """
        Return the divisor of zeros corresponding to the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2+Y+x+1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I.divisor_of_zeros()
            2*Place (x + 1, x*y)
        """
        data = {prime.place(): multiplicity for prime, multiplicity in self.factor() if multiplicity > 0}
        return FunctionFieldDivisor(self.ring().fraction_field(), data)

    def divisor_of_poles(self):
        """
        Return the divisor of poles corresponding to the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I.divisor_of_poles()
            Place (x, x*y)
        """
        data = {prime.place(): - multiplicity for prime, multiplicity in self.factor() if multiplicity < 0}
        return FunctionFieldDivisor(self.ring().fraction_field(), data)


    def valuation(self, ideal):
        """
        Return the valuation of the ideal with respect to the prime ideal.

        INPUT:

        - ``ideal`` -- fractional ideal

        EXAMPLES::

            sage: K.<x>=FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y>=K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(y)
            sage: [f.valuation(I) for f,_ in I.factor()]
            [-1]
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        return self._ideal.valuation(self.ring()._to_iF(ideal))

class IdealMonoid(Parent):
    r"""
    The monoid of ideals in orders of function fields.

    TESTS::

        sage: K.<x> = FunctionField(GF(2))
        sage: O = K.maximal_order()
        sage: M = O.ideal_monoid()
        sage: # TestSuite(M).run()
    """
    def __init__(self, R):
        """
        Initialize the ideal monoid.

        INPUT:

        - ``R`` -- order

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: O = K.maximal_order()
            sage: M = O.ideal_monoid(); M
            Monoid of ideals of Maximal order of Rational function field in x over Finite Field of size 2
        """
        Parent.__init__(self, base = R.fraction_field(), category = Monoids())
        self.__R = R
        self._populate_coercion_lists_()

    def _repr_(self):
        """
        Return a string representation of the ideal monoid.

        TESTS::

            sage: K.<x> = FunctionField(GF(2))
            sage: O = K.maximal_order()
            sage: M = O.ideal_monoid(); M._repr_()
            'Monoid of ideals of Maximal order of Rational function field in x over Finite Field of size 2'
        """
        return "Monoid of ideals of %s"%self.__R

    def ring(self):
        """
        Return the ring of which this is the ideal monoid.

        EXAMPLE::

            sage: K.<x> = FunctionField(GF(2))
            sage: O = K.maximal_order()
            sage: M = O.ideal_monoid(); M.ring() is O
            True
        """
        return self.__R

    def _element_constructor_(self, x):
        """
        Create an ideal in the monoid from ``x``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: O = K.maximal_order()
            sage: M = O.ideal_monoid()
            sage: M(x)
            Ideal (x) of Maximal order of Rational function field in x over Finite Field of size 2
            sage: M([x-4, 1/x])
            Ideal (1/x) of Maximal order of Rational function field in x over Finite Field of size 2
        """
        try: # x is an ideal
            x = x.gens()
        except AttributeError:
            pass
        return self.__R.ideal(x)

    def _coerce_map_from_(self, x):
        """
        Used by coercion framework.

        EXAMPLE::

            sage: K.<x> = FunctionField(GF(2))
            sage: O = K.maximal_order()
            sage: M = O.ideal_monoid()
            sage: M.has_coerce_map_from(O) # indirect doctest
            True
            sage: M.has_coerce_map_from(O.ideal_monoid())
            True
        """
        if isinstance(x, IdealMonoid):
            return self.ring().has_coerce_map_from(x.ring())
        else:
            return self.ring().has_coerce_map_from(x)

    def _an_element_(self):
        """
        Return an element of the ideal monoid.

        EXAMPLE::

            sage: K.<x> = FunctionField(GF(2))
            sage: O = K.maximal_order()
            sage: M = O.ideal_monoid()
            sage: M.an_element() # indirect doctest; random
            Ideal (x) of Maximal order of Rational function field in x
            over Finite Field of size 2
        """
        x = self.__R.an_element()
        return self.__R.ideal([x])
