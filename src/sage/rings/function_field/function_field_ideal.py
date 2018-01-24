r"""
Ideals

Ideals of an order of a function field include all fractional ideals of the order.
Sage provides basic arithmetic with fractional ideals.

The fractional ideals of the maximal order of a global function field forms a multiplicative
monoid. Sage allows advanced arithmetic with the fractional ideals. For example, an ideal
of the maximal order can be factored into a product of prime ideals.

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

from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent
from sage.structure.element import AlgebraElement
from sage.structure.richcmp import richcmp
from sage.structure.factorization import Factorization

from sage.modules.free_module_element import vector

from sage.categories.monoids import Monoids

import sage.rings.integer_ring
from sage.rings.infinity import infinity

lazy_import('sage.matrix.constructor', 'matrix')

def is_Ideal(x):
    """
    Return ``True`` if ``x`` is an ideal in a function field.
    """
    return isinstance(x, FunctionFieldIdeal)

class FunctionFieldIdeal(AlgebraElement):
    """
    Fractional ideals of function fields.
    """
    def __init__(self, ring):
        """
        Initialize.
        """
        AlgebraElement.__init__(self, ring.ideal_monoid())
        self._ring = ring

    def gens_reduced(self):
        r"""
        Return reduced generators. This just returns the generators for now.

        This method is provided so that ideals in funtion fields have the method
        :meth:`gens_reduced()`, just like ideals of number fields. Sage linear algebra
        machinery sometimes requires this.

        """
        return self.gens()

class FunctionFieldIdeal_module(FunctionFieldIdeal):
    """
    A fractional ideal specified by a finitely generated module over
    the integers of the base field.

    EXAMPLES:

    An ideal in an extension of a rational function field::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x^3 - 1)
        sage: O = L.equation_order()
    """
    def __init__(self, ring, module):
        """
        Initialize.

        INPUT:

        - ``ring`` -- order in a function field

        - ``module`` -- module
        """
        FunctionFieldIdeal.__init__(self, ring)

        self._module = module
        self._structure = ring.fraction_field().vector_space()

        V, from_V, to_V = self._structure
        self._gens = tuple([from_V(a) for a in module.basis()])

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is in the ideal.
        """
        return self._structure[2](x) in self._module

    def __hash__(self):
        """
        Return the hash of the ideal.
        """
        return hash((self._ring,self._module))

    def _richcmp_(self, other, op):
        """
        Compare this ideal with the other ideal with respect to ``op``.
        """
        return richcmp(self.module(), other.module(), op)

    def module(self):
        """
        Return module over the maximal order of the base field that
        underlies the ideal.

        The formation of the module is compatible with the vector
        space corresponding to the function field.

        OUTPUT:

        - a module over the maximal order of the base field of the ideal

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

        EXAMPLES::

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

        EXAMPLES::

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

        EXAMPLES::

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

        EXAMPLES::

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

        EXAMPLES::

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

        The number of generators of ``self * other`` will be
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

        EXAMPLES::

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

