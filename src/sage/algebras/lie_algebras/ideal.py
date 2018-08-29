r"""
Ideals of Lie algebras

AUTHORS:

- Eero Hakavuori (2018-08-29): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Eero Hakavuori <eero.hakavuori@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                 https://www.gnu.org/licenses/
# ****************************************************************************

from sage.algebras.lie_algebras.subalgebra import LieSubalgebra_finite_dimensional_with_basis
from sage.categories.lie_algebras import LieAlgebras
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.sets.family import Family


class LieIdeal_finite_dimensional_with_basis(LieSubalgebra_finite_dimensional_with_basis):
    r"""
    An ideal of a finite dimensional Lie algebra with basis.

    INPUT:

    - ``ambient`` -- the Lie algebra containing the ideal
    - ``gens`` -- a list of generators of the ideal
    - ``category`` -- (optional) a subcategory of subobjects of finite
      dimensional Lie algebras with basis

    EXAMPLES:

    An ideal is defined by giving a list of generators::

        sage: L = lie_algebras.Heisenberg(QQ, 1)
        sage: X, Y, Z = L.basis()
        sage: I =  L.ideal([X, Z]); I
        Ideal (p1, z) of Heisenberg algebra of rank 1 over Rational Field
        sage: I.basis()
        Family (p1, z)

    Different generating sets can lead to the same basis::

        sage: I2 = L.ideal(X); I2
        Ideal (p1) of Heisenberg algebra of rank 1 over Rational Field
        sage: I2.basis()
        Family (p1, z)

    Elements of the ambient Lie algebra can be reduced modulo the ideal::

        sage: I.reduce(X + 2*Y + 3*Z)
        2*q1

    The reduction gives elements in a fixed complementary subspace to the ideal.
    The complementary subspace is spanned by those basis elements which are not
    leading supports of the basis of the ideal::

        sage: L.<X,Y,Z> = LieAlgebra(SR, {('X','Y'): {'Z': 1}})
        sage: I =  L.ideal(X + Y)
        sage: I.basis()
        Family (X + Y, Z)
        sage: el = var('x')*X + var('y')*Y + var('z')*Z; el
        x*X + y*Y + z*Z
        sage: I.reduce(el)
        (x-y)*X

    It is possible to compute with ideals of ideals::

        sage: sc = {('X','Y'): {'Z': 1}, ('X','Z'): {'W': 1}}
        sage: L.<X,Y,Z,W> = LieAlgebra(QQ, sc)
        sage: I = L.ideal(Z); I
        Ideal (Z) of Lie algebra on 4 generators (X, Y, Z, W) over Rational Field
        sage: z, w = I.basis()
        sage: J = I.ideal(w); J
        Ideal (W) of Ideal (Z) of Lie algebra on 4 generators (X, Y, Z, W) over Rational Field
        sage: J.basis()
        Family (W,)
        sage: J.reduce(z + w)
        Z

    The zero ideal can be created by giving 0 as a generator or with an empty
    list of generators::

        sage: L.<X,Y,Z> = LieAlgebra(QQ, {('X','Y'): {'Z': 1}})
        sage: I1 = L.ideal(0)
        sage: I2 = L.ideal([])
        sage: I1 is I2
        True
        sage: I1.basis()
        Family ()

    TESTS:

    A test suite::

        sage: I =  L.ideal(X + Y)
        sage: TestSuite(I).run()
    """

    @staticmethod
    def __classcall_private__(cls, ambient, gens, category=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: L.<X,Y> = LieAlgebra(QQ, {('X','Y'): {'X': 1}})
            sage: I1 = L.ideal(X)
            sage: I2 = L.ideal((X,))
            sage: I3 = L.ideal([X])
            sage: I1 is I2 and I2 is I3
            True
        """
        if not isinstance(gens, (list, tuple)):
            gens = [gens]
        gens = tuple(ambient(gen) for gen in gens)

        gens = tuple(gens)
        if len(gens) == 0:
            gens = (ambient.zero(),)

        cat = LieAlgebras(ambient.base_ring()).FiniteDimensional().WithBasis()
        category = cat.Subobjects().or_subcategory(category)

        sup = super(LieIdeal_finite_dimensional_with_basis, cls)
        return sup.__classcall__(cls, ambient, gens, category)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L.<X,Y> = LieAlgebra(QQ, abelian=True)
            sage: L.ideal([X, Y])
            Ideal (X, Y) of Abelian Lie algebra on 2 generators (X, Y) over Rational Field
        """
        gens = self.gens()
        if len(gens) == 1:
            gens = "(%s)" % gens[0]
        return "Ideal %s of %s" % (gens, self.ambient())

    # for submodule computations, the order of the basis is reversed so that
    # the pivot elements in the echelon form are the leading terms
    def _to_m(self, X):
        r"""
        Return the reversed vector of an element of the ambient Lie algebra.

        INPUT:

        - ``X`` -- an element of the ambient Lie algebra

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z': 1}})
            sage: I =  L.ideal([x, z])
            sage: el = x + 2*y + 3*z
            sage: el.to_vector()
            (1, 2, 3)
            sage: I._to_m(el)
            (3, 2, 1)
        """
        return vector(self.ambient().base_ring(), reversed(X.to_vector()))

    def _from_m(self, v):
        r"""
        Return the element of the ambient Lie algebra from a reversed vector.

        INPUT:

        - ``v`` -- a vector

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z': 1}})
            sage: I =  L.ideal([x, z])
            sage: I._from_m([3, 2, 1])
            x + 2*y + 3*z
        """
        R = self.ambient().base_ring()
        return self.ambient().from_vector(vector(R, reversed(v)))

    @cached_method
    def basis(self):
        r"""
        Return a basis of ``self``.

        EXAMPLES:

            sage: sc = {('x','y'): {'z': 1}, ('x','z'): {'w': 1}}
            sage: L.<x,y,z,w> = LieAlgebra(QQ, sc)
            sage: L.ideal([x + y + z + w]).basis()
            Family (x + y, z, w)
        """
        L = self.ambient()
        m = L.module()
        B = m.basis()

        sm = m.submodule([self._to_m(X) for X in self.gens()])
        d = 0

        while sm.dimension() > d:
            d = sm.dimension()
            SB = sm.basis()
            sm = m.submodule(sm.basis() +
                             [self._to_m(L.bracket(v, self._from_m(w)))
                              for v in B for w in SB])

        return Family(reversed([self.element_class(self, self._from_m(v))
                                for v in sm.echelonized_basis()]))

    def reduce(self, X):
        r"""
        Reduce an element of the ambient Lie algebra modulo ``self``.

        INPUT:

        - ``X`` -- an element of the ambient Lie algebra

        OUTPUT:

        An element `Y` of the ambient Lie algebra that is contained in a fixed
        complementary submodule `V` to ``self`` such that `X = Y` mod ``self``.

        When the base ring of ``self`` is a field, the complementary submodule
        `V` is spanned by the elements of the basis that are not the leading
        supports of the basis of ``self``.

        EXAMPLES:

        An example reduction in a 6 dimensional Lie algebra::

            sage: sc = {('a','b'): {'d': 1}, ('a','c'): {'e': 1},
            ....:       ('b','c'): {'f': 1}}
            sage: L.<a,b,c,d,e,f> = LieAlgebra(QQ, sc)
            sage: I =  L.ideal(c)
            sage: I.reduce(a + b + c + d + e + f)
            a + b + d

        The reduction of an element is zero if and only if the element belongs
        to the ideal::

            sage: I.reduce(c + e)
            0
            sage: c + e in I
            True

        Over non-fields, the complementary submodule may not be spanned by
        a subset of the basis of the ambient Lie algebra::

            sage: L.<X,Y,Z> = LieAlgebra(ZZ, {('X','Y'): {'Z': 3}})
            sage: I = L.ideal(Y)
            sage: I.basis()
            Family (Y, 3*Z)
            sage: I.reduce(3*Z)
            0
            sage: I.reduce(Y + 14*Z)
            2*Z
        """
        R = self.base_ring()
        for Y in self.basis():
            Y = Y.value
            k, c = Y.leading_item()

            if R.is_field():
                X = X - X[k] / c * Y
            else:
                try:
                    q, r = X[k].quo_rem(c)
                    X = X - q * Y
                except AttributeError:
                    pass

        return X
