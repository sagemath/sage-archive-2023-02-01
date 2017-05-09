"""
Group algebras

This functionality has been moved to :mod:`sage.categories.group_algebras`.
"""

#*****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#                     2008 David Loeffler <d.loeffler.01@cantab.net>
#                     2009 Martin Raum <mraum@mpim-bonn.mpg.de>
#                     2011 John Palmieri <palmieri@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import IntegerRing
from sage.misc.cachefunc import cached_method
from sage.categories.pushout import ConstructionFunctor
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.all import Rings, Hom
from sage.categories.fields import Fields
from sage.categories.groups import Groups
from sage.categories.additive_groups import AdditiveGroups
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.morphism import SetMorphism
from sage.categories.sets_cat import Sets


import six

# TODO: Move to the categories/group_algebras.py file
class GroupAlgebraFunctor(ConstructionFunctor):
    r"""
    For a fixed group, a functor sending a commutative ring to the
    corresponding group algebra.

    INPUT:

    - ``group`` -- the group associated to each group algebra under
      consideration

    EXAMPLES::

        sage: from sage.algebras.group_algebra import GroupAlgebraFunctor
        sage: F = GroupAlgebraFunctor(KleinFourGroup())
        sage: loads(dumps(F)) == F
        True
        sage: GroupAlgebra(SU(2, GF(4, 'a')), IntegerModRing(12)).category()
        Category of finite group algebras over Ring of integers modulo 12
    """
    def __init__(self, group):
        r"""
        See :class:`GroupAlgebraFunctor` for full documentation.

        EXAMPLES::

            sage: from sage.algebras.group_algebra import GroupAlgebraFunctor
            sage: GroupAlgebra(SU(2, GF(4, 'a')), IntegerModRing(12)).category()
            Category of finite group algebras over Ring of integers modulo 12
        """
        self.__group = group

        ConstructionFunctor.__init__(self, Rings(), Rings())

    def group(self):
        r"""
        Return the group which is associated to this functor.

        EXAMPLES::

            sage: from sage.algebras.group_algebra import GroupAlgebraFunctor
            sage: GroupAlgebraFunctor(CyclicPermutationGroup(17)).group() == CyclicPermutationGroup(17)
            True
         """
        return self.__group

    def _apply_functor(self, base_ring):
        r"""
        Create the group algebra with given base ring over ``self.group()``.

        INPUT :

        - ``base_ring`` - the base ring of the group algebra.

        OUTPUT:

        A group algebra.

        EXAMPLES::

            sage: from sage.algebras.group_algebra import GroupAlgebraFunctor
            sage: F = GroupAlgebraFunctor(CyclicPermutationGroup(17))
            sage: F(QQ)
            Group algebra of Cyclic group of order 17 as a permutation group
             over Rational Field
        """
        return self.__group.algebra(base_ring)

    def _apply_functor_to_morphism(self, f):
        r"""
        Lift a homomorphism of rings to the corresponding homomorphism
        of the group algebras of ``self.group()``.

        INPUT:

        - ``f`` -- a morphism of rings

        OUTPUT:

        A morphism of group algebras.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: A = GroupAlgebra(G, ZZ)
            sage: h = sage.categories.morphism.SetMorphism(Hom(ZZ, GF(5), Rings()), lambda x: GF(5)(x))
            sage: hh = A.construction()[0](h)
            sage: hh(A.0 + 5 * A.1)
            (1,2,3)
        """
        codomain = self(f.codomain())
        return SetMorphism(Hom(self(f.domain()), codomain, Rings()),
                           lambda x: sum(codomain(g) * f(c) for (g, c) in six.iteritems(dict(x))))

def group_algebra(G, R=None):
    """
    Construct a group algebra.

    INPUT:

    - ``group`` -- a group
    - ``base_ring`` -- (default: `\ZZ`) a commutative ring
    """
    if R is None:
        R = IntegerRing()
    if R not in Rings():
        G,R = R,G
    return G.algebra(R)

class GroupAlgebra(CombinatorialFreeModule):
    r"""
    Create the given group algebra.

    EXAMPLES::

        sage: GroupAlgebra(GL(3, GF(7)))
        Group algebra of General Linear Group of degree 3 over Finite Field of size 7
         over Integer Ring
        sage: GroupAlgebra(GL(3, GF(7)), QQ)
        Group algebra of General Linear Group of degree 3 over Finite Field of size 7
         over Rational Field
        sage: GroupAlgebra(1)
        Traceback (most recent call last):
        ...
        TypeError: "1" is not a group

        sage: GroupAlgebra(SU(2, GF(4, 'a')), IntegerModRing(12)).category()
        Category of finite group algebras over Ring of integers modulo 12
        sage: GroupAlgebra(KleinFourGroup()) is GroupAlgebra(KleinFourGroup())
        True

    The one of the group indexes the one of this algebra::

        sage: A = GroupAlgebra(DihedralGroup(6), QQ)
        sage: A.one_basis()
        ()
        sage: A.one()
        ()

    The product of two basis elements is induced by the product of the
    corresponding elements of the group::

        sage: A = GroupAlgebra(DihedralGroup(3), QQ)
        sage: (a, b) = A.group().gens()
        sage: a*b
        (1,2)
        sage: A.product_on_basis(a, b)
        (1,2)

    The basis elements are group-like for the coproduct:
    `\Delta(g) = g \otimes g`::

        sage: A = GroupAlgebra(DihedralGroup(3), QQ)
        sage: (a, b) = A.group().gens()
        sage: A.coproduct_on_basis(a)
        (1,2,3) # (1,2,3)

    The counit on the basis elements is 1::

        sage: A = GroupAlgebra(DihedralGroup(6), QQ)
        sage: (a, b) = A.group().gens()
        sage: A.counit_on_basis(a)
        1

    The antipode on basis elements is given by `\chi(g) = g^{-1}`::

        sage: A = GroupAlgebra(DihedralGroup(3), QQ)
        sage: (a, b) = A.group().gens(); a
        (1,2,3)
        sage: A.antipode_on_basis(a)
        (1,3,2)

    TESTS::

        sage: A = GroupAlgebra(GL(3, GF(7)))
        sage: A.has_coerce_map_from(GL(3, GF(7)))
        True
        sage: G = SymmetricGroup(5)
        sage: x,y = G.gens()
        sage: A = GroupAlgebra(G)
        sage: A( A(x) )
        (1,2,3,4,5)

    These tests are to be moved!!!

    EXAMPLES::

        sage: GroupAlgebra(GL(3, GF(7)))
        Group algebra of General Linear Group of degree 3 over Finite Field of size 7
         over Integer Ring

    TESTS::

        sage: GroupAlgebra(AbelianGroup(1)) == GroupAlgebra(AbelianGroup(1))
        True
        sage: GroupAlgebra(AbelianGroup(1), QQ) == GroupAlgebra(AbelianGroup(1), ZZ)
        False
        sage: GroupAlgebra(AbelianGroup(2)) == GroupAlgebra(AbelianGroup(1))
        False
        sage: A = GroupAlgebra(KleinFourGroup(), ZZ)
        sage: B = GroupAlgebra(KleinFourGroup(), QQ)
        sage: A == B
        False
        sage: A == A
        True

        sage: GroupAlgebra(SymmetricGroup(2)).is_commutative()
        True
        sage: GroupAlgebra(SymmetricGroup(3)).is_commutative()
        False
    """
    # coercion methods:

    def _coerce_map_from_(self, S):
        r"""
        True if there is a coercion from ``S`` to ``self``, False otherwise.
        The actual coercion is done by the :meth:`_element_constructor_`
        method.

        INPUT:

        -  ``S`` - a Sage object

        The objects that coerce into a group algebra `k[G]` are:

        - any group algebra `R[H]` as long as `R` coerces into `k` and
          `H` coerces into `G`.

        - any ring `R` which coerces into `k`

        - any group `H` which coerces into either `k` or `G`.

        Note that if `H` is a group which coerces into both `k` and
        `G`, then Sage will always use the map to `k`.  For example,
        if `\ZZ` is the ring (or group) of integers, then `\ZZ` will
        coerce to any `k[G]`, by sending `\ZZ` to `k`.

        EXAMPLES::

            sage: A = GroupAlgebra(SymmetricGroup(4), QQ)
            sage: B = GroupAlgebra(SymmetricGroup(3), ZZ)
            sage: A.has_coerce_map_from(B)
            True
            sage: B.has_coerce_map_from(A)
            False
            sage: A.has_coerce_map_from(ZZ)
            True
            sage: A.has_coerce_map_from(CC)
            False
            sage: A.has_coerce_map_from(SymmetricGroup(5))
            False
            sage: A.has_coerce_map_from(SymmetricGroup(2))
            True
        """
        from sage.rings.ring import is_Ring
        from sage.groups.old import Group
        k = self.base_ring()
        G = self.group()
        if isinstance(S, GroupAlgebra):
            return (k.has_coerce_map_from(S.base_ring())
                    and G.has_coerce_map_from(S.group()))
        if is_Ring(S):
            return k.has_coerce_map_from(S)
        if isinstance(S, Group):
            return k.has_coerce_map_from(S) or G.has_coerce_map_from(S)

    def _element_constructor_(self, x):
        r"""
        Try to turn ``x`` into an element of ``self``.

        INPUT:

        - ``x`` -- an element of some group algebra or of a
          ring or of a group

        OUTPUT:

        ``x`` as a member of ``self``.

        EXAMPLES::

            sage: G = KleinFourGroup()
            sage: f = G.gen(0)
            sage: ZG = GroupAlgebra(G)
            sage: ZG(f)  # indirect doctest
            (3,4)
            sage: ZG(1) == ZG(G(1))
            True
            sage: G = AbelianGroup(1)
            sage: ZG = GroupAlgebra(G)
            sage: f = ZG.group().gen()
            sage: ZG(FormalSum([(1,f), (2, f**2)]))
            f + 2*f^2
            sage: G = GL(2,7)
            sage: OG = GroupAlgebra(G, ZZ[sqrt(5)])
            sage: OG(2)
            2*[1 0]
            [0 1]
            sage: OG(G(2)) # conversion is not the obvious one
            [2 0]
            [0 2]
            sage: OG(FormalSum([ (1, G(2)), (2, RR(0.77)) ]) )
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer
            sage: OG(OG.base_ring().basis()[1])
            sqrt5*[1 0]
            [0 1]
        """
        k = self.base_ring()

        #Coerce ints to Integers
        if isinstance(x, int):
            x = Integer(x)

        if x in k:
            if x == 0:
                return self.zero()
            else:
                return k(x) * self.one()

        G = self.group()
        S = x.parent()
        if S is G:
            return self.monomial(x)

        from sage.rings.ring import is_Ring
        from sage.structure.formal_sum import FormalSum

        if isinstance(S, GroupAlgebra):
            if self.has_coerce_map_from(S):
                # coerce monomials, coerce coefficients, reassemble
                d = x.monomial_coefficients()
                new_d = {}
                for g in d:
                    g1 = G(g)
                    if g1 in new_d:
                        new_d[g1] += k(d[g]) + new_d[g1]
                    else:
                        new_d[g1] = k(d[g])
                return self._from_dict(new_d)

        elif is_Ring(S):
            # coerce to multiple of identity element
            return k(x) * self.one()

        elif isinstance(x, FormalSum) and k.has_coerce_map_from(S.base_ring()):
            y = [(G(g), k(coeff)) for coeff,g in x]
            return self.sum_of_terms(y)

        # Check whether group coerces to base_ring first.
        if k.has_coerce_map_from(S):
            return k(x) * self.one()
        if G.has_coerce_map_from(S):
            return self.monomial(G(x))

        raise TypeError("Don't know how to create an element of %s from %s" % \
                             (self, x))

