r"""
Group Algebras

This module implements the category of group algebras for arbitrary
groups over arbitrary commutative rings. This also provides much of
the functionality for the implementation of a group algebra using
:class:`CombinatorialFreeModule`.

EXAMPLES::

    sage: D4 = DihedralGroup(4)
    sage: kD4 = GroupAlgebra(D4, GF(7))
    sage: a = kD4.an_element(); a
    () + 4*(1,2,3,4) + 2*(1,4)(2,3)
    sage: a * a
    5*() + (2,4) + (1,2,3,4) + (1,3) + 2*(1,3)(2,4) + 4*(1,4)(2,3)

Given the group and the base ring, the corresponding group
algebra is unique::

    sage: A = GroupAlgebra(GL(3, QQ), ZZ)
    sage: B = GroupAlgebra(GL(3, QQ), ZZ)
    sage: A is B
    True
    sage: C = GroupAlgebra(GL(3, QQ), QQ)
    sage: A == C
    False

As long as there is no natural map from the group to the base ring,
you can easily convert elements of the group to the group algebra::

    sage: A = GroupAlgebra(DihedralGroup(2), ZZ)
    sage: g = DihedralGroup(2).gen(0); g
    (3,4)
    sage: A(g)
    (3,4)
    sage: A(2) * g
    2*(3,4)

Since there is a natural inclusion from the dihedral group `D_2` of
order 4 into the symmetric group `S_4` of order 4!, and since there is
a natural map from the integers to the rationals, there is a natural
map from `\ZZ[D_2]` to `\QQ[S_4]`::

    sage: A = GroupAlgebra(DihedralGroup(2), ZZ)
    sage: B = GroupAlgebra(SymmetricGroup(4), QQ)
    sage: a = A.an_element(); a
    () + 3*(3,4) + 3*(1,2)
    sage: b = B.an_element(); b
    () + 2*(1,2) + 4*(1,2,3,4)
    sage: B(a)
    () + 3*(3,4) + 3*(1,2)
    sage: a * b  # a is automatically converted to an element of B
    7*() + 3*(3,4) + 5*(1,2) + 6*(1,2)(3,4) + 12*(1,2,3) + 4*(1,2,3,4) + 12*(1,3,4)
    sage: parent(a * b)
    Group algebra of Symmetric group of order 4! as a permutation group
     over Rational Field

    sage: G = GL(3, GF(7))
    sage: ZG = GroupAlgebra(G)
    sage: c, d = G.random_element(), G.random_element()
    sage: zc, zd = ZG(c), ZG(d)
    sage: zc * d == zc * zd  # d is automatically converted to an element of ZG
    True

There is no obvious map in the other direction, though::

    sage: A(b)
    Traceback (most recent call last):
    ...
    TypeError: Don't know how to create an element of
     Group algebra of Dihedral group of order 4 as a permutation group over Integer Ring
     from () + 2*(1,2) + 4*(1,2,3,4)

Group algebras have the structure of Hopf algebras::

    sage: a = kD4.an_element(); a
    () + 4*(1,2,3,4) + 2*(1,4)(2,3)
    sage: a.antipode()
    () + 4*(1,4,3,2) + 2*(1,4)(2,3)
    sage: a.coproduct()
    () # () + 4*(1,2,3,4) # (1,2,3,4) + 2*(1,4)(2,3) # (1,4)(2,3)

.. NOTE::

    As alluded to above, it is problematic to make group algebras fit
    nicely with Sage's coercion model. The problem is that (for
    example) if G is the additive group `(\ZZ,+)`, and `R = \ZZ[G]` is
    its group ring, then the integer 2 can be coerced into R in two
    ways -- via G, or via the base ring -- and *the answers are
    different*. In practice we get around this by preventing elements
    of a group `H` from coercing automatically into a group ring
    `k[G]` if `H` coerces into both `k` and `G`.  This is unfortunate,
    but it seems like the most sensible solution in this ambiguous
    situation.

AUTHOR:

- David Loeffler (2008-08-24): initial version
- Martin Raum (2009-08): update to use new coercion model -- see
  :trac:`6670`.
- John Palmieri (2011-07): more updates to coercion, categories, etc.,
  group algebras constructed using CombinatorialFreeModule -- see
  :trac:`6670`.
"""

#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.algebra_functor import AlgebrasCategory

class GroupAlgebras(AlgebrasCategory):
    r"""
    The category of group algebras over a given base ring.

    EXAMPLES::

        sage: C = GroupAlgebras(IntegerRing()); C
        Category of group algebras over Integer Ring
        sage: C.super_categories()
        [Category of hopf algebras with basis over Integer Ring,
         Category of monoid algebras over Integer Ring]

    We can also construct group algebras by::

        sage: C is Groups().Algebras(ZZ)
        True

    Here is how to create the group algebra of a group `G`::

        sage: G = DihedralGroup(5)
        sage: QG = G.algebra(QQ); QG
        Group algebra of Dihedral group of order 10 as a permutation group over Rational Field

    and an example of computation::

        sage: g = G.an_element(); g
        (1,2,3,4,5)
        sage: (QG.term(g) + 1)**3
        () + 3*(1,2,3,4,5) + 3*(1,3,5,2,4) + (1,4,2,5,3)

    .. TODO::

        - Check which methods would be better located in
          ``Monoid.Algebras`` or ``Groups.Finite.Algebras``.

    TESTS::

        sage: A = GroupAlgebras(QQ).example(GL(3, GF(11)))
        sage: A.one_basis()
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: A = SymmetricGroupAlgebra(QQ,4)
        sage: x = Permutation([4,3,2,1])
        sage: A.product_on_basis(x,x)
        [1, 2, 3, 4]

        sage: C = GroupAlgebras(ZZ)
        sage: TestSuite(C).run()
    """
    def extra_super_categories(self):
        """
        Implement the fact that the algebra of a group is a Hopf
        algebra.

        EXAMPLES::

            sage: C = Groups().Algebras(QQ)
            sage: C.extra_super_categories()
            [Category of hopf algebras over Rational Field]
            sage: sorted(C.super_categories(), key=str)
            [Category of hopf algebras with basis over Rational Field,
             Category of monoid algebras over Rational Field]
        """
        from sage.categories.hopf_algebras import HopfAlgebras
        return [HopfAlgebras(self.base_ring())]

    def example(self, G=None):
        """
        Return an example of group algebra.

        EXAMPLES::

            sage: GroupAlgebras(QQ['x']).example()
            Group algebra of Dihedral group of order 8 as a permutation group over Univariate Polynomial Ring in x over Rational Field

        An other group can be specified as optional argument::

            sage: GroupAlgebras(QQ).example(AlternatingGroup(4))
            Group algebra of Alternating group of order 4!/2 as a permutation group over Rational Field
        """
        from sage.groups.perm_gps.permgroup_named import DihedralGroup
        if G is None:
            G = DihedralGroup(4)
        return G.algebra(self.base_ring())

    class ParentMethods:
        def __init_extra__(self):
            """
            Enable coercion from the defining group.

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
            if not self.base_ring().has_coerce_map_from(self.group()):
                ## some matrix groups assume that coercion is only valid to
                ## other matrix groups. This is a workaround
                ## call _element_constructor_ to coerce group elements
                #try :
                self._populate_coercion_lists_(coerce_list=[self.group()])

        def _repr_(self):
            r"""
            Return the string representation of `self`.

            EXAMPLES::

                sage: A = Groups().example().algebra(QQ); A
                Group algebra of General Linear Group of degree 4 over Rational Field
                 over Rational Field
                sage: A._name = "foo"
                sage: A
                foo over Rational Field
                sage: A = GroupAlgebra(KleinFourGroup(), ZZ)
                sage: A
                Group algebra of The Klein 4 group of order 4, as a permutation group
                 over Integer Ring
            """
            if hasattr(self, "_name"):
                return self._name + " over {}".format(self.base_ring())
            else:
                return 'Group algebra of {} over {}'.format(self.basis().keys(),
                                                            self.base_ring())

        def _latex_(self):
            r"""
            Latex string of ``self``.

            EXAMPLES::

                sage: A = GroupAlgebra(KleinFourGroup(), ZZ)
                sage: latex(A) # indirect doctest
                \Bold{Z}[\langle (3,4), (1,2) \rangle]
            """
            from sage.misc.all import latex
            return "%s[%s]" % (latex(self.base_ring()), latex(self.group()))

        def group(self):
            r"""
            Return the underlying group of the group algebra.

            EXAMPLES::

                sage: GroupAlgebras(QQ).example(GL(3, GF(11))).group()
                General Linear Group of degree 3 over Finite Field of size 11
                sage: SymmetricGroup(10).algebra(QQ).group()
                Symmetric group of order 10! as a permutation group
            """
            return self.basis().keys()

        @cached_method
        def algebra_generators(self):
            r"""
            Return generators of this group algebra (as an algebra).

            EXAMPLES::

                sage: GroupAlgebras(QQ).example(AlternatingGroup(10)).algebra_generators()
                Finite family {(8,9,10): (8,9,10), (1,2,3,4,5,6,7,8,9): (1,2,3,4,5,6,7,8,9)}

                sage: A = GroupAlgebra(DihedralGroup(3), QQ); A
                Group algebra of Dihedral group of order 6 as a permutation group
                 over Rational Field
                sage: A.algebra_generators()
                Finite family {(1,3): (1,3), (1,2,3): (1,2,3)}
            """
            from sage.sets.family import Family
            return Family(self.group().gens(), self.monomial)

        def ngens(self):
            r"""
            Return the number of generators of ``self``.

            EXAMPLES::

                sage: GroupAlgebra(SL2Z).ngens()
                2
                sage: GroupAlgebra(DihedralGroup(4), RR).ngens()
                2
            """
            return self.algebra_generators().cardinality()

        def gen(self, i=0):
            r"""
            Return the ``i``-th generator of ``self``.

            EXAMPLES::

                sage: A = GroupAlgebra(GL(3, GF(7)))
                sage: A.gen(0)
                [3 0 0]
                [0 1 0]
                [0 0 1]
            """
            return self.monomial(self.group().gen(i))

        def construction(self):
            r"""
            Return the functorial construction of ``self``.

            EXAMPLES::

                sage: A = GroupAlgebra(KleinFourGroup(), QQ)
                sage: A.construction()
                (GroupAlgebraFunctor, Rational Field)
            """
            from sage.algebras.group_algebra import GroupAlgebraFunctor
            return GroupAlgebraFunctor(self.group()), self.base_ring()

        @cached_method
        def center_basis(self):
            r"""
            Return a basis of the center of the group algebra.

            The canonical basis of the center of the group algebra
            is the family `(f_\sigma)_{\sigma\in C}`, where `C` is
            any collection of representatives of the conjugacy
            classes of the group, and `f_\sigma` is the sum of the
            elements in the conjugacy class of `\sigma`.

            OUTPUT:

            - ``tuple`` of elements of ``self``

            .. WARNING::

                - This method requires the underlying group to
                  have a method ``conjugacy_classes``
                  (every permutation group has one, thanks GAP!).

            EXAMPLES::

                sage: SymmetricGroup(3).algebra(QQ).center_basis()
                ((), (2,3) + (1,2) + (1,3), (1,2,3) + (1,3,2))

            .. SEEALSO::

                - :meth:`Groups.Algebras.ElementMethods.central_form`
                - :meth:`Monoids.Algebras.ElementMethods.is_central`
            """
            return tuple([self.sum_of_monomials(conj) for conj  in
                          self.basis().keys().conjugacy_classes()])

        # Coalgebra structure

        def coproduct_on_basis(self, g):
            r"""
            Return the coproduct of the element ``g`` of the basis.

            Each basis element ``g`` is group-like. This method is
            used to compute the coproduct of any element.

            EXAMPLES::

                sage: A = CyclicPermutationGroup(6).algebra(ZZ); A
                Group algebra of Cyclic group of order 6 as a permutation group over Integer Ring
                sage: g = CyclicPermutationGroup(6).an_element(); g
                (1,2,3,4,5,6)
                sage: A.coproduct_on_basis(g)
                (1,2,3,4,5,6) # (1,2,3,4,5,6)
                sage: a = A.an_element(); a
                () + 3*(1,2,3,4,5,6) + 3*(1,3,5)(2,4,6)
                sage: a.coproduct()
                () # () + 3*(1,2,3,4,5,6) # (1,2,3,4,5,6) + 3*(1,3,5)(2,4,6) # (1,3,5)(2,4,6)
            """
            from sage.categories.tensor import tensor
            g = self.term(g)
            return tensor([g, g])

        def antipode_on_basis(self,g):
            r"""
            Return the antipode of the element ``g`` of the basis.

            Each basis element ``g`` is group-like, and so has
            antipode `g^{-1}`. This method is used to compute the
            antipode of any element.

            EXAMPLES::

                sage: A = CyclicPermutationGroup(6).algebra(ZZ); A
                Group algebra of Cyclic group of order 6 as a permutation group over Integer Ring
                sage: g = CyclicPermutationGroup(6).an_element();g
                (1,2,3,4,5,6)
                sage: A.antipode_on_basis(g)
                (1,6,5,4,3,2)
                sage: a = A.an_element(); a
                () + 3*(1,2,3,4,5,6) + 3*(1,3,5)(2,4,6)
                sage: a.antipode()
                () + 3*(1,5,3)(2,6,4) + 3*(1,6,5,4,3,2)
            """
            return self.term(~g)

        def counit_on_basis(self,g):
            r"""
            Return the counit of the element ``g`` of the basis.

            Each basis element ``g`` is group-like, and so has
            counit `1`. This method is used to compute the
            counit of any element.

            EXAMPLES::

                sage: A=CyclicPermutationGroup(6).algebra(ZZ);A
                Group algebra of Cyclic group of order 6 as a permutation group over Integer Ring
                sage: g=CyclicPermutationGroup(6).an_element();g
                (1,2,3,4,5,6)
                sage: A.counit_on_basis(g)
                1
            """
            return self.base_ring().one()

        def counit(self,x):
            r"""
            Return the counit of the element ``x`` of the group
            algebra.

            This is the sum of all coefficients of ``x`` with respect
            to the standard basis of the group algebra.

            EXAMPLES::

                sage: A = CyclicPermutationGroup(6).algebra(ZZ); A
                Group algebra of Cyclic group of order 6 as a permutation group over Integer Ring
                sage: a = A.an_element(); a
                () + 3*(1,2,3,4,5,6) + 3*(1,3,5)(2,4,6)
                sage: a.counit()
                7
            """
            return self.base_ring().sum(x.coefficients())

        def is_field(self, proof=True):
            r"""
            Return ``True`` if ``self`` is a field.

            This is always false unless ``self.group()`` is trivial
            and ``self.base_ring()`` is a field.

            EXAMPLES::

                sage: GroupAlgebra(SymmetricGroup(2)).is_field()
                False
                sage: GroupAlgebra(SymmetricGroup(1)).is_field()
                False
                sage: GroupAlgebra(SymmetricGroup(1), QQ).is_field()
                True
            """
            if not self.base_ring().is_field(proof):
                return False
            return (self.group().order() == 1)

        def is_finite(self):
            r"""
            Return ``True`` if ``self`` is finite, which is true if and only
            if ``self.group()`` and ``self.base_ring()`` are both finite.

            EXAMPLES::

                sage: GroupAlgebra(SymmetricGroup(2), IntegerModRing(10)).is_finite()
                True
                sage: GroupAlgebra(SymmetricGroup(2)).is_finite()
                False
                sage: GroupAlgebra(AbelianGroup(1), IntegerModRing(10)).is_finite()
                False
            """
            return (self.base_ring().is_finite() and self.group().is_finite())

        def is_integral_domain(self, proof=True):
            r"""
            Return ``True`` if ``self`` is an integral domain.

            This is false unless ``self.base_ring()`` is an integral
            domain, and even then it is false unless ``self.group()``
            has no nontrivial elements of finite order. I don't know
            if this condition suffices, but it obviously does if the
            group is abelian and finitely generated.

            EXAMPLES::

                sage: GroupAlgebra(SymmetricGroup(2)).is_integral_domain()
                False
                sage: GroupAlgebra(SymmetricGroup(1)).is_integral_domain()
                True
                sage: GroupAlgebra(SymmetricGroup(1), IntegerModRing(4)).is_integral_domain()
                False
                sage: GroupAlgebra(AbelianGroup(1)).is_integral_domain()
                True
                sage: GroupAlgebra(AbelianGroup(2, [0,2])).is_integral_domain()
                False
                sage: GroupAlgebra(GL(2, ZZ)).is_integral_domain()
                False
            """
            from sage.sets.set import Set
            ans = False
            try:
                if self.base_ring().is_integral_domain():
                    if self.group().is_finite():
                        if self.group().order() > 1:
                            ans = False
                        else:
                            ans = True
                    else:
                        if self.group().is_abelian():
                            invs = self.group().invariants()
                            if Set(invs) != Set([0]):
                                ans = False
                            else:
                                ans = True
                        else:
                            raise NotImplementedError
                else:
                    ans = False
            except (AttributeError, NotImplementedError):
                if proof:
                    raise NotImplementedError("cannot determine whether self is an integral domain")

            return ans

        # I haven't written is_noetherian(), because I don't know when group
        # algebras are noetherian, and I haven't written is_prime_field(), because
        # I don't know if that means "is canonically isomorphic to a prime field"
        # or "is identical to a prime field".

        def random_element(self, n=2):
            r"""
            Return a 'random' element of ``self``.

            INPUT:

            - ``n`` -- integer (default: 2); number of summands

            ALGORITHM:

            Return a sum of ``n`` terms, each of which is formed by
            multiplying a random element of the base ring by a random
            element of the group.

            EXAMPLES::

                sage: GroupAlgebra(DihedralGroup(6), QQ).random_element()
                -1/95*() - 1/2*(1,4)(2,5)(3,6)
                sage: GroupAlgebra(SU(2, 13), QQ).random_element(1)
                1/2*[       0 4*a + 11]
                [2*a + 12        4]
            """
            a = self(0)
            for i in range(n):
                a += self.term(self.group().random_element(),
                               self.base_ring().random_element())
            return a

    class ElementMethods:

        def central_form(self):
            r"""
            Return ``self`` expressed in the canonical basis of the center
            of the group algebra.

            INPUT:

            - ``self`` -- an element of the center of the group algebra

            OUTPUT:

            - A formal linear combination of the conjugacy class
              representatives representing its coordinates in the
              canonical basis of the center. See
              :meth:`Groups.Algebras.ParentMethods.center_basis` for
              details.

            .. WARNING::

                - This method requires the underlying group to
                  have a method ``conjugacy_classes_representatives``
                  (every permutation group has one, thanks GAP!).
                - This method does not check that the element is
                  indeed central. Use the method
                  :meth:`Monoids.Algebras.ElementMethods.is_central`
                  for this purpose.
                - This function has a complexity linear in the
                  number of conjugacy classes of the group. One
                  could easily implement a function whose
                  complexity is linear in the size of the support
                  of ``self``.

            EXAMPLES::

                sage: QS3 = SymmetricGroup(3).algebra(QQ)
                sage: A = QS3([2,3,1]) + QS3([3,1,2])
                sage: A.central_form()
                B[(1,2,3)]
                sage: QS4 = SymmetricGroup(4).algebra(QQ)
                sage: B = sum(len(s.cycle_type())*QS4(s) for s in Permutations(4))
                sage: B.central_form()
                4*B[()] + 3*B[(1,2)] + 2*B[(1,2)(3,4)] + 2*B[(1,2,3)] + B[(1,2,3,4)]

                sage: QG = GroupAlgebras(QQ).example(PermutationGroup([[(1,2,3),(4,5)],[(3,4)]]))
                sage: sum(i for i in QG.basis()).central_form()
                B[()] + B[(4,5)] + B[(3,4,5)] + B[(2,3)(4,5)] + B[(2,3,4,5)] + B[(1,2)(3,4,5)] + B[(1,2,3,4,5)]

            .. SEEALSO::

                - :meth:`Groups.Algebras.ParentMethods.center_basis`
                - :meth:`Monoids.Algebras.ElementMethods.is_central`
            """
            from sage.combinat.free_module import CombinatorialFreeModule
            conj_classes_reps = self.parent().basis().keys().conjugacy_classes_representatives()
            Z = CombinatorialFreeModule(self.base_ring(), conj_classes_reps)
            return sum(self[i] * Z.basis()[i] for i in Z.basis().keys())

