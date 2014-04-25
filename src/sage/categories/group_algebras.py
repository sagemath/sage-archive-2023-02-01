r"""
GroupAlgebras
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.sets.family import Family

class GroupAlgebras(Category_over_base_ring):
    """
    The category of group algebras over the base ring.

    EXAMPLES::

        sage: GroupAlgebras(IntegerRing())
        Category of group algebras over Integer Ring
        sage: GroupAlgebras(IntegerRing()).super_categories()
        [Category of monoid algebras over Integer Ring]

    Here is how to create the group algebra of a group G (this has to
    be changed if there are no bugs in this category)::

        sage: G=DihedralGroup(5)
        sage: QG=GroupAlgebras(QQ).example(G); QG
        The group algebra of the Dihedral group of order 10 as a permutation group over Rational Field

    and an example of computation::

        sage: g=DihedralGroup(5).an_element(); g
        (1,2,3,4,5)
        sage: (QG.term(g) + 1)**3
        B[()] + 3*B[(1,2,3,4,5)] + 3*B[(1,3,5,2,4)] + B[(1,4,2,5,3)]

    .. TODO::

        - Some methods can be transferred to MonoidAlgebras and some to a new category FiniteDimensionalGroupAlgebras
        - The Hopf algebra structure has not yet been implemented

    TESTS::

        sage: C = GroupAlgebras(ZZ)
        sage: TestSuite(C).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GroupAlgebras(QQ).super_categories()
            [Category of monoid algebras over Rational Field]
        """
        from monoid_algebras import MonoidAlgebras
        return [MonoidAlgebras(self.base_ring())] # TODO: could become Kac algebras / Hopf algebras

    def example(self, G = None):
        """
        Returns an example of group algebra

        EXAMPLES::

            sage: GroupAlgebras(QQ[x]).example()
            The group algebra of the Dihedral group of order 8 as a permutation group over Univariate Polynomial Ring in x over Rational Field

        An other group can be specified as optional argument::

            sage: GroupAlgebras(QQ).example(SymmetricGroup(4))
            The group algebra of the Symmetric group of order 4! as a permutation group over Rational Field

        """
        from sage.categories.examples.group_algebras import MyGroupAlgebra
        from sage.groups.perm_gps.permgroup_named import DihedralGroup
        if G is None:
            G = DihedralGroup(4)
        return MyGroupAlgebra(self.base_ring(), G)

    class ParentMethods:

        @abstract_method
        def group(self):
            r"""
            Returns the underlying group of the group algebra

            EXAMPLES::

                sage: GroupAlgebras(QQ).example(GL(3, GF(11))).group()
                General Linear Group of degree 3 over Finite Field of size 11
                sage: SymmetricGroupAlgebra(QQ,10).group()
                Symmetric group of order 10! as a permutation group
            """

        @cached_method
        def one_basis(self):
            """
            Returns the unit

            EXAMPLES::

                sage: A = GroupAlgebras(QQ).example(GL(3, GF(11)))
                sage: A.one_basis()
                [1 0 0]
                [0 1 0]
                [0 0 1]
            """
            return self.group().one()

        def product_on_basis(self, g1, g2):
            """
            Returns the product of two elements

            EXAMPLES::

                sage: A = SymmetricGroupAlgebra(QQ,4)
                sage: x = Permutation([4,3,2,1])
                sage: A.product_on_basis(x,x)
                [1, 2, 3, 4]
            """
            return self.basis()[g1 * g2]

    #    def coproduct_on_basis(self, g):
    #        g = self.term(g)
    #        return tensor([g, g])
    # The Hopf algebra structure of Hopf algebra is not yet implemented

        def algebra_generators(self):
            r"""
            Returns generators of this group algebra (as an algebra).

            EXAMPLES::

                sage: GroupAlgebras(QQ).example(SymmetricGroup(10)).algebra_generators()
                Finite family {(1,2): B[(1,2)], (1,2,3,4,5,6,7,8,9,10): B[(1,2,3,4,5,6,7,8,9,10)]}

            .. NOTE::

                This function is overloaded for SymmetricGroupAlgebras to return Permutations and not Elements of the symmetric group
            """
            return Family(self.group().gens(),self.term)

        def _conjugacy_classes_representatives_underlying_group(self):
            r"""
            Returns a complete list of representatives of conjugacy classes

            This works only for permutations group. The ordering is that given by GAP.

            EXAMPLES::

                sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3,4)]])
                sage: SG = GroupAlgebras(QQ).example(G)
                sage: SG._conjugacy_classes_representatives_underlying_group()
                [(), (2,4), (1,2)(3,4), (1,2,3,4), (1,3)(2,4)]

            .. NOTE::

                This function is overloaded for SymmetricGroupAlgebras to return Permutations and not Elements of the symmetric group::

                sage: SymmetricGroupAlgebra(ZZ,3)._conjugacy_classes_representatives_underlying_group()
                [[2, 3, 1], [2, 1, 3], [1, 2, 3]]
            """
            return self.group().conjugacy_classes_representatives()

        def center(self):
            r"""
            Returns the center of the group algebra.

            The element of the canonical basis (sum of the elements of the group in the same conjugacy classes) are indexed by one element of the class.

            EXAMPLES::

                sage: SymmetricGroupAlgebra(ZZ,3).center()
                Free module generated by {[2, 3, 1], [2, 1, 3], [1, 2, 3]} over Integer Ring

            .. NOTE::

                The product has not been implemented yet.
            """
            from sage.categories.all import Algebras
            from sage.combinat.free_module import CombinatorialFreeModule
            return CombinatorialFreeModule(self.base_ring(),self._conjugacy_classes_representatives_underlying_group())

    class ElementMethods:
        def is_central(self):
            r"""
            Returns True if the element is central and False otherwise.

            EXAMPLES::

                sage: SG4=SymmetricGroupAlgebra(ZZ,4)
                sage: SG4(1).is_central()
                True
                sage: SG4(Permutation([1,3,2,4])).is_central()
                False
                sage: A=GroupAlgebras(QQ).example(); A
                The group algebra of the Dihedral group of order 8 as a permutation group over Rational Field
                sage: sum(i for i in A.basis()).is_central()
                True
            """
            return all([i*self == self*i for i in self.parent().algebra_generators()])

        def central_form(self):
            r"""
            Returns ``self`` as an element of a center of the group algebra in the canonical basis of the center of the group algebra (family of sums of elements in a conjugacy class).

            Note that the element of the canonical basis of the center (sum of the elements of the group in the same conjugacy classes) are indexed by one element of the class.

            WARNINGS::

                - This function needs the underlying group to have a method conjugacy_classes_representatives (every permutation group has one, thanks GAP!).
                - It does not check that the element is indeed central. Use the method :meth:`.is_central` for this purpose.

            EXAMPLES::

                sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
                sage: A=QS3([2,3,1])+QS3([3,1,2])
                sage: A.central_form()
                B[[2, 3, 1]]
                sage: QS4 = SymmetricGroupAlgebra(QQ, 4)
                sage: B=sum(len(s.cycle_type())*QS4(s) for s in Permutations(4))
                sage: B.central_form()
                4*B[[1, 2, 3, 4]] + 3*B[[2, 1, 3, 4]] + 2*B[[2, 1, 4, 3]] + 2*B[[2, 3, 1, 4]] + B[[2, 3, 4, 1]]
                sage: QG=GroupAlgebras(QQ).example(PermutationGroup([[(1,2,3),(4,5)],[(3,4)]]))
                sage: sum(i for i in QG.basis()).central_form()
                B[()] + B[(4,5)] + B[(3,4,5)] + B[(2,3)(4,5)] + B[(2,3,4,5)] + B[(1,2)(3,4,5)] + B[(1,2,3,4,5)]

            .. NOTE::

                This function has a complexity linear in the number of conjugacy classes of the group. One could have done easily a function, whose complexity is linear in the size of the support of ``self``.
            """
            Z=self.parent().center()
            return sum(self[i] * Z.basis()[i] for i in Z.basis().keys())

