r"""
Groups
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.monoids import Monoids
from sage.categories.algebra_functor import AlgebrasCategory

class Groups(CategoryWithAxiom):
    """
    The category of (multiplicative) groups, i.e. monoids with
    inverses.

    EXAMPLES::

        sage: Groups()
        Category of groups
        sage: Groups().super_categories()
        [Category of monoids, Category of inverse unital magmas]

    TESTS::

        sage: TestSuite(Groups()).run()
    """
    _base_category_class_and_axiom = (Monoids, "Inverse")

    def example(self):
        """
        EXAMPLES::

            sage: Groups().example()
            General Linear Group of degree 4 over Rational Field
        """
        from sage.rings.rational_field import QQ
        from sage.groups.matrix_gps.linear import GL
        return GL(4,QQ)

    class ParentMethods:

        def group_generators(self):
            """
            Returns group generators for self.

            This default implementation calls :meth:`.gens`, for
            backward compatibility.

            EXAMPLES::

                sage: A = AlternatingGroup(4)
                sage: A.group_generators()
                Family ((2,3,4), (1,2,3))
            """
            from sage.sets.family import Family
            return Family(self.gens())

        def _test_inverse(self, **options):
            """
            Run generic tests on the method :meth:`.__invert__`.

            See also: :class:`TestSuite`.

            EXAMPLES::

                sage: G = SymmetricGroup(3)
                sage: G._test_inverse()
            """
            tester = self._tester(**options)
            for x in tester.some_elements():
                tester.assertEquals(x * ~x, self.one())
                tester.assertEquals(~x * x, self.one())

        def semidirect_product(self, N, mapping, check = True):
            r"""
            The semi-direct product of two groups

            EXAMPLES::

                sage: G = Groups().example()
                sage: G.semidirect_product(G,Morphism(G,G))
                Traceback (most recent call last):
                ...
                NotImplementedError: semidirect product of General Linear Group of degree 4 over Rational Field and General Linear Group of degree 4 over Rational Field not yet implemented
            """
            raise NotImplementedError("semidirect product of %s and %s not yet implemented"%(self, N))

        def holomorph(self):
            r"""
            The holomorph of a group

            The holomorph of a group `G` is the semidirect product
            `G \rtimes_{id} Aut(G)`, where `id` is the identity function
            on `Aut(G)`, the automorphism group of `G`.

            See :wikipedia:`Holomorph (mathematics)`

            EXAMPLES::

                sage: G = Groups().example()
                sage: G.holomorph()
                Traceback (most recent call last):
                ...
                NotImplementedError: holomorph of General Linear Group of degree 4 over Rational Field not yet implemented
            """
            raise NotImplementedError("holomorph of %s not yet implemented"%self)

        def cayley_table(self, names='letters', elements=None):
            r"""
            Returns the "multiplication" table of this multiplicative group,
            which is also known as the "Cayley table".

            .. note:: The order of the elements in the row and column
              headings is equal to the order given by the table's
              :meth:`~sage.matrix.operation_table.OperationTable.column_keys`
              method.  The association between the actual elements and the
              names/symbols used in the table can also be retrieved as
              a dictionary with the
              :meth:`~sage.matrix.operation_table.OperationTable.translation`
              method.

            For groups, this routine should behave identically to the
            :meth:`~sage.categories.magmas.Magmas.ParentMethods.multiplication_table`
            method for magmas, which applies in greater generality.

            INPUT:

            - ``names`` - the type of names used, values are:

              * ``'letters'`` - lowercase ASCII letters are used
                for a base 26 representation of the elements'
                positions in the list given by :meth:`list`,
                padded to a common width with leading 'a's.
              * ``'digits'`` - base 10 representation of the
                elements' positions in the list given by
                :meth:`~sage.matrix.operation_table.OperationTable.column_keys`,
                padded to a common width with leading zeros.
              * ``'elements'`` - the string representations
                of the elements themselves.
              * a list - a list of strings, where the length
                of the list equals the number of elements.

            - ``elements`` - default = ``None``.  A list of
              elements of the group, in forms that can be
              coerced into the structure, eg. their string
              representations. This may be used to impose an
              alternate ordering on the elements, perhaps when
              this is used in the context of a particular structure.
              The default is to use whatever ordering is provided by the
              the group, which is reported by the
              :meth:`~sage.matrix.operation_table.OperationTable.column_keys`
              method.  Or the ``elements`` can be a subset
              which is closed under the operation. In particular,
              this can be used when the base set is infinite.

            OUTPUT:
            An object representing the multiplication table.  This is
            an :class:`~sage.matrix.operation_table.OperationTable` object
            and even more documentation can be found there.


            EXAMPLES:

            Permutation groups, matrix groups and abelian groups
            can all compute their multiplication tables.  ::

                sage: G = DiCyclicGroup(3)
                sage: T = G.cayley_table()
                sage: T.column_keys()
                ((), (5,6,7), (5,7,6)...(1,4,2,3)(5,7))
                sage: T
                *  a b c d e f g h i j k l
                 +------------------------
                a| a b c d e f g h i j k l
                b| b c a e f d i g h l j k
                c| c a b f d e h i g k l j
                d| d e f a b c j k l g h i
                e| e f d b c a l j k i g h
                f| f d e c a b k l j h i g
                g| g h i j k l d e f a b c
                h| h i g k l j f d e c a b
                i| i g h l j k e f d b c a
                j| j k l g h i a b c d e f
                k| k l j h i g c a b f d e
                l| l j k i g h b c a e f d

            ::

                sage: M=SL(2,2)
                sage: M.cayley_table()
                *  a b c d e f
                 +------------
                a| a b c d e f
                b| b a d c f e
                c| c f e b a d
                d| d e f a b c
                e| e d a f c b
                f| f c b e d a
                <BLANKLINE>

            ::

                sage: A=AbelianGroup([2,3])
                sage: A.cayley_table()
                *  a b c d e f
                 +------------
                a| a b c d e f
                b| b c a e f d
                c| c a b f d e
                d| d e f a b c
                e| e f d b c a
                f| f d e c a b

            Lowercase ASCII letters are the default symbols used
            for the table, but you can also specify the use of
            decimal digit strings, or provide your own strings
            (in the proper order if they have meaning).
            Also, if the elements themselves are not too complex,
            you can choose to just use the string representations
            of the elements themselves.  ::

                sage: C=CyclicPermutationGroup(11)
                sage: C.cayley_table(names='digits')
                 *  00 01 02 03 04 05 06 07 08 09 10
                  +---------------------------------
                00| 00 01 02 03 04 05 06 07 08 09 10
                01| 01 02 03 04 05 06 07 08 09 10 00
                02| 02 03 04 05 06 07 08 09 10 00 01
                03| 03 04 05 06 07 08 09 10 00 01 02
                04| 04 05 06 07 08 09 10 00 01 02 03
                05| 05 06 07 08 09 10 00 01 02 03 04
                06| 06 07 08 09 10 00 01 02 03 04 05
                07| 07 08 09 10 00 01 02 03 04 05 06
                08| 08 09 10 00 01 02 03 04 05 06 07
                09| 09 10 00 01 02 03 04 05 06 07 08
                10| 10 00 01 02 03 04 05 06 07 08 09

            ::

                sage: G=QuaternionGroup()
                sage: names=['1', 'I', '-1', '-I', 'J', '-K', '-J', 'K']
                sage: G.cayley_table(names=names)
                 *   1  I -1 -I  J -K -J  K
                  +------------------------
                 1|  1  I -1 -I  J -K -J  K
                 I|  I -1 -I  1  K  J -K -J
                -1| -1 -I  1  I -J  K  J -K
                -I| -I  1  I -1 -K -J  K  J
                 J|  J -K -J  K -1 -I  1  I
                -K| -K -J  K  J  I -1 -I  1
                -J| -J  K  J -K  1  I -1 -I
                 K|  K  J -K -J -I  1  I -1

            ::

                sage: A=AbelianGroup([2,2])
                sage: A.cayley_table(names='elements')
                    *      1    f1    f0 f0*f1
                     +------------------------
                    1|     1    f1    f0 f0*f1
                   f1|    f1     1 f0*f1    f0
                   f0|    f0 f0*f1     1    f1
                f0*f1| f0*f1    f0    f1     1

            The :meth:`~sage.matrix.operation_table.OperationTable.change_names`
            routine behaves similarly, but changes an existing table "in-place."
            ::

                sage: G=AlternatingGroup(3)
                sage: T=G.cayley_table()
                sage: T.change_names('digits')
                sage: T
                *  0 1 2
                 +------
                0| 0 1 2
                1| 1 2 0
                2| 2 0 1

            For an infinite group, you can still work with finite sets of
            elements, provided the set is closed under multiplication.
            Elements will be coerced into the group as part of setting
            up the table.  ::

                sage: G=SL(2,ZZ)
                sage: G
                Special Linear Group of degree 2 over Integer Ring
                sage: identity = matrix(ZZ, [[1,0], [0,1]])
                sage: G.cayley_table(elements=[identity, -identity])
                *  a b
                 +----
                a| a b
                b| b a

            The
            :class:`~sage.matrix.operation_table.OperationTable`
            class provides even greater flexibility, including changing
            the operation.  Here is one such example, illustrating the
            computation of commutators.  ``commutator`` is defined as
            a function of two variables, before being used to build
            the table. From this, the commutator subgroup seems obvious,
            and creating a Cayley table with just these three elements
            confirms that they form a closed subset in the group.
            ::

                sage: from sage.matrix.operation_table import OperationTable
                sage: G=DiCyclicGroup(3)
                sage: commutator = lambda x, y: x*y*x^-1*y^-1
                sage: T=OperationTable(G, commutator)
                sage: T
                .  a b c d e f g h i j k l
                 +------------------------
                a| a a a a a a a a a a a a
                b| a a a a a a c c c c c c
                c| a a a a a a b b b b b b
                d| a a a a a a a a a a a a
                e| a a a a a a c c c c c c
                f| a a a a a a b b b b b b
                g| a b c a b c a c b a c b
                h| a b c a b c b a c b a c
                i| a b c a b c c b a c b a
                j| a b c a b c a c b a c b
                k| a b c a b c b a c b a c
                l| a b c a b c c b a c b a
                sage: trans = T.translation()
                sage: comm = [trans['a'], trans['b'],trans['c']]
                sage: comm
                [(), (5,6,7), (5,7,6)]
                sage: P=G.cayley_table(elements=comm)
                sage: P
                *  a b c
                 +------
                a| a b c
                b| b c a
                c| c a b

            TODO:

            Arrange an ordering of elements into cosets of a normal
            subgroup close to size `\sqrt{n}`.  Then the quotient
            group structure is often apparent in the table.  See
            comments on Trac #7555.

            AUTHOR:

            - Rob Beezer (2010-03-15)

            """
            from sage.matrix.operation_table import OperationTable
            import operator
            return OperationTable(self, operation=operator.mul, names=names, elements=elements)

    class ElementMethods:
        ## inv(x), x/y
        pass

    Finite = LazyImport('sage.categories.finite_groups', 'FiniteGroups')
    #Algebras = LazyImport('sage.categories.group_algebras', 'GroupAlgebras')

    class Algebras(AlgebrasCategory):
        r"""
        The category of group algebras over a given base ring.

        EXAMPLES::

            sage: GroupAlgebras(IntegerRing())
            Category of group algebras over Integer Ring
            sage: GroupAlgebras(IntegerRing()).super_categories()
            [Category of hopf algebras with basis over Integer Ring,
             Category of monoid algebras over Integer Ring]

        Here is how to create the group algebra of a group `G`::

            sage: G = DihedralGroup(5)
            sage: QG = G.algebra(QQ); QG
            Group algebra of Dihedral group of order 10 as a permutation group over Rational Field

        and an example of computation::

            sage: g = G.an_element(); g
            (1,2,3,4,5)
            sage: (QG.term(g) + 1)**3
            B[()] + 3*B[(1,2,3,4,5)] + 3*B[(1,3,5,2,4)] + B[(1,4,2,5,3)]

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

        def example(self, G = None):
            """
            Return an example of group algebra.

            EXAMPLES::

                sage: GroupAlgebras(QQ[x]).example()
                Group algebra of Dihedral group of order 8 as a permutation group over Univariate Polynomial Ring in x over Rational Field

            An other group can be specified as optional argument::

                sage: GroupAlgebras(QQ).example(SymmetricGroup(4))
                Group algebra of Symmetric group of order 4! as a permutation group over Rational Field
            """
            from sage.groups.perm_gps.permgroup_named import DihedralGroup
            if G is None:
                G = DihedralGroup(4)
            return G.algebra(self.base_ring())

        class ParentMethods:

            def _repr_(self):
                r"""
                Return the string representation of `self`.

                EXAMPLES::

                    sage: A = Groups().example().algebra(QQ); A
                    Group algebra of General Linear Group of degree 4 over Rational Field over Rational Field
                    sage: A._name= "foo"
                    sage: A
                    foo over Rational Field
                """
                if hasattr(self, "_name"):
                    return self._name + " over {}".format(self.base_ring())
                else:
                    return 'Group algebra of {} over {}'.format(self.basis().keys(),
                                                                self.base_ring())

            def group(self):
                r"""
                Return the underlying group of the group algebra.

                EXAMPLES::

                    sage: GroupAlgebras(QQ).example(GL(3, GF(11))).group()
                    General Linear Group of degree 3 over Finite Field of size 11
                    sage: SymmetricGroupAlgebra(QQ,10).group()
                    Symmetric group of order 10! as a permutation group
                """
                return self.basis().keys()

            def algebra_generators(self):
                r"""
                Return generators of this group algebra (as an algebra).

                EXAMPLES::

                    sage: GroupAlgebras(QQ).example(SymmetricGroup(10)).algebra_generators()
                    Finite family {(1,2): B[(1,2)], (1,2,3,4,5,6,7,8,9,10): B[(1,2,3,4,5,6,7,8,9,10)]}

                .. NOTE::

                    This function is overloaded for SymmetricGroupAlgebras
                    to return Permutations and not Elements of the
                    symmetric group.
                """
                from sage.sets.family import Family
                return Family(self.group().gens(), self.term)

            def _conjugacy_classes_representatives_underlying_group(self):
                r"""
                Return a complete list of representatives of conjugacy
                classes of the underlying group.

                This works only for permutation groups. The ordering is
                that given by GAP.

                EXAMPLES::

                    sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3,4)]])
                    sage: SG = GroupAlgebras(QQ).example(G)
                    sage: SG._conjugacy_classes_representatives_underlying_group()
                    [(), (2,4), (1,2)(3,4), (1,2,3,4), (1,3)(2,4)]

                .. NOTE::

                    This function is overloaded for SymmetricGroupAlgebras to
                    return Permutations and not Elements of the symmetric group::

                    sage: SymmetricGroupAlgebra(ZZ,3)._conjugacy_classes_representatives_underlying_group()
                    [[2, 3, 1], [2, 1, 3], [1, 2, 3]]
                """
                return self.group().conjugacy_classes_representatives()

            def center(self):
                r"""
                Return the center of the group algebra.

                The canonical basis of the center of the group algebra
                is the family `(f_\sigma)_{\sigma\in C}`, where `C` is
                any collection of representatives of the conjugacy
                classes of the group, and `f_\sigma` is the sum of the
                elements in the conjugacy class of `\sigma`.

                OUTPUT:

                - A free module `V` indexed by conjugacy class
                  representatives of the group; its elements represent
                  formal linear combinations of the canonical basis
                  elements.

                .. WARNING::

                    - This method requires the underlying group to
                      have a method ``conjugacy_classes_representatives``
                      (every permutation group has one, thanks GAP!).
                    - The product has not been implemented yet.

                EXAMPLES::

                    sage: SymmetricGroupAlgebra(ZZ,3).center()
                    Free module generated by {[2, 3, 1], [2, 1, 3], [1, 2, 3]} over Integer Ring

                .. SEEALSO::

                    - :meth:`Groups.Algebras.ElementMethods.central_form`
                    - :meth:`Monoids.Algebras.ElementMethods.is_central`
                """
                I = self._conjugacy_classes_representatives_underlying_group()
                from sage.combinat.free_module import CombinatorialFreeModule
                return CombinatorialFreeModule(self.base_ring(), I)

            # Coalgebra structure

            def coproduct_on_basis(self, g):
                r"""
                Return the coproduct of the element ``g`` of the basis.

                Each basis element ``g`` is group-like. This method is
                used to compute the coproduct of any element.

                EXAMPLES::

                    sage: A=CyclicPermutationGroup(6).algebra(ZZ);A
                    Group algebra of Cyclic group of order 6 as a permutation group over Integer Ring
                    sage: g=CyclicPermutationGroup(6).an_element();g
                    (1,2,3,4,5,6)
                    sage: A.coproduct_on_basis(g)
                    B[(1,2,3,4,5,6)] # B[(1,2,3,4,5,6)]
                    sage: a=A.an_element();a
                    B[()] + 3*B[(1,2,3,4,5,6)] + 3*B[(1,3,5)(2,4,6)]
                    sage: a.coproduct()
                    B[()] # B[()] + 3*B[(1,2,3,4,5,6)] # B[(1,2,3,4,5,6)] + 3*B[(1,3,5)(2,4,6)] # B[(1,3,5)(2,4,6)]
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

                    sage: A=CyclicPermutationGroup(6).algebra(ZZ);A
                    Group algebra of Cyclic group of order 6 as a permutation group over Integer Ring
                    sage: g=CyclicPermutationGroup(6).an_element();g
                    (1,2,3,4,5,6)
                    sage: A.antipode_on_basis(g)
                    B[(1,6,5,4,3,2)]
                    sage: a=A.an_element();a
                    B[()] + 3*B[(1,2,3,4,5,6)] + 3*B[(1,3,5)(2,4,6)]
                    sage: a.antipode()
                    B[()] + 3*B[(1,5,3)(2,6,4)] + 3*B[(1,6,5,4,3,2)]
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

                    sage: A=CyclicPermutationGroup(6).algebra(ZZ);A
                    Group algebra of Cyclic group of order 6 as a permutation group over Integer Ring
                    sage: a=A.an_element();a
                    B[()] + 3*B[(1,2,3,4,5,6)] + 3*B[(1,3,5)(2,4,6)]
                    sage: a.counit()
                    7
                """
                return self.base_ring().sum(x.coefficients())

        class ElementMethods:

            def central_form(self):
                r"""
                Return ``self`` in the canonical basis of the center
                of the group algebra.

                INPUT:

                - ``self`` -- a central element of the group algebra

                OUTPUT:

                - A formal linear combination of the conjugacy class
                  representatives representing its coordinates in the
                  canonical basis of the center. See
                  :meth:`Groups.Algebras.ParentMethods.center` for
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

                .. SEEALSO::

                    - :meth:`Groups.Algebras.ParentMethods.center`
                    - :meth:`Monoids.Algebras.ElementMethods.is_central`
                """
                Z = self.parent().center()
                return sum(self[i] * Z.basis()[i] for i in Z.basis().keys())


