r"""
Additive semigroups
"""
#*****************************************************************************
#  Copyright (C) 2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom_singleton
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.algebra_functor import AlgebrasCategory
from sage.categories.additive_magmas import AdditiveMagmas

class AdditiveSemigroups(CategoryWithAxiom_singleton):
    """
    The category of additive semigroups.

    An *additive semigroup* is an associative :class:`additive magma
    <AdditiveMagmas>`, that is a set endowed with an operation `+`
    which is associative.

    EXAMPLES::

        sage: from sage.categories.additive_semigroups import AdditiveSemigroups
        sage: C = AdditiveSemigroups(); C
        Category of additive semigroups
        sage: C.super_categories()
        [Category of additive magmas]
        sage: C.all_super_categories()
        [Category of additive semigroups,
         Category of additive magmas,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

        sage: C.axioms()
        frozenset(['AdditiveAssociative'])
        sage: C is AdditiveMagmas().AdditiveAssociative()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (AdditiveMagmas, "AdditiveAssociative")

    AdditiveCommutative = LazyImport('sage.categories.commutative_additive_semigroups', 'CommutativeAdditiveSemigroups', at_startup=True)
    AdditiveUnital = LazyImport('sage.categories.additive_monoids', 'AdditiveMonoids', at_startup=True)

    class ParentMethods:
        def _test_additive_associativity(self, **options):
            r"""
            Test associativity for (not necessarily all) elements of this
            additive semigroup.

            INPUT:

            - ``options`` -- any keyword arguments accepted by :meth:`_tester`

            EXAMPLES:

            By default, this method tests only the elements returned by
            ``self.some_elements()``::

                sage: S = CommutativeAdditiveSemigroups().example()
                sage: S._test_additive_associativity()

            However, the elements tested can be customized with the
            ``elements`` keyword argument::

                sage: (a,b,c,d) = S.additive_semigroup_generators()
                sage: S._test_additive_associativity(elements = (a, b+c, d))

            See the documentation for :class:`TestSuite` for more information.
            """
            tester = self._tester(**options)
            S = tester.some_elements()
            from sage.combinat.cartesian_product import CartesianProduct
            for x,y,z in tester.some_elements(CartesianProduct(S,S,S)):
                tester.assert_((x + y) + z == x + (y + z))

    class CartesianProducts(CartesianProductsCategory):

        def extra_super_categories(self):
            """
            Implement the fact that a cartesian product of additive semigroups
            is an additive semigroup.

            EXAMPLES::

                sage: from sage.categories.additive_semigroups import AdditiveSemigroups
                sage: C = AdditiveSemigroups().CartesianProducts()
                sage: C.extra_super_categories()
                [Category of additive semigroups]
                sage: C.axioms()
                frozenset(['AdditiveAssociative'])
            """
            return [AdditiveSemigroups()]

    class Algebras(AlgebrasCategory):

        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: from sage.categories.additive_semigroups import AdditiveSemigroups
                sage: AdditiveSemigroups().Algebras(QQ).extra_super_categories()
                [Category of semigroups]
                sage: CommutativeAdditiveSemigroups().Algebras(QQ).super_categories()
                [Category of additive semigroup algebras over Rational Field,
                 Category of additive commutative additive magma algebras over Rational Field]
            """
            from sage.categories.semigroups import Semigroups
            return [Semigroups()]

        class ParentMethods:

            @cached_method
            def algebra_generators(self):
                r"""
                Return the generators of this algebra, as per
                :meth:`MagmaticAlgebras.ParentMethods.algebra_generators()
                <.magmatic_algebras.MagmaticAlgebras.ParentMethods.algebra_generators>`.

                They correspond to the generators of the additive semigroup.

                EXAMPLES::

                    sage: S = CommutativeAdditiveSemigroups().example(); S
                    An example of a commutative monoid: the free commutative monoid generated by ('a', 'b', 'c', 'd')
                    sage: A = S.algebra(QQ)
                    sage: A.algebra_generators()
                    Finite family {0: B[a], 1: B[b], 2: B[c], 3: B[d]}
                """
                return self.basis().keys().additive_semigroup_generators().map(self.monomial)

            def product_on_basis(self, g1, g2):
                r"""
                Product, on basis elements, as per
                :meth:`MagmaticAlgebras.WithBasis.ParentMethods.product_on_basis()
                <sage.categories.magmatic_algebras.MagmaticAlgebras.WithBasis.ParentMethods.product_on_basis>`.

                The product of two basis elements is induced by the
                addition of the corresponding elements of the group.

                EXAMPLES::

                    sage: S = CommutativeAdditiveSemigroups().example(); S
                    An example of a commutative monoid: the free commutative monoid generated by ('a', 'b', 'c', 'd')
                    sage: A = S.algebra(QQ)
                    sage: a,b,c,d = A.algebra_generators()
                    sage: a * b + b * d * c
                    B[c + b + d] + B[a + b]
                """
                return self.monomial(g1 + g2)

