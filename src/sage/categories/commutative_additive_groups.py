r"""
Commutative additive groups
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_types import AbelianCategory
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.algebra_functor import AlgebrasCategory
from sage.categories.additive_groups import AdditiveGroups

class CommutativeAdditiveGroups(CategoryWithAxiom, AbelianCategory):
    """
    The category of abelian groups, i.e. additive abelian monoids
    where each element has an inverse.

    EXAMPLES::

        sage: C = CommutativeAdditiveGroups(); C
        Category of commutative additive groups
        sage: C.super_categories()
        [Category of additive groups, Category of commutative additive monoids]
        sage: sorted(C.axioms())
        ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveInverse', 'AdditiveUnital']
        sage: C is CommutativeAdditiveMonoids().AdditiveInverse()
        True
        sage: from sage.categories.additive_groups import AdditiveGroups
        sage: C is AdditiveGroups().AdditiveCommutative()
        True

    .. NOTE::

        This category is currently empty. It's left there for backward
        compatibility and because it is likely to grow in the future.

    TESTS::

        sage: TestSuite(CommutativeAdditiveGroups()).run()
    """
    _base_category_class_and_axiom = (AdditiveGroups, "AdditiveCommutative")

    class Algebras(AlgebrasCategory):
        pass

    class CartesianProducts(CartesianProductsCategory):

        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
                sage: CommutativeAdditiveGroups().CartesianProducts().extra_super_categories();
                [Category of commutative additive groups]
                sage: CommutativeAdditiveGroups().CartesianProducts().super_categories()
                [Category of commutative additive groups, Category of Cartesian products of sets]
            """
            return [CommutativeAdditiveGroups()]

        class ElementMethods:
            # TODO: move to AdditiveMagmas.CartesianProducts.ElementMethods
            def _add_(self, right):
                r"""
                EXAMPLES::

                    sage: G5=GF(5); G8=GF(4,'x'); GG = G5.cartesian_product(G8)
                    sage: e = GG((G5(1),G8.primitive_element())); e
                    (1, x)
                    sage: e+e
                    (2, 0)
                    sage: e=groups.misc.AdditiveCyclic(8)
                    sage: x=e.cartesian_product(e)((e(1),e(2)))
                    sage: x
                    (1, 2)
                    sage: 4*x
                    (4, 0)
                """
                return self.parent()._cartesian_product_of_elements(
                    x+y for x,y in zip(self.cartesian_factors(),
                                       right.cartesian_factors()))

            # TODO: move to AdditiveGroups.CartesianProducts.ElementMethods
            def _neg_(self):
                r"""
                EXAMPLES::

                    sage: G=GF(5); GG = G.cartesian_product(G)
                    sage: oneone = GG([GF(5)(1),GF(5)(1)])
                    sage: -oneone
                    (4, 4)
                """
                return self.parent()._cartesian_product_of_elements(
                    -x for x in self.cartesian_factors())

        class ParentMethods:
            # TODO: move to AdditiveMagmas.AdditiveUnital.CartesianProducts.ElementMethods
            def zero(self):
                r"""
                Returns the zero of this group

                EXAMPLE::

                    sage: GF(8,'x').cartesian_product(GF(5)).zero()
                    (0, 0)
                """
                return self._cartesian_product_of_elements(
                    _.zero() for _ in self._sets)
