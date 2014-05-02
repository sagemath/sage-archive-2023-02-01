r"""
Commutative additive groups
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_types import AbelianCategory
from sage.categories.commutative_additive_monoids import CommutativeAdditiveMonoids
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.structure.sage_object import have_same_parent

class CommutativeAdditiveGroups(AbelianCategory):
    """
    The category of abelian groups, i.e. additive abelian monoids
    where each element has an inverse.

    EXAMPLES::

        sage: CommutativeAdditiveGroups()
        Category of commutative additive groups
        sage: CommutativeAdditiveGroups().super_categories()
        [Category of commutative additive monoids]
        sage: CommutativeAdditiveGroups().all_super_categories()
        [Category of commutative additive groups, Category of commutative additive monoids, Category of commutative additive semigroups, Category of additive magmas, Category of sets, Category of sets with partial maps, Category of objects]

    TESTS::

        sage: TestSuite(CommutativeAdditiveGroups()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: CommutativeAdditiveGroups().super_categories()
            [Category of commutative additive monoids]
        """
        return [CommutativeAdditiveMonoids()]

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
                return self.parent()(tuple([x+y for x,y in zip(self.summand_split(),right.summand_split())]))

            def _neg_(self):
                r"""
                EXAMPLES::

                    sage: G=GF(5); GG = G.cartesian_product(G)
                    sage: oneone = GG([GF(5)(1),GF(5)(1)])
                    sage: -oneone
                    (4, 4)
                """
                return self.parent()(tuple([-x for x in self.summand_split()]))

        class ParentMethods:
            def zero(self):
                r"""
                Returns the zero of this group

                EXAMPLE::

                    sage: GF(8,'x').cartesian_product(GF(5)).zero()
                    (0, 0)
                """
                return self(tuple([0 for _ in self._sets]))

    class ParentMethods:
        pass

    class ElementMethods:
        ##def -x, -(x,y):
        def __sub__(left, right):
            """
            Top-level subtraction operator
            See extensive documentation at the top of element.pyx.

            EXAMPLES::

                sage: F = CombinatorialFreeModule(QQ, ['a','b'])
                sage: a,b = F.basis()
                sage: a - b
                B['a'] - B['b']
            """
            if have_same_parent(left, right) and hasattr(left, "_sub_"):
                return left._sub_(right)
            from sage.structure.element import get_coercion_model
            import operator
            return get_coercion_model().bin_op(left, right, operator.sub)

        ##################################################
        # Negation
        ##################################################

        def __neg__(self):
            """
            Top-level negation operator for elements of abelian
            monoids, which may choose to implement _neg_ rather than
            __neg__ for consistancy.

            EXAMPLES::

                sage: F = CombinatorialFreeModule(QQ, ['a','b'])
                sage: a,b = F.basis()
                sage: - b
                -B['b']
            """
            return self._neg_()
