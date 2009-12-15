r"""
CommutativeAdditiveGroups
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
        [Category of commutative additive groups, Category of commutative additive monoids, Category of commutative additive semigroups, Category of sets, Category of sets with partial maps, Category of objects]

    TESTS::

        sage: TestSuite(CommutativeAdditiveGroups()).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: CommutativeAdditiveGroups().super_categories()
            [Category of commutative additive monoids]
        """
        return [CommutativeAdditiveMonoids()]

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
            if left.parent() == right.parent() and hasattr(left, "_sub_"):
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
