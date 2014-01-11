r"""
Additive monoids
"""
#*****************************************************************************
#  Copyright (C) 2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom_singleton
from sage.categories.additive_semigroups import AdditiveSemigroups

class AdditiveMonoids(CategoryWithAxiom_singleton):
    """
    The category of additive monoids, i.e. sets with an associative
    operation ``+`` which is associative and admits a zero.

    EXAMPLES::

        sage: from sage.categories.additive_monoids import AdditiveMonoids
        sage: C = AdditiveMonoids(); C
        Category of additive monoids
        sage: C.super_categories()
        [Category of additive unital additive magmas, Category of additive semigroups]
        sage: sorted(C.axioms())
        ['AdditiveAssociative', 'AdditiveUnital']
        sage: from sage.categories.additive_semigroups import AdditiveSemigroups
        sage: C is AdditiveSemigroups().AdditiveUnital()
        True

    TESTS::

        sage: C.Algebras(QQ).is_subcategory(AlgebrasWithBasis(QQ))
        True
        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (AdditiveSemigroups, "AdditiveUnital")

    AdditiveCommutative = LazyImport('sage.categories.commutative_additive_monoids', 'CommutativeAdditiveMonoids', at_startup=True)
    AdditiveInverse = LazyImport('sage.categories.additive_groups', 'AdditiveGroups', at_startup=True)

    class ParentMethods:
        def sum(self, args):
            r"""
            n-ary sum

            INPUT:

             - ``args`` -- a list (or iterable) of elements of ``self``

            Returns the sum of the elements in `args`, as an element of `self`.

            EXAMPLES::

                sage: S = CommutativeAdditiveMonoids().example()
                sage: (a,b,c,d) = S.additive_semigroup_generators()
                sage: S.sum((a,b,a,c,a,b))
                3*a + c + 2*b
                sage: S.sum(())
                0
                sage: S.sum(()).parent() == S
                True
            """
            return sum(args, self.zero())
