r"""
Additive groups
"""
#*****************************************************************************
#  Copyright (C) 2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom_singleton
from sage.categories.additive_monoids import AdditiveMonoids
from sage.structure.sage_object import have_same_parent

class AdditiveGroups(CategoryWithAxiom_singleton):
    r"""
    The category of additive groups.

    An *additive group* is a set with an internal binary operation `+` which
    is associative, admits a zero, and where every element can be negated.

    EXAMPLES::

        sage: from sage.categories.additive_groups import AdditiveGroups
        sage: from sage.categories.additive_monoids import AdditiveMonoids
        sage: AdditiveGroups()
        Category of additive groups
        sage: AdditiveGroups().super_categories()
        [Category of additive inverse additive unital additive magmas,
         Category of additive monoids]
        sage: AdditiveGroups().all_super_categories()
        [Category of additive groups,
         Category of additive inverse additive unital additive magmas,
         Category of additive monoids,
         Category of additive unital additive magmas,
         Category of additive semigroups,
         Category of additive magmas,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

        sage: AdditiveGroups().axioms()
        frozenset(['AdditiveAssociative', 'AdditiveUnital', 'AdditiveInverse'])
        sage: AdditiveGroups() is AdditiveMonoids().AdditiveInverse()
        True

    TESTS::

        sage: C = AdditiveGroups()
        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (AdditiveMonoids, "AdditiveInverse")

    AdditiveCommutative = LazyImport('sage.categories.commutative_additive_groups', 'CommutativeAdditiveGroups', at_startup=True)
