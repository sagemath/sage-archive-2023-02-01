r"""
Lie Groups
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import LazyImport
from sage.categories.category import Category
from sage.categories.category_types import Category_over_base_ring
from sage.categories.sets_cat import Sets

class Manifolds(Category_over_base_ring):
    r"""
    The category of Lie groups.

    A Lie group is a topological group with a smooth differentiable
    manifold structure.

    INPUT:

    - ``k`` -- (default: ``RR``) the field `k`

    EXAMPLES::

        sage: from sage.categories.lie_groups import LieGroups
        sage: C = LieGroups(); C
        Category of manifolds over Real Field
        sage: C.super_categories()
        [Category of additive commutative additive associative additive unital distributive magmas and additive magmas,
         Category of modules over Integer Ring]
        sage: C = LieGroups(CC); C
        Category of manifolds over Complex Field

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.lie_groups import LieGroups
            sage: LieGroups().super_categories()
        """
        R = self.base_ring()
        # TODO: Make this smooth differentiable manifolds
        return [Groups().Topological(), Manifolds(R)]

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of Lie groups defines no new
        structure: a morphism of topological spaces and of smooth
        differentiable manifolds is a morphism as Lie groups.

        .. SEEALSO:: :meth:`Category.additional_structure`

        EXAMPLES::

            sage: from sage.categories.lie_groups import LieGroups
            sage: LieGroups().additional_structure()
        """
        return None

