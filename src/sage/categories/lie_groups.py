r"""
Lie Groups
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

#from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.groups import Groups
from sage.categories.manifolds import Manifolds

class LieGroups(Category_singleton):
    r"""
    The category of Lie groups.

    A Lie group is a topological group with a smooth manifold structure.

    EXAMPLES::

        sage: from sage.categories.lie_groups import LieGroups
        sage: C = LieGroups(); C
        Category of Lie groups

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.lie_groups import LieGroups
            sage: LieGroups().super_categories()
            [Category of topological groups, Category of smooth real manifolds]
        """
        return [Groups().Topological(), Manifolds().Real().Smooth()]

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of Lie groups defines no new
        structure: a morphism of topological spaces and of smooth
        manifolds is a morphism as Lie groups.

        .. SEEALSO:: :meth:`Category.additional_structure`

        EXAMPLES::

            sage: from sage.categories.lie_groups import LieGroups
            sage: LieGroups().additional_structure()
        """
        return None

    # Because Lie is a name that deserves to be capitalized
    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: from sage.categories.lie_groups import LieGroups
            sage: LieGroups() # indirect doctest
            Category of Lie groups
        """
        return "Lie groups"

