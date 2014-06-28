r"""
Permutation groups
"""
#*****************************************************************************
#  Copyright (C) 2011      Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category import Category
from sage.categories.groups import Groups
from sage.misc.lazy_import import LazyImport

class PermutationGroups(Category):
    r"""
    The category of permutation groups.

    A *permutation group* is a group whose elements are concretely
    represented by permutations of some set.

    .. TODO:: shall we accept only permutations with finite support or not?

    EXAMPLES::

        sage: PermutationGroups()
        Category of permutation groups
        sage: PermutationGroups().super_categories()
        [Category of groups]

    TESTS::

        sage: C = PermutationGroups()
        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        Return a list of the immediate super categories of ``self``.

        EXAMPLES::

            sage: PermutationGroups().super_categories()
            [Category of groups]
        """
        return [Groups()]

    Finite = LazyImport('sage.categories.finite_permutation_groups', 'FinitePermutationGroups')
