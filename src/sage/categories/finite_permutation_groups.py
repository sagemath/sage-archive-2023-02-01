r"""
FinitePermutationGroups
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.misc import attrcall
from sage.categories.category import Category
from sage.categories.finite_groups import FiniteGroups

class FinitePermutationGroups(Category):
    r"""
    The category of finite permutation groups, i.e. groups concretely
    represented as groups of permutations acting on a finite set.

    EXAMPLES::

        sage: FinitePermutationGroups()
        Category of finite permutation groups
        sage: FinitePermutationGroups().super_categories()
        [Category of finite groups]

        sage: FinitePermutationGroups().example()
        Dihedral group of order 6 as a permutation group

    TESTS::

        sage: C = FinitePermutationGroups()
        sage: TestSuite(C).run(verbose = True)
        running ._test_category() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass

        sage: G = FinitePermutationGroups().example()
        sage: TestSuite(G).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_inverse() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
    """
    @cached_method
    def super_categories(self):
        r"""
        Returns a list of the (immediate) super categories of ``self``.

        EXAMPLES::

            sage: FinitePermutationGroups().super_categories()
            [Category of finite groups]
        """
        return [FiniteGroups()]

    def example(self):
        """
        Returns an example of finite permutation group, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: G = FinitePermutationGroups().example(); G
            Dihedral group of order 6 as a permutation group
        """
        from sage.groups.perm_gps.permgroup_named import DihedralGroup
        return DihedralGroup(3)

    class ParentMethods:
        # TODO (Cf MuPAD-Combinat, lib/DOMAINS/CATEGORIES/PermutationGroup.mu)
        #  - cycle_indicator
        #  - strong_generating_set
        #  - ...
        pass

    class ElementMethods:
        # TODO: put abstract_methods for
        #  - cycle_type
        #  - orbit
        #  - ...
        pass
