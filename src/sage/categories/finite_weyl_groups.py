r"""
Finite Weyl Groups
"""
#*****************************************************************************
#  Copyright (C) 2009    Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category import Category
from sage.categories.finite_coxeter_groups import FiniteCoxeterGroups
from sage.categories.weyl_groups import WeylGroups

class FiniteWeylGroups(Category):
    """
    The category of finite Weyl groups.

    EXAMPLES::

        sage: C = FiniteWeylGroups()
        sage: C
        Category of finite weyl groups
        sage: C.super_categories()
        [Category of weyl groups, Category of finite coxeter groups]
        sage: C.example()
        The symmetric group on {0, ..., 3}

    TESTS::

        sage: W = FiniteWeylGroups().example()
        sage: TestSuite(W).run(verbose = "True")
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_element_pickling() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_has_descent() . . . pass
        running ._test_inverse() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_reduced_word() . . . pass
        running ._test_simple_projections() . . . pass
        running ._test_some_elements() . . . pass
    """

    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: FiniteWeylGroups().super_categories()
            [Category of weyl groups, Category of finite coxeter groups]
        """
        return [WeylGroups(), FiniteCoxeterGroups()]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
