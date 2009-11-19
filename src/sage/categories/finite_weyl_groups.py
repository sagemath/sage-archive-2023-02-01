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
        running ._test_an_element() ... done
        running ._test_associativity() ... done
        running ._test_element_pickling() ... done
        running ._test_inverse() ... done
        running ._test_one() ... done
        running ._test_pickling() ... done
        running ._test_prod() ... done
        running ._test_reduced_word() ... done
        running ._test_some_elements() ... done
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
