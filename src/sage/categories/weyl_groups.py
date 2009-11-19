r"""
Weyl Groups
"""
#*****************************************************************************
#  Copyright (C) 2009    Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category import Category
from sage.categories.coxeter_groups import CoxeterGroups
import sage.combinat.ranker

class WeylGroups(Category):
    r"""
    The category of Weyl groups

    See: http://en.wikipedia.org/wiki/Weyl_groups

    EXAMPLES::

        sage: WeylGroups()                      # todo: uppercase for Weyl
        Category of weyl groups
        sage: WeylGroups().super_categories()
        [Category of coxeter groups]

    Here are some examples::

        sage: WeylGroups().example()            # todo: not implemented
        sage: FiniteWeylGroups().example()
        The symmetric group on {0, ..., 3}
        sage: AffineWeylGroups().example()      # todo: not implemented
        sage: WeylGroup(["B", 3])
        Weyl Group of type ['B', 3] (as a matrix group acting on the ambient space)

    This one will eventually be also in this category::

        sage: SymmetricGroup(4)
        Symmetric group of order 4! as a permutation group

    TESTS::

        sage: C = WeylGroups()
        sage: TestSuite(C).run()
    """

    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: WeylGroups().super_categories()
            [Category of coxeter groups]
        """
        return [CoxeterGroups()]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
