r"""
PartiallyOrderedSets
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.categories.basic import Sets
from sage.misc.cachefunc import cached_method

class PartiallyOrderedSets(Category):
    """
    The category of partially ordered sets


    EXAMPLES::

        sage: PartiallyOrderedSets()
        Category of partially ordered sets
        sage: PartiallyOrderedSets().super_categories()
        [Category of sets]

    TESTS::

        sage: TestSuite(PartiallyOrderedSets()).run()

    TODO: add appropriate abstract_methods to specify the operations
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: PartiallyOrderedSets().super_categories()
            [Category of sets]
        """
        return [Sets()]

    class ParentMethods:
        pass

    class ElementMethods:
        ##def <(x,y):
        pass
