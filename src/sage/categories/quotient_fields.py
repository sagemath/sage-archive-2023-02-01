r"""
QuotientFields
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.categories.all import Fields
from sage.misc.cachefunc import cached_method

class QuotientFields(Category):
    """
    The category of quotient fields over an integral domain

    EXAMPLES::

        sage: QuotientFields()
        Category of quotient fields
        sage: QuotientFields().super_categories()
        [Category of fields]

    TESTS::

        sage: TestSuite(QuotientFields()).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: QuotientFields().super_categories()
            [Category of fields]
        """
        return [Fields()]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
