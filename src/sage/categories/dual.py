"""
Dual functorial construction
"""
#*****************************************************************************
#  Copyright (C) 2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from category import Category

class DualityCategory(Category):
    """
    TESTS::

        sage: TestSuite(DualityCategory()).run() # mostly to silence sage -coverage on this abstract class
        Traceback (most recent call last):
        ...
        AssertionError: Not implemented method: dual
    """
    @abstract_method
    def dual(self):
        """
        Returns the dual category

        EXAMPLES::

            sage: Algebras(QQ).dual()
            Category of coalgebras over Rational Field
        """

# could do SelfDualCategory
