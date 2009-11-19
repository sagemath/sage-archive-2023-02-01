r"""
Domains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.rings import Rings
from sage.categories.category import Category
from sage.misc.cachefunc import cached_method

class Domains(Category):
    """
    The category of domains

    An domain (or non-commutative integral domains), is a not
    necessarily commutative ring which has no zero divisors.

    EXAMPLES::

        sage: Domains()
        Category of domains
        sage: Domains().super_categories()
        [Category of rings]

    TESTS::

        sage: TestSuite(Domains()).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Domains().super_categories()
            [Category of rings]
        """
        return [Rings()]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
