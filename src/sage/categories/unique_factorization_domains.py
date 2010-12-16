r"""
Unique factorization domains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.misc.cachefunc import cached_method

class UniqueFactorizationDomains(Category):
    """
    The category of unique factorization domains
    constructive unique factorization domains, i.e. where one can constructively
    factor members into a product of a finite number of irreducible elements

    EXAMPLES::

        sage: UniqueFactorizationDomains()
        Category of unique factorization domains
        sage: UniqueFactorizationDomains().super_categories()
        [Category of gcd domains]

    TESTS::

        sage: TestSuite(UniqueFactorizationDomains()).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: UniqueFactorizationDomains().super_categories()
            [Category of gcd domains]
        """
        from sage.categories.gcd_domains import GcdDomains
        return [GcdDomains()]

    class ParentMethods:
        pass

    class ElementMethods:
        # prime?
        # squareFree
        # factor
        pass
