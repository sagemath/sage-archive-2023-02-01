r"""
IntegralDomains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.misc.cachefunc import cached_method

class IntegralDomains(Category):
    """
    The category of integral domains
    commutative rings with no zero divisors

    EXAMPLES::

        sage: IntegralDomains()
        Category of integral domains
        sage: IntegralDomains().super_categories()
        [Category of commutative rings, Category of domains]

    TESTS::

        sage: TestSuite(IntegralDomains()).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: IntegralDomains().super_categories()
            [Category of commutative rings, Category of domains]
        """
        from sage.categories.basic import CommutativeRings, Domains
        return [CommutativeRings(), Domains()]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
