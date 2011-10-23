r"""
Integral domains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.categories.category_singleton import Category_singleton
from sage.misc.cachefunc import cached_method
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.domains import Domains

class IntegralDomains(Category_singleton):
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

    def super_categories(self):
        """
        EXAMPLES::

            sage: IntegralDomains().super_categories()
            [Category of commutative rings, Category of domains]
        """
        return [CommutativeRings(), Domains()]

    class ParentMethods:
        def is_integral_domain(self):
            """
            Return True, since this in an object of the category of integral domains.

            EXAMPLES::

                sage: Parent(QQ,category=IntegralDomains()).is_integral_domain()
                True

            """
            return True

    class ElementMethods:
        pass
