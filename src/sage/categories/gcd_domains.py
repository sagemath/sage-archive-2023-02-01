r"""
GcdDomains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.misc.cachefunc import cached_method

class GcdDomains(Category):
    """
    The category of gcd domains
    domains where gcd can be computed but where there is no guarantee of
    factorisation into irreducibles

    EXAMPLES::

        sage: GcdDomains()
        Category of gcd domains
        sage: GcdDomains().super_categories()
        [Category of integral domains]

    TESTS::

        sage: TestSuite(GcdDomains()).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GcdDomains().super_categories()
            [Category of integral domains]
        """
        from sage.categories.integral_domains import IntegralDomains
        return [IntegralDomains()]

    class ParentMethods:
        pass

    class ElementMethods:
        # gcd(x,y)
        # lcm(x,y)
        pass
