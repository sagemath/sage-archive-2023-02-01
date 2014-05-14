r"""
Gcd domains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.categories.integral_domains import IntegralDomains

class GcdDomains(Category_singleton):
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

    def super_categories(self):
        """
        EXAMPLES::

            sage: GcdDomains().super_categories()
            [Category of integral domains]
        """
        return [IntegralDomains()]

    def is_structure_category(self):
        """
        Return whether ``self`` is a structure category.

        .. SEEALSO:: :meth:`Category.is_structure_category`


        The category of gcd domains does not define new structrure: a
        ring-morphism between gcd domains is a gcd domain-morphism.

        EXAMPLES::

            sage: GcdDomains().is_structure_category()
            False
        """
        return False

    class ParentMethods:
        pass

    class ElementMethods:
        # gcd(x,y)
        # lcm(x,y)
        pass
