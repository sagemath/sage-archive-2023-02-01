r"""
Unique factorization domains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.categories.gcd_domains import GcdDomains

class UniqueFactorizationDomains(Category_singleton):
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

    def super_categories(self):
        """
        EXAMPLES::

            sage: UniqueFactorizationDomains().super_categories()
            [Category of gcd domains]
        """
        return [GcdDomains()]

    def is_structure_category(self):
        """
        Return whether ``self`` is a structure category.

        .. SEEALSO:: :meth:`Category.is_structure_category`

        The category of unique factorization domains does not define
        new structrure: a ring-morphism between unique factorization
        domains is a unique factorization domain-morphism.

        EXAMPLES::

            sage: UniqueFactorizationDomains().is_structure_category()
            False
        """
        return False

    class ParentMethods:
        def is_unique_factorization_domain(self):
            """
            Return True, since this in an object of the category of unique factorization domains.

            EXAMPLES::

                sage: Parent(QQ,category=UniqueFactorizationDomains()).is_unique_factorization_domain()
                True

            """
            return True

    class ElementMethods:
        # prime?
        # squareFree
        # factor
        pass
