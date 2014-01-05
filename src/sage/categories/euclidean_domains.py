r"""
Euclidean domains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.categories.principal_ideal_domains import PrincipalIdealDomains

class EuclideanDomains(Category_singleton):
    """
    The category of euclidean domains
    constructive euclidean domain, i.e. one can divide producing a quotient and a
    remainder where the remainder is either zero or is "smaller" than the divisor

    EXAMPLES::

      sage: EuclideanDomains()
      Category of euclidean domains
      sage: EuclideanDomains().super_categories()
      [Category of principal ideal domains]

    TESTS::

        sage: TestSuite(EuclideanDomains()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: EuclideanDomains().super_categories()
            [Category of principal ideal domains]
        """
        return [PrincipalIdealDomains()]

    class ParentMethods:
        def is_euclidean_domain(self):
            """
            Return True, since this in an object of the category of Euclidean domains.

            EXAMPLES::

                sage: Parent(QQ,category=EuclideanDomains()).is_euclidean_domain()
                True

            """
            return True

    class ElementMethods:
        pass
