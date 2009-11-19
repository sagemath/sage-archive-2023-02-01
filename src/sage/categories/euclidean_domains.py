r"""
EuclideanDomains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.categories.basic import PrincipalIdealDomains
from sage.misc.cachefunc import cached_method

class EuclideanDomains(Category):
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

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: EuclideanDomains().super_categories()
            [Category of principal ideal domains]
        """
        return [PrincipalIdealDomains()]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
