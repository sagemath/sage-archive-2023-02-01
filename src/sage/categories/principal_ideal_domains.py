r"""
Principal ideal domains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains

class PrincipalIdealDomains(Category_singleton):
    """
    The category of (constructive) principal ideal domains

    By constructive, we mean that a single generator can be
    constructively found for any ideal given by a finite set of
    generators. Note that this constructive definition only implies
    that finitely generated ideals are principal. It is not clear what
    we would mean by an infinitely generated ideal.

    EXAMPLES::

      sage: PrincipalIdealDomains()
      Category of principal ideal domains
      sage: PrincipalIdealDomains().super_categories()
      [Category of unique factorization domains]

    See also: http://en.wikipedia.org/wiki/Principal_ideal_domain

    TESTS::

        sage: TestSuite(PrincipalIdealDomains()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: PrincipalIdealDomains().super_categories()
            [Category of unique factorization domains]
        """
        return [UniqueFactorizationDomains()]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
