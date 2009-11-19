r"""
PrincipalIdealDomains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.categories.basic import GcdDomains
from sage.misc.cachefunc import cached_method

class PrincipalIdealDomains(Category):
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
      [Category of gcd domains]

    See also: http://en.wikipedia.org/wiki/Principal_ideal_domain

    TESTS::

        sage: TestSuite(PrincipalIdealDomains()).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: PrincipalIdealDomains().super_categories()
            [Category of gcd domains]
        """
        return [GcdDomains()]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
