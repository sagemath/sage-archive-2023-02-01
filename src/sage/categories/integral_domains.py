r"""
Integral domains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2012 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.domains import Domains

class IntegralDomains(CategoryWithAxiom):
    """
    The category of integral domains

    An integral domain is commutative ring with no zero divisors, or
    equivalently a commutative domain.

    EXAMPLES::

        sage: C = IntegralDomains(); C
        Category of integral domains
        sage: sorted(C.super_categories(), key=str)
        [Category of commutative rings, Category of domains]
        sage: C is Domains().Commutative()
        True
        sage: C is Rings().Commutative().NoZeroDivisors()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (Domains, "Commutative")

    class ParentMethods:
        def is_integral_domain(self):
            """
            Return True, since this in an object of the category of integral domains.

            EXAMPLES::

                sage: QQ.is_integral_domain()
                True
                sage: Parent(QQ,category=IntegralDomains()).is_integral_domain()
                True

            """
            return True

    class ElementMethods:
        pass
