r"""
Domains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2012 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.rings import Rings

class Domains(CategoryWithAxiom):
    """
    The category of domains

    A domain (or non-commutative integral domain), is a ring, not
    necessarily commutative, with no nonzero zero divisors.

    EXAMPLES::

        sage: C = Domains(); C
        Category of domains
        sage: C.super_categories()
        [Category of rings]
        sage: C is Rings().NoZeroDivisors()
        True

    TESTS::

        sage: TestSuite(C).run()
    """

    _base_category_class_and_axiom = (Rings, "NoZeroDivisors")

    def super_categories(self):
        """
        EXAMPLES::

            sage: Domains().super_categories()
            [Category of rings]
        """
        return [Rings()]

    Commutative = LazyImport('sage.categories.integral_domains', 'IntegralDomains', at_startup=True)

    class ParentMethods:
        pass

    class ElementMethods:
        pass
