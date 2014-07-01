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
        def _test_zero_divisors(self, **options):
            """
            Check to see that there are no zero divisors.

            EXAMPLES::

                sage: ZZ._test_zero_divisors()
                sage: Zp(5)._test_zero_divisors()
            """
            # This try-except is probably not needed as is_exact()
            #   should be a method for every commutative ring
            try:
                if not self.is_exact():
                    return # Can't check on inexact rings
            except AttributeError:
                # Assume a ring is exact if it doesn't have the method
                pass

            tester = self._tester(**options)

            # Filter out zero
            S = [s for s in tester.some_elements() if not s.is_zero()]

            from sage.combinat.cartesian_product import CartesianProduct
            for a,b in tester.some_elements(CartesianProduct(S,S)):
                p = a * b
                tester.assertFalse(p.is_zero())

    class ElementMethods:
        pass
