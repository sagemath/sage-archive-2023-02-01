r"""
Euclidean domains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.categories.category_singleton import Category_singleton
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.misc.cachefunc import cached_method
from sage.structure.element import coerce_binop

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
        @coerce_binop
        def gcd(self, other):
            """
            Return the greatest common divisor of this element and ``other``.

            INPUT:

                - ``other`` -- an element in the same ring as ``self``

            ALGORITHM:

            Algorithm 3.2.1 in [Coh1996].

            REFERENCES:

            .. [Coh1996] Henri Cohen. A Course in Computational Algebraic Number Theory. Springer, 1996.

            EXAMPLES::

                sage: EuclideanDomains().ElementMethods().gcd(6,4)
                2

            """
            A = self
            B = other
            while not B.is_zero():
                Q, R = A.quo_rem(B)
                A = B
                B = R
            return A
