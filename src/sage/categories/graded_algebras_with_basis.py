r"""
Graded algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import GradedAlgebras, GradedModulesWithBasis, AlgebrasWithBasis
from sage.misc.cachefunc import cached_method

class GradedAlgebrasWithBasis(Category_over_base_ring):
    """
    The category of graded algebras with a distinguished basis

    EXAMPLES::

        sage: GradedAlgebrasWithBasis(ZZ)
        Category of graded algebras with basis over Integer Ring
        sage: GradedAlgebrasWithBasis(ZZ).super_categories()
        [Category of graded modules with basis over Integer Ring, Category of graded algebras over Integer Ring, Category of algebras with basis over Integer Ring]

    TESTS::

        sage: TestSuite(GradedAlgebrasWithBasis(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: GradedAlgebrasWithBasis(QQ).super_categories()
            [Category of graded modules with basis over Rational Field, Category of graded algebras over Rational Field, Category of algebras with basis over Rational Field]
        """
        R = self.base_ring()
        return [GradedModulesWithBasis(R),GradedAlgebras(R), AlgebrasWithBasis(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        def is_homogeneous(self):
            """
            Return whether this element is homogeneous.

            EXAMPLES::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: (x, y) = (S[2], S[3])
                sage: (3*x).is_homogeneous()
                True
                sage: (x^3 - y^2).is_homogeneous()
                True
                sage: ((x + y)^2).is_homogeneous()
                False
            """
            degree_on_basis = self.parent().degree_on_basis
            degree = None
            for m in self.support():
                if degree is None:
                    degree = degree_on_basis(m)
                else:
                    if degree != degree_on_basis(m):
                        return False
            return True

        def homogeneous_degree(self):
            """
            The degree of this element.

            .. note::

               This raises an error if the element is not homogeneous.
               To obtain the maximum of the degrees of the homogeneous
               summands, use :meth:`maximal_degree`

            .. seealso: :meth:`maximal_degree`

            EXAMPLES::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: (x, y) = (S[2], S[3])
                sage: x.homogeneous_degree()
                2
                sage: (x^3 + 4*y^2).homogeneous_degree()
                6
                sage: ((1 + x)^3).homogeneous_degree()
                Traceback (most recent call last):
                ...
                ValueError: Element is not homogeneous.

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S.zero().degree()
                Traceback (most recent call last):
                ...
                ValueError: The zero element does not have a well-defined degree.
            """
            if self.is_zero():
                raise ValueError("The zero element does not have a well-defined degree.")
            try:
                assert self.is_homogeneous()
                return self.parent().degree_on_basis(self.leading_support())
            except AssertionError:
                raise ValueError("Element is not homogeneous.")

        # default choice for degree; will be overridden as necessary
        degree = homogeneous_degree

        def maximal_degree(self):
            """
            The maximum of the degrees of the homogeneous summands.

            .. seealso: :meth:`homogeneous_degree`

            EXAMPLES::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: (x, y) = (S[2], S[3])
                sage: x.maximal_degree()
                2
                sage: (x^3 + 4*y^2).maximal_degree()
                6
                sage: ((1 + x)^3).maximal_degree()
                6

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S.zero().degree()
                Traceback (most recent call last):
                ...
                ValueError: The zero element does not have a well-defined degree.
            """
            if self.is_zero():
                raise ValueError("The zero element does not have a well-defined degree.")
            else:
                degree_on_basis = self.parent().degree_on_basis
                return max(degree_on_basis(m) for m in self.support())
