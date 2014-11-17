r"""
Filtered Algebras With Basis

A filtered algebra with basis over a commutative ring `R`
is a filtered algebra over `R` endowed with the structure
of a filtered module with basis (with the same underlying
filtered-module structure). See
:class:`~sage.categories.filtered_algebras.FilteredAlgebras` and
:class:`~sage.categories.filtered_modules_with_basis.FilteredModulesWithBasis`
for these two notions.
"""
#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.filtered_modules import FilteredModulesCategory

class FilteredAlgebrasWithBasis(FilteredModulesCategory):
    """
    The category of filtered algebras with a distinguished
    homogeneous basis.

    EXAMPLES::

        sage: C = AlgebrasWithBasis(ZZ).Filtered(); C
        Category of filtered algebras with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of algebras with basis over Integer Ring,
         Category of filtered algebras over Integer Ring,
         Category of filtered modules with basis over Integer Ring]

    TESTS::

        sage: TestSuite(C).run()
    """
    class ParentMethods:
        def graded_algebra(self):
            """
            Return the associated graded algebra to ``self``.

            See :class:`~sage.algebras.associated_graded.AssociatedGradedAlgebra`
            for the definition and the properties of this.

            EXAMPLES::

                sage: A = AlgebrasWithBasis(ZZ).Filtered().example()
                sage: A.graded_algebra()
                Graded Algebra of An example of a filtered algebra with basis:
                 the universal enveloping algebra of
                 Lie algebra of RR^3 with cross product over Integer Ring
            """
            from sage.algebras.associated_graded import AssociatedGradedAlgebra
            return AssociatedGradedAlgebra(self)

    class ElementMethods:

        def is_homogeneous(self):
            r"""
            Return whether ``self`` is homogeneous.

            EXAMPLES:

            Here is a case where the algebra is graded::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: (x, y) = (S[2], S[3])
                sage: (3*x).is_homogeneous()
                True
                sage: (x^3 - y^2).is_homogeneous()
                True
                sage: ((x + y)^2).is_homogeneous()
                False

            Let us now test a filtered algebra (but remember that the
            notion of homogeneity now depends on the choice of a
            basis)::

                sage: A = AlgebrasWithBasis(QQ).Filtered().example()
                sage: x,y,z = A.algebra_generators()
                sage: (x*y).is_homogeneous()
                True
                sage: (y*x).is_homogeneous()
                False
                sage: A.one().is_homogeneous()
                True
                sage: A.zero().is_homogeneous()
                True
                sage: (A.one()+x).is_homogeneous()
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
            The degree of a nonzero homogeneous element ``self`` in the
            filtered module.

            .. NOTE::

               This raises an error if the element is not homogeneous.
               To obtain the maximum of the degrees of the homogeneous
               summands, use :meth:`maximal_degree`.

            .. SEEALSO:: :meth:`maximal_degree`

            EXAMPLES:

            First, an example where the algebra is graded::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: (x, y) = (S[2], S[3])
                sage: x.homogeneous_degree()
                2
                sage: (x^3 + 4*y^2).homogeneous_degree()
                6
                sage: ((1 + x)^3).homogeneous_degree()
                Traceback (most recent call last):
                ...
                ValueError: element is not homogeneous

            Let us now test a filtered algebra (but remember that the
            notion of homogeneity now depends on the choice of a
            basis)::

                sage: A = AlgebrasWithBasis(QQ).Filtered().example()
                sage: x,y,z = A.algebra_generators()
                sage: (x*y).homogeneous_degree()
                2
                sage: (y*x).homogeneous_degree()
                Traceback (most recent call last):
                ...
                ValueError: element is not homogeneous
                sage: A.one().homogeneous_degree()
                0

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S.zero().degree()
                Traceback (most recent call last):
                ...
                ValueError: the zero element does not have a well-defined degree
            """
            if self.is_zero():
                raise ValueError("the zero element does not have a well-defined degree")
            if not self.is_homogeneous():
                raise ValueError("element is not homogeneous")
            return self.parent().degree_on_basis(self.leading_support())

        # default choice for degree; will be overridden as necessary
        degree = homogeneous_degree

        def maximal_degree(self):
            """
            The maximum of the degrees of the homogeneous components
            of ``self``.

            This is also the smallest `i` such that ``self`` belongs
            to `F_i`. Hence, it does not depend on the basis of the
            parent of ``self``.

            .. SEEALSO:: :meth:`homogeneous_degree`

            EXAMPLES:

            First, we test this on a graded algebra::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: (x, y) = (S[2], S[3])
                sage: x.maximal_degree()
                2
                sage: (x^3 + 4*y^2).maximal_degree()
                6
                sage: ((1 + x)^3).maximal_degree()
                6

            Let us now test a filtered algebra::

                sage: A = AlgebrasWithBasis(QQ).Filtered().example()
                sage: x,y,z = A.algebra_generators()
                sage: (x*y).maximal_degree()
                2
                sage: (y*x).maximal_degree()
                2
                sage: A.one().maximal_degree()
                0
                sage: A.zero().maximal_degree()
                Traceback (most recent call last):
                ...
                ValueError: the zero element does not have a well-defined degree
                sage: (A.one()+x).maximal_degree()
                1

            TESTS::

                sage: S = NonCommutativeSymmetricFunctions(QQ).S()
                sage: S.zero().degree()
                Traceback (most recent call last):
                ...
                ValueError: the zero element does not have a well-defined degree
            """
            if self.is_zero():
                raise ValueError("the zero element does not have a well-defined degree")
            degree_on_basis = self.parent().degree_on_basis
            return max(degree_on_basis(m) for m in self.support())

