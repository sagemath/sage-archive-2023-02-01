"""
Lie Conformal Algebras With Basis

AUTHORS:

- Reimundo Heluani (2019-10-05): Initial implementation.
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.misc.abstract_method import abstract_method
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.super_modules import SuperModulesCategory

class LieConformalAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The subcategory of Lie conformal algebras with basis.

    EXAMPLES::

        sage: LieConformalAlgebras(QQbar).WithBasis()
        Category of Lie conformal algebras with basis over Algebraic Field
    """
    class ElementMethods:

        def index(self):
            """
            The index of this basis element.

            EXAMPLES::

                sage: V = lie_conformal_algebras.NeveuSchwarz(QQ)
                sage: V.inject_variables()
                Defining L, G, C
                sage: G.T(3).index()
                ('G', 3)
                sage: v = V.an_element(); v
                L + G + C
                sage: v.index()
                Traceback (most recent call last):
                ...
                ValueError: index can only be computed for monomials, got L + G + C
            """
            if self.is_zero():
                return None
            if not self.is_monomial():
                raise ValueError ("index can only be computed for "\
                                  "monomials, got {}".format(self))

            return next(iter(self.monomial_coefficients()))

    class Super(SuperModulesCategory):
        """
        The subcategory of super Lie conformal algebras with basis.

        EXAMPLES::

            sage: LieConformalAlgebras(AA).WithBasis().Super()
            Category of super Lie conformal algebras with basis over Algebraic Real Field
        """
        def extra_super_categories(self):
            """
            The extra super categories of this category

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).WithBasis().Super().super_categories()
                [Category of super modules with basis over Rational Field,
                 Category of Lie conformal algebras with basis over Rational Field,
                 Category of super Lie conformal algebras over Rational Field]
            """
            return [self.base_category(),]

    class Graded(GradedModulesCategory):
        """
        The subcategory of H-graded Lie conformal algebras with basis.

        EXAMPLES::

            sage: LieConformalAlgebras(QQbar).WithBasis().Graded()
            Category of H-graded Lie conformal algebras with basis over Algebraic Field
        """
        def _repr_object_names(self):
            """
            The names of the objects of this category

            EXAMPLES::

                sage: LieConformalAlgebras(QQbar).WithBasis().Graded()
                Category of H-graded Lie conformal algebras with basis over Algebraic Field
            """
            return "H-graded {}".format(self.base_category().\
                                        _repr_object_names())

        class Super(SuperModulesCategory):
            """
            The subcategory of super H-graded Lie conformal algebras
            with basis.

            EXAMPLES::

                sage: C = LieConformalAlgebras(QQbar).WithBasis()
                sage: C.Graded().Super()
                Category of super H-graded Lie conformal algebras with basis over Algebraic Field
                sage: C.Graded().Super() is C.Super().Graded()
                True
            """
            def extra_super_categories(self):
                """
                The extra super categories of this category.

                EXAMPLES::

                    sage: LieConformalAlgebras(QQ).WithBasis().Graded().Super().super_categories()
                    [Category of super Lie conformal algebras with basis over Rational Field,
                     Category of H-graded Lie conformal algebras with basis over Rational Field,
                     Category of super H-graded Lie conformal algebras over Rational Field]
                """
                return [self.base_category(),]

    class FinitelyGeneratedAsLieConformalAlgebra(CategoryWithAxiom_over_base_ring):
        """
        The subcategory of finitely generated Lie conformal
        algebras with basis.

        EXAMPLES::

            sage: C = LieConformalAlgebras(QQbar)
            sage: C.WithBasis().FinitelyGenerated()
            Category of finitely generated Lie conformal algebras with basis over Algebraic Field
            sage: C.WithBasis().FinitelyGenerated() is C.FinitelyGenerated().WithBasis()
            True
        """
        class Super(SuperModulesCategory):
            """
            The subcategory of super finitely generated Lie conformal
            algebras with basis.

            EXAMPLES::

                sage: LieConformalAlgebras(AA).WithBasis().FinitelyGenerated().Super()
                Category of super finitely generated Lie conformal algebras with basis over Algebraic Real Field
            """
            def extra_super_categories(self):
                """
                The extra super categories of this category

                EXAMPLES::

                    sage: C = LieConformalAlgebras(QQ).WithBasis().FinitelyGenerated().Super()
                    sage: C.super_categories()
                    [Category of super Lie conformal algebras with basis over Rational Field,
                     Category of finitely generated Lie conformal algebras with basis over Rational Field,
                     Category of super finitely generated Lie conformal algebras over Rational Field]
                """
                return [self.base_category(),]


        class Graded(GradedModulesCategory):
            """
            The subcategory of H-graded finitely generated Lie conformal
            algebras with basis.

            EXAMPLES::

                sage: LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated().Graded()
                Category of H-graded finitely generated Lie conformal algebras with basis over Algebraic Field
            """
            def _repr_object_names(self):
                """
                The names of the objects of this category

                EXAMPLES::

                    sage: LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated().Graded()
                    Category of H-graded finitely generated Lie conformal algebras with basis over Algebraic Field
                """
                return "H-graded {}".format(self.base_category().\
                                            _repr_object_names())

            class Super(SuperModulesCategory):
                """
                The subcategory of super H-graded finitely generated
                Lie conformal algebras with basis.

                EXAMPLES::

                    sage: C = LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated()
                    sage: C.Graded().Super()
                    Category of super H-graded finitely generated Lie conformal algebras with basis over Algebraic Field
                    sage: C.Graded().Super() is C.Super().Graded()
                    True
                """
                def extra_super_categories(self):
                    """
                    The extra super categories of this category.

                    EXAMPLES::

                        sage: C = LieConformalAlgebras(QQ).FinitelyGenerated().Graded().Super()
                        sage: C.super_categories()
                        [Category of super finitely generated Lie conformal algebras over Rational Field,
                         Category of super H-graded Lie conformal algebras over Rational Field,
                         Category of H-graded finitely generated Lie conformal algebras over Rational Field]
                    """
                    return [self.base_category(),]


