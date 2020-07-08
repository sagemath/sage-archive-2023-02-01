"""
Finitely Generated Lie Conformal Algebras

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
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.super_modules import SuperModulesCategory

class FinitelyGeneratedAsLieConformalAlgebra(CategoryWithAxiom_over_base_ring):
    """
    The category of finitely generated Lie conformal algebras.

    EXAMPLES::

        sage: LieConformalAlgebras(QQbar).FinitelyGenerated()
        Category of finitely generated Lie conformal algebras over Algebraic Field
    """
    class ParentMethods:
        def ngens(self):
            r"""
            The number of generators of this Lie conformal algebra.

            EXAMPLES::

                sage: Vir = lie_conformal_algebras.Virasoro(QQ)
                sage: Vir.ngens()
                2

                sage: V = lie_conformal_algebras.Affine(QQ, 'A2')
                sage: V.ngens()
                9
            """
            return len(self.gens())

        def gen(self,i):
            r"""
            The ``i``-th generator of this Lie conformal algebra.

            EXAMPLES::

                sage: V = lie_conformal_algebras.Affine(QQ, 'A1')
                sage: V.gens()
                (B[alpha[1]], B[alphacheck[1]], B[-alpha[1]], B['K'])
                sage: V.gen(0)
                B[alpha[1]]
                sage: V.1
                B[alphacheck[1]]
            """
            return self.gens()[i]

        def some_elements(self):
            """
            Some elements of this Lie conformal algebra.

            This method returns a list with elements containing at
            least the generators.

            EXAMPLES::

                sage: V = lie_conformal_algebras.Affine(QQ, 'A1', names=('e', 'h', 'f'))
                sage: V.some_elements()
                [e, h, f, K, Th + 4*T^(2)e, 4*T^(2)h, Te + 4*T^(2)e, Te + 4*T^(2)h]
            """
            S = list(self.gens())
            from sage.misc.misc import some_tuples
            for x,y in some_tuples(S, 2, 0, max_samples=self.ngens()):
                S.append(x.T() + 2*y.T(2))
            return S

    class Super(SuperModulesCategory):
        """
        The category of super finitely generated Lie conformal algebras.

        EXAMPLES::

            sage: LieConformalAlgebras(AA).FinitelyGenerated().Super()
            Category of super finitely generated Lie conformal algebras over Algebraic Real Field
        """
        def extra_super_categories(self):
            """
            The extra super categories of this category

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).FinitelyGenerated().Super().super_categories()
                [Category of super Lie conformal algebras over Rational Field,
                 Category of finitely generated Lie conformal algebras over Rational Field]
            """
            return [self.base_category(),]

    class Graded(GradedModulesCategory):
        """
        The category of H-graded finitely generated Lie conformal algebras.

        EXAMPLES::

            sage: LieConformalAlgebras(QQbar).FinitelyGenerated().Graded()
            Category of H-graded finitely generated Lie conformal algebras over Algebraic Field
        """
        def _repr_object_names(self):
            """
            The names of the objects of this category

            EXAMPLES::

                sage: LieConformalAlgebras(QQbar).FinitelyGenerated().Graded()
                Category of H-graded finitely generated Lie conformal algebras over Algebraic Field
            """
            return "H-graded {}".format(self.base_category().\
                                        _repr_object_names())

        class Super(SuperModulesCategory):
            """
            The subcategory of super H-graded finitely generated Lie
            conformal algebras.

            EXAMPLES::

                sage: LieConformalAlgebras(QQbar).FinitelyGenerated().Graded()
                Category of H-graded finitely generated Lie conformal algebras over Algebraic Field
            """
            def extra_super_categories(self):
                """
                The extra super categories of this category.

                EXAMPLES::

                    sage: C = LieConformalAlgebras(QQ).FinitelyGenerated().Graded()
                    sage: C.super_categories()
                    [Category of H-graded Lie conformal algebras over Rational Field,
                     Category of finitely generated Lie conformal algebras over Rational Field]
                """
                return [self.base_category(),]
