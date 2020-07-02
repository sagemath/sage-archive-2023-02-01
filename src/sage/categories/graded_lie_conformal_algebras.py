r"""
Graded Lie Conformal Algebras

A Lie conformal algebra `L` is called *H-Graded* [DSK2006]_ if there exists
a decomposition `L = \oplus L_n` such that the
`\lambda`-bracket becomes graded of degree `-1`, that is:

.. MATH::

    a_{(n)} b \in L_{p + q -n -1} \qquad
    a \in L_p, \: b \in L_q, \: n \geq 0.

In particular this implies that the action of `T` increases
degree by `1`.

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

from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.super_modules import SuperModulesCategory
from sage.misc.abstract_method import abstract_method

class GradedLieConformalAlgebras(GradedModulesCategory):
    """
    The subcategory of H-graded Lie conformal algebras.

    EXAMPLES::

        sage: LieConformalAlgebras(QQbar).Graded()
        Category of H-graded Lie conformal algebras over Algebraic Field
    """

    def _repr_object_names(self):
        """
        The names of the objects of this category

        EXAMPLES::

            sage: LieConformalAlgebras(QQbar).Super().Graded()
            Category of super H-graded Lie conformal algebras over Algebraic Field
            sage: LieConformalAlgebras(ZZ).Graded().FinitelyGenerated().WithBasis()
            Category of finitely generated H-graded Lie conformal algebras with basis over Integer Ring
        """
        return "H-graded {}".format(self.base_category().\
                                    _repr_object_names())


    class SubcategoryMethods:

        def Super(self, base_ring=None):
            """
            The subcategory of H-graded super Lie conformal algebras

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).Super().Graded()
                Category of super H-graded Lie conformal algebras over Rational Field
            """
            assert base_ring is None or base_ring is self.base_ring()
            if isinstance(self,CategoryWithAxiom_over_base_ring):
                axioms_whitelist = frozenset(["WithBasis",
                                  "FinitelyGeneratedAsLieConformalAlgebra"])
                axioms = axioms_whitelist.intersection(self.axioms())
                return self.base_category().Super()._with_axioms(axioms)
            return SuperModulesCategory.category_of(self)

    class WithBasis(CategoryWithAxiom_over_base_ring):
        """
        The subcategory of H-graded Lie conformal algebras with
        basis.

        EXAMPLES::

            sage: LieConformalAlgebras(ZZ).Graded().WithBasis()
            Category of H-graded Lie conformal algebras with basis over Integer Ring
        """

        class FinitelyGeneratedAsLieConformalAlgebra(
                                        CategoryWithAxiom_over_base_ring):
            """
            The subcategory of finitely generated H-graded
            Lie conformal algebras with basis.

            EXAMPLES::

                sage: LieConformalAlgebras(ZZ).Graded().FinitelyGenerated().WithBasis()
                Category of finitely generated H-graded Lie conformal algebras with basis over Integer Ring
            """
            pass

    class FinitelyGeneratedAsLieConformalAlgebra(
                                        CategoryWithAxiom_over_base_ring):
        """
        The subcategory of finitely generated H-graded Lie
        conformal algebras.

        EXAMPLES::

            sage: C = LieConformalAlgebras(ZZ).Graded().FinitelyGenerated(); C
            Category of finitely generated H-graded Lie conformal algebras over Integer Ring
            sage: C is LieConformalAlgebras(ZZ).FinitelyGenerated().Graded()
            True
        """
        pass

    class Super(SuperModulesCategory):
        """
        The subcategory of H-graded super Lie conformal algebras.

        EXAMPLES::

            sage: C = LieConformalAlgebras(QQbar).Graded().Super(); C
            Category of super H-graded Lie conformal algebras over Algebraic Field
            sage: C is LieConformalAlgebras(QQbar).Super().Graded()
            True
        """
        def extra_super_categories(self):
            """
            The extra super categories of this category.

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).FinitelyGenerated().Graded().Super().super_categories()
                [Category of finitely generated super Lie conformal algebras over Rational Field,
                 Category of finitely generated H-graded Lie conformal algebras over Rational Field,
                 Category of super H-graded Lie conformal algebras over Rational Field]
            """
            return [self.base_category(),]


    class ElementMethods:
        @abstract_method
        def degree(self):
            """
            The degree of this element.

            EXAMPLES::

                sage: Vir = lie_conformal_algebras.Virasoro(QQ); L = Vir.0
                sage: L.degree()
                2
                sage: L.T(3).degree()
                5
            """
            raise NotImplementedError("Not implemented")
