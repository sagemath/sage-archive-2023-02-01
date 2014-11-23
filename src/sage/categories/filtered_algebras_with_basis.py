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

        # Maps

        def to_graded_conversion(self):
            r"""
            Return the canonical `R`-module isomorphism
            `A \to \operatorname{gr} A` induced by the basis of `A`.

            This is an isomorphism of `R`-modules, not of algebras. See
            the class documentation :class:`AssociatedGradedAlgebra`.

            .. SEEALSO::

                :meth:`from_graded_conversion`

            EXAMPLES::

                sage: A = Algebras(QQ).WithBasis().Filtered().example()
                sage: p = A.an_element() + A.algebra_generators()['x'] + 2; p
                U['x']^2*U['y']^2*U['z']^3 + U['x'] + 2
                sage: q = A.to_graded_conversion()(p); q
                bar(U['x']^2*U['y']^2*U['z']^3) + bar(U['x']) + 2*bar(1)
                sage: q.parent() is A.graded_algebra()
                True
            """
            base_one = self.base_ring().one()
            return self.module_morphism(diagonal=lambda x: base_one,
                                        codomain=self.graded_algebra())

        def from_graded_conversion(self):
            r"""
            Return the inverse of the canonical `R`-module isomorphism
            `A \to \operatorname{gr} A` induced by the basis of `A`
            (that is a map `\operatorname{gr} A \to A`).

            This is an isomorphism of `R`-modules, not of algebras. See
            the class documentation :class:`AssociatedGradedAlgebra`.

            .. SEEALSO::

                :meth:`to_graded_conversion`

            EXAMPLES::

                sage: A = Algebras(QQ).WithBasis().Filtered().example()
                sage: p = A.an_element() + A.algebra_generators()['x'] + 2; p
                U['x']^2*U['y']^2*U['z']^3 + U['x'] + 2
                sage: q = A.to_graded_conversion()(p)
                sage: A.from_graded_conversion()(q) == p
                True
                sage: q.parent() is A.graded_algebra()
                True
            """
            base_one = self.base_ring().one()
            return self.graded_algebra().module_morphism(diagonal=lambda x: base_one,
                                                         codomain=self)

        def projection(self, i):
            r"""
            Return the `i`-th projection `p_i : F_i \to G_i` (in the
            notations of the class documentation
            :class:`AssociatedGradedAlgebra`).

            This method actually does not return the map `p_i` itself,
            but an extension of `p_i` to the whole `R`-module `A`.
            This extension is the composition of the `R`-module
            isomorphism `A \to \operatorname{gr} A` with the canonical
            projection of the graded `R`-module `\operatorname{gr} A`
            onto its `i`-th graded component `G_i`. The codomain of
            this map is `\operatorname{gr} A`, although its actual
            image is `G_i`.

            EXAMPLES::

                sage: A = Algebras(QQ).WithBasis().Filtered().example()
                sage: p = A.an_element() + A.algebra_generators()['x'] + 2; p
                U['x']^2*U['y']^2*U['z']^3 + U['x'] + 2
                sage: q = A.projection(7)(p); q
                bar(U['x']^2*U['y']^2*U['z']^3)
                sage: q.parent() is A.graded_algebra()
                True
                sage: A.projection(8)(p)
                0
            """
            base_zero = self.base_ring().zero()
            base_one = self.base_ring().one()
            grA = self.graded_algebra()
            proj = lambda x: (base_one if grA.degree_on_basis(x) == i
                              else base_zero)
            return self.module_morphism(diagonal=proj, codomain=grA)

    class ElementMethods:
        pass

