r"""
Graded and stratified Lie algebras

AUTHORS:

- Eero Hakavuori (2018-08-16): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Eero Hakavuori <eero.hakavuori@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.category import Category
from sage.categories.category_types import Category_over_base_ring
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory
from sage.categories.graded_modules import GradedModulesCategory


class GradedLieAlgebras(GradedModulesCategory):
    r"""
    Category of graded Lie algebras.

    TESTS::

        sage: C = LieAlgebras(QQ).Graded()
        sage: TestSuite(C).run()

    """
    pass


class FiniteDimensionalGradedLieAlgebrasWithBasis(GradedModulesCategory):
    r"""
    Category of finite dimensional graded Lie algebras with a basis.
    
    A grading of a Lie algebra `\mathfrak{g}` is a direct sum decomposition
    `\mathfrak{g} = \bigoplus_{i} V_i` such that `[V_i,V_j] \subset V_{i+j}`.

    TESTS::

        sage: C = LieAlgebras(QQ).FiniteDimensional().WithBasis().Graded()
        sage: TestSuite(C).run()

    """

    class ParentMethods:

        def homogeneous_component_as_submodule(self, d):
            """
            Return the ``d``-th homogeneous component of ``self`` as a submodule.

            EXAMPLES::

                sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra
                sage: L = NilpotentLieAlgebra(QQ, 2, step=2)
                sage: L.homogeneous_component_as_submodule(2)
                Sparse vector space of degree 3 and dimension 1 over Rational Field
                Basis matrix:
                [0 0 1]
            """
            B = self.homogeneous_component_basis(d)
            return self.module().submodule([X.to_vector() for X in B])

        def _test_grading(self, **options):
            r"""
            Tests that the Lie bracket respects the grading.

            INPUT:

            - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

            EXAMPLES::

                sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra
                sage: L.<X,Y,Z> = NilpotentLieAlgebra(ZZ, 2, step=2)
                sage: L._test_grading()
                sage: L._basis_degrees[Y] = 2
                sage: L._test_grading()
                Traceback (most recent call last):
                ...
                AssertionError: Lie bracket [X, Y] is not in the homogeneous component of degree 3

            See the documentation for :class:`TestSuite` for more information.
            """
            tester = self._tester(**options)

            for X_ind in self.indices():
                X = self.basis()[X_ind]
                i = self.degree_on_basis(X_ind)
                for Y_ind in self.indices():
                    Y = self.basis()[Y_ind]
                    j = self.degree_on_basis(Y_ind)
                    Z = self.bracket(X, Y).to_vector()
                    tester.assertTrue(
                        Z in self.homogeneous_component_as_submodule(i + j),
                        msg="Lie bracket [%s, %s] is not in the "
                        "homogeneous component of degree %d" % (X, Y, i + j))


class StratifiedLieAlgebrasCategory(RegressiveCovariantConstructionCategory,
                                    Category_over_base_ring):

    def __init__(self, base_category):
        super(StratifiedLieAlgebrasCategory, self).__init__(base_category,
                                                            base_category.base_ring())

    _functor_category = "Stratified"

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: LieAlgebras(QQ).Stratified()  # indirect doctest
            Category of stratified Lie algebras over Rational Field
        """
        return "stratified {}".format(self.base_category()._repr_object_names())

    @classmethod
    def default_super_categories(cls, category, *args):
        r"""
        Return the default super categories of ``category.Stratified()``.

        Mathematical meaning: every stratified Lie algebra is also graded.

        INPUT:

        - ``cls`` -- the class ``StratifiedLieAlgebrasCategory``
        - ``category`` -- a category

        OUTPUT: a (join) category

        In practice, this returns ``category.Stratified()``, joined
        together with the result of the method
        :meth:`RegressiveCovariantConstructionCategory.default_super_categories() <sage.categories.covariant_functorial_construction.RegressiveCovariantConstructionCategory.default_super_categories>`
        (that is the join of ``category.Graded()`` and ``cat`` for
        each ``cat`` in the super categories of ``category``).

        EXAMPLES:

        Consider ``category=LieAlgebras(QQ)``. Then, a stratification of a Lie 
        algebra `L` is also a grading::

            sage: LieAlgebras(QQ).Stratified().super_categories()
            [Category of graded Lie algebras over Rational Field]

        This resulted from the following call::

            sage: sage.categories.graded_lie_algebras.StratifiedLieAlgebrasCategory.default_super_categories(LieAlgebras(QQ))
            Category of graded Lie algebras over Rational Field
        """
        P = super(StratifiedLieAlgebrasCategory, cls)
        cat = P.default_super_categories(category, *args)
        return Category.join([category.Graded(), cat])


class StratifiedLieAlgebras(StratifiedLieAlgebrasCategory):
    r"""
    Category of stratified Lie algebras.

    A stratified Lie algebra is a graded Lie algebra that is generated as a Lie
    algebra by its homogeneous component of degree 1.

    TESTS::

        sage: C = LieAlgebras(QQ).Stratified()
        sage: TestSuite(C).run()

    """
    pass


class FiniteDimensionalStratifiedLieAlgebrasWithBasis(StratifiedLieAlgebrasCategory):
    r"""
    Category of finite dimensional stratified Lie algebras with a basis.

    A stratified Lie algebra is a graded Lie algebra that is generated as a Lie
    algebra by its homogeneous component of degree 1.

    TESTS::

        sage: C = LieAlgebras(QQ).FiniteDimensional().WithBasis().Stratified()
        sage: TestSuite(C).run()

    """

    def extra_super_categories(self):
        r"""
        Implements the fact that a finite dimensional stratified Lie algebra
        is necessarily nilpotent.
        """
        R = self.base_ring()
        from sage.categories.lie_algebras import LieAlgebras
        return [LieAlgebras(R).FiniteDimensional().WithBasis().Nilpotent()]

    class ParentMethods:

        def degree_on_basis(self, m):
            r"""
            Returns the degree of the basis element ``m``

            INPUT:

            - ``m`` -- an element in ``self.indices()`` or in ``self.basis()``

            If the degrees of the basis elements are not defined, they will
            be computed. By assumption the stratification
            `V_1 \oplus \dots \oplus V_s` of ``self`` is such that each
            component `V_k` is spanned by some subset of the basis.

            The degree of a basis element `X` is therefore the largest index
            `k`such that `X \in V_k\oplus\dots\oplus V_s`. The space 
            `V_k\oplus\dots\oplus V_s` is by assumption the `k`th term of the 
            lower central series. 

            EXAMPLES::
 
                sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra
                sage: L.<X,Y,Z> = NilpotentLieAlgebra(QQ, 2, step=2)
                sage: L.degree_on_basis(X)
                1
                sage: L.degree_on_basis(Y)
                1
                sage: L[X, Y]
                Z
                sage: L.degree_on_basis(Z)
                2
            """
            if not hasattr(self, '_basis_degrees'):
                lcs = self.lower_central_series(submodule=True)
                self._basis_degrees = {}

                for k in reversed(range(len(lcs) - 1)):
                    for X in self.basis():
                        if X in self._basis_degrees:
                            continue
                        if X.to_vector() in lcs[k]:
                            self._basis_degrees[X] = k + 1

            if m in self.indices():
                m = self.basis()[m]
            return self._basis_degrees[m]

        def _test_generated_by_degree_one(self, **options):
            r"""
            Tests that the Lie algebra is generated by the homogeneous component
            of degree one.

            INPUT:

            - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

            EXAMPLES::

                sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra
                sage: L = NilpotentLieAlgebra(QQ, 2, step=3)
                sage: L._test_generated_by_degree_one()
                sage: L.inject_variables()
                Defining X_1, X_2, X_12, X_112, X_122
                sage: L._basis_degrees[X_2] = 2
                sage: L._test_generated_by_degree_one()
                Traceback (most recent call last):
                ...
                AssertionError: [X_1] does not generate Nilpotent Lie algebra on
                5 generators (X_1, X_2, X_12, X_112, X_122) over Rational Field

            See the documentation for :class:`TestSuite` for more information.
            """
            tester = self._tester(**options)

            V1 = self.homogeneous_component_as_submodule(1)
            B1 = V1.basis()
            m = self.module()

            V = V1
            d = 0
            while V.dimension() > d:
                B = V.basis()
                d = V.dimension()
                V = m.submodule(B + [self.bracket(X, Y).to_vector()
                                 for X in B1 for Y in B])

            tester.assertEqual(V, m,
                msg="%s does not generate %s" % ([self(X) for X in B1], self))
