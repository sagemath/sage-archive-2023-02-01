r"""
Finite dimensional modules with basis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring

class FiniteDimensionalModulesWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The category of finite dimensional modules with a distinguished basis

    EXAMPLES::

      sage: C = FiniteDimensionalModulesWithBasis(ZZ); C
      Category of finite dimensional modules with basis over Integer Ring
      sage: sorted(C.super_categories(), key=str)
      [Category of finite dimensional modules over Integer Ring,
       Category of modules with basis over Integer Ring]
      sage: C is Modules(ZZ).WithBasis().FiniteDimensional()
      True

    TESTS::

        sage: TestSuite(C).run()
    """

    class ParentMethods:
        pass

    class ElementMethods:
        pass

    class MorphismMethods:
        def matrix(self, base_ring=None):
            """
            Returns the matrix of self in the distinguished basis
            of the domain and codomain.

            EXAMPLES::

                sage: category = FiniteDimensionalModulesWithBasis(ZZ)
                sage: X = CombinatorialFreeModule(ZZ, [1,2], category = category); X.rename("X"); x = X.basis()
                sage: Y = CombinatorialFreeModule(ZZ, [3,4], category = category); Y.rename("Y"); y = Y.basis()
                sage: phi = X.module_morphism(on_basis = {1: y[3] + 3*y[4], 2: 2*y[3] + 5*y[4]}.__getitem__,
                ...                           codomain = Y, category = category)
                sage: phi.matrix()
                [1 2]
                [3 5]
                sage: phi.matrix().parent()
                Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
                sage: phi.matrix().is_mutable()
                False

                sage: category = FiniteDimensionalModulesWithBasis(QQ)
                sage: X = CombinatorialFreeModule(QQ, [1, 3, 7], category = category); X.rename("X"); x = X.basis()
                sage: Y = CombinatorialFreeModule(QQ, [2, 4], category = category); Y.rename("Y"); y = Y.basis()
                sage: Hom(X,Y).zero().matrix()
                [0 0 0]
                [0 0 0]
            """
            from sage.matrix.constructor import matrix
            if base_ring is None:
                base_ring = self.domain().base_ring()

            on_basis = self.on_basis()
            basis_keys = self.domain().basis().keys()
            m = matrix(base_ring,
                       [on_basis(x).to_vector() for x in basis_keys]).transpose()
            m.set_immutable()
            return m

        def _from_matrix(self, hom, m):
            """
            Construct a morphism in the homset ``hom`` from the matrix ``m``

            EXAMPLES::

                sage: category = FiniteDimensionalModulesWithBasis(ZZ)
                sage: X = CombinatorialFreeModule(ZZ, [1,2], category = category); X.rename("X"); x = X.basis()
                sage: Y = CombinatorialFreeModule(ZZ, [3,4], category = category); Y.rename("Y"); y = Y.basis()
                sage: H = Hom(X,Y)

            This is a static method which can be called from any
            element of ``hom``::

                sage: C = H.zero()
                sage: phi = C._from_matrix(H, matrix([[1,2],[3,5]]))
                sage: phi.parent()
                Set of Morphisms from X to Y in Category of finite dimensional modules with basis over Integer Ring
                sage: phi(x[1])
                B[3] + 3*B[4]
                sage: phi(x[2])
                2*B[3] + 5*B[4]

            .. TODO::

                Design a better API. We would want to do something
                like ``hom(matrix=)`` which would call
                ``hom._from_matrix`` or
                ``domain.morphism_from_matrix``. The thing is to
                handle proper inheritance since this should be
                available to homsets in any subcategory of
                ModulesWithBasis.FiniteDimensional. See also
                discussion on trac:`10668`.
            """
            X = hom.domain()
            Y = hom.codomain()
            import sage.combinat.ranker
            rankX  = sage.combinat.ranker.rank_from_list(tuple(X.basis().keys()))
            d = dict( (xt, Y.from_vector(m.column(rankX(xt)))) for xt in X.basis().keys() )
            return hom(on_basis = d.__getitem__)

        def __invert__(self):
            """
            Returns the inverse morphism of ``self``, or raise an
            error if ``self`` is not invertible (should this
            return None instead???). This is achieved by inverting
            the ``self.matrix()``.

            EXAMPLES::

                sage: category = FiniteDimensionalModulesWithBasis(ZZ)
                sage: X = CombinatorialFreeModule(ZZ, [1,2], category = category); X.rename("X"); x = X.basis()
                sage: Y = CombinatorialFreeModule(ZZ, [3,4], category = category); Y.rename("Y"); y = Y.basis()
                sage: phi = X.module_morphism(on_basis = {1: y[3] + 3*y[4], 2: 2*y[3] + 5*y[4]}.__getitem__,
                ...                           codomain = Y, category = category)
                sage: psi = ~phi
                sage: psi
                Generic morphism:
                  From: Y
                  To:   X
                sage: psi.parent()
                Set of Morphisms from Y to X in Category of finite dimensional modules with basis over Integer Ring
                sage: psi(y[3])
                -5*B[1] + 3*B[2]
                sage: psi(y[4])
                2*B[1] - B[2]
                sage: psi.matrix()
                [-5  2]
                [ 3 -1]
                sage: psi(phi(x[1])), psi(phi(x[2]))
                (B[1], B[2])
                sage: phi(psi(y[3])), phi(psi(y[4]))
                (B[3], B[4])

            We check that this function complains if the morphism is not invertible::

                sage: phi = X.module_morphism(on_basis = {1: y[3] + y[4], 2: y[3] + y[4]}.__getitem__,
                ...                           codomain = Y, category = category)
                sage: ~phi
                Traceback (most recent call last):
                ...
                RuntimeError: morphism is not invertible

                sage: phi = X.module_morphism(on_basis = {1: y[3] + y[4], 2: y[3] + 5*y[4]}.__getitem__,
                ...                           codomain = Y, category = category)
                sage: ~phi
                Traceback (most recent call last):
                ...
                RuntimeError: morphism is not invertible
            """
            from sage.categories.homset import Hom
            mat = self.matrix()
            try:
                inv_mat = mat.parent()(~mat)
            except (ZeroDivisionError, TypeError):
                raise RuntimeError, "morphism is not invertible"
            return self._from_matrix(Hom(self.codomain(), self.domain(), category = self.category_for()),
                                     inv_mat)

