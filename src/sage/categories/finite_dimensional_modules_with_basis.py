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
        def matrix(self, base_ring=None, side="left"):
            r"""
            Return the matrix of this morphism in the distinguished
            bases of the domain and codomain.

            INPUT:

            - ``base_ring`` -- a ring (default: ``None``, meaning the
              base ring of the codomain)

            - ``side`` -- "left" or "right" (default: "left")

            If ``side`` is "left", this morphism is considered as
            acting on the left; i.e. each column of the matrix
            represents the image of an element of the basis of the
            domain.

            The order of the rows and columns matches with the order
            in which the bases are enumerated.

            .. SEEALSO:: :func:`_from_matrix`

            EXAMPLES::

                sage: X = CombinatorialFreeModule(ZZ, [1,2]); x = X.basis()
                sage: Y = CombinatorialFreeModule(ZZ, [3,4]); y = Y.basis()
                sage: phi = X.module_morphism(on_basis = {1: y[3] + 3*y[4], 2: 2*y[3] + 5*y[4]}.__getitem__,
                ...                           codomain = Y)
                sage: phi.matrix()
                [1 2]
                [3 5]
                sage: phi.matrix(side="right")
                [1 3]
                [2 5]

                sage: phi.matrix().parent()
                Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
                sage: phi.matrix(QQ).parent()
                Full MatrixSpace of 2 by 2 dense matrices over Rational Field

            The resulting matrix is immutable::

                sage: phi.matrix().is_mutable()
                False

            The zero morphism has a zero matrix::

                sage: Hom(X,Y).zero().matrix()
                [0 0]
                [0 0]

            .. TODO::

                Add support for morphisms where the codomain has a
                different base ring than the domain::

                    sage: Y = CombinatorialFreeModule(QQ, [3,4]); y = Y.basis()
                    sage: phi = X.module_morphism(on_basis = {1: y[3] + 3*y[4], 2: 2*y[3] + 5/2*y[4]}.__getitem__,
                    ...                           codomain = Y)
                    sage: phi.matrix().parent()          # todo: not implemented
                    Full MatrixSpace of 2 by 2 dense matrices over Rational Field

                This currently does not work because, in this case,
                the morphism is just in the category of commutative
                additive groups (i.e. the intersection of the
                categories of modules over `\ZZ` and over `\QQ`)::

                    sage: phi.parent().homset_category()
                    Category of commutative additive semigroups
                    sage: phi.parent().homset_category() # todo: not implemented
                    Category of finite dimensional modules with basis over Integer Ring
            """
            from sage.matrix.constructor import matrix
            if base_ring is None:
                base_ring = self.codomain().base_ring()

            on_basis = self.on_basis()
            basis_keys = self.domain().basis().keys()
            m = matrix(base_ring,
                       [on_basis(x).to_vector() for x in basis_keys])
            if side == "left":
                m = m.transpose()
            m.set_immutable()
            return m

        def __invert__(self):
            """
            Return the inverse morphism of ``self``

            This is achieved by inverting the ``self.matrix()``. An
            error is raised if ``self`` is not invertible.

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
            return self.codomain().module_morphism(
                matrix=inv_mat,
                codomain=self.domain(), category=self.category_for())

