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

import operator
from sage.misc.cachefunc import cached_method
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

        @cached_method
        def vectors_parent(self):
            """
            Returns the parent of the vectors created with
            ``x.to_vector()``.

            EXAMPES::

                sage: F = Modules(QQ).WithBasis().example()
                sage: F.vectors_parent()
                
            """
            return self.zero().to_vector().parent()

        def span(self, gens, check = True, already_echelonized = False):
            return self.vectors_parent().span(
                [g.to_vector() for g in gens],
                check=check,
                already_echelonized = already_echelonized)

        def annihilator(self, S, action=operator.mul, side='right', category=None):
            r"""
            INPUT:
             - ``S`` -- a finite set of objects
             - ``side`` -- 'left' or 'right' (default: 'right')

             - ``action`` -- a function (default: `operator.mul`)

            Assumptions: ``action`` takes elements of ``self`` as
            first argument, and elements of ``S`` as second
            argument, and is linear on its first argument
            (typically it is bilinear).

            Returns the supspace of the elements `x` of ``self``
            such that `action(x,s) = 0` for all `s\in S`.

            If ``side`` is 'left', then the order of the arguments
            of ``action`` is reversed.

            TODO: double check the convention for ``left/right``.

            .. seealso:: :meth:`annihilator_basis` for lots of examples.

            EXAMPLES::

                sage: F = FiniteDimensionalAlgebrasWithBasis(QQ).example(); F
                An example of a finite dimensional algebra with basis: the path algebra of the Kronecker quiver (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: x,y,a,b = F.basis()
                sage: A = F.annihilator([a + 3*b + 2*y]); A
                Vector space of degree 4 and dimension 1 over Rational Field
                Basis matrix:
                [   1    0 -1/2 -3/2]

            Note: currently `A` forgot everything about `F`; its
            ambient space is a :class:`FreeModule` of same dimension
            as `F`. But one can, for example, compute inclusions::

                sage: A   = F.annihilator([]);    A  .rename("A")
                sage: Ax  = F.annihilator([x]);   Ax .rename("Ax")
                sage: Ay  = F.annihilator([y]);   Ay .rename("Ay")
                sage: Axy = F.annihilator([x,y]); Axy.rename("Axy")
                sage: P = Poset(([A, Ax, Ay, Axy], attrcall("is_subspace")))
                sage: P.cover_relations()
                [[Axy, Ay], [Axy, Ax], [Ay, A], [Ax, A]]

            """
            return self.span(self.annihilator_basis(S, action, side), already_echelonized=True)

        def annihilator_basis(self, S, action=operator.mul, side='right'):
            """

            EXAMPLES:

            By default, the action is the standard `*`
            operation. So our first example is about an algebra::

                sage: F = FiniteDimensionalAlgebrasWithBasis(QQ).example(); F
                An example of a finite dimensional algebra with basis: the path algebra of the Kronecker quiver (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: x,y,a,b = F.basis()

            In this algebra, multiplication on the right by `x`
            kills all basis elements but `x`::

                sage: x*x, y*x, a*x, b*x
                (x, 0, 0, 0)

            So the annihilator is the subspace spanned by
            `y`,`a`,and `b`::

                sage: F.annihilator_basis([x])
                [y, a, b]

            The same holds for `a` and `b`::

                sage: x*a, y*a, a*a, b*a
                (a, 0, 0, 0)
                sage: F.annihilator_basis([a])
                [y, a, b]

            On the other hand, `y` kills only `x`::

                sage: F.annihilator_basis([y])
                [x]

            Here is a non trivial annihilator::

                sage: F.annihilator_basis([a + 3*b + 2*y])
                [-1/2*a - 3/2*b + x]

            Let's check it::

                sage: (-1/2*a - 3/2*b + x) * (a + 3*b + 2*y)
                0

            Doing the same calculations on the left exchanges the
            roles of `x` and `y`::

                sage: F.annihilator_basis([y], side="left")
                [x, a, b]
                sage: F.annihilator_basis([a], side="left")
                [x, a, b]
                sage: F.annihilator_basis([b], side="left")
                [x, a, b]
                sage: F.annihilator_basis([x], side="left")
                [y]
                sage: F.annihilator_basis([a+3*b+2*x], side="left")
                [-1/2*a - 3/2*b + y]

            By specifying an inner product, this method can be
            used to compute the orthogonal of a subspace (TODO:
            add an example).

            By specifying instead the standard Lie bracket as
            action, one can compute the commutator of a subspace
            of `F`::

                sage: F.annihilator_basis([a+b], action = F.bracket)
                [x + y, a, b]

            In particular one can computer the center of the
            algebra. In our example, it is reduced to the
            identity::

                sage: F.annihilator_basis(F.algebra_generators(), action = F.bracket)
                [x + y]

            .. seealso:: :meth:`FiniteAlgebrasWithBasis.ParentMethods.center_basis`.
            """
            # TODO: optimize this!
            from sage.matrix.constructor import matrix
            if side == 'right':
                action_left = action
                action = lambda b,s: action_left(s, b)

            mat = matrix(self.base_ring(), self.dimension(), 0)
            for s in S:
                mat = mat.augment(matrix(self.base_ring(),
                                         [action(s,b).to_vector() for b in self.basis()]))
            return map(self.from_vector, mat.left_kernel().basis())

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

