r"""
Finite Dimensional Lie Algebras With Basis

AUTHORS:

- Travis Scrimshaw (07-15-2013): Initial implementation
"""

#*****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.lie_algebras import LieAlgebras
from sage.categories.subobjects import SubobjectsCategory
from sage.algebras.free_algebra import FreeAlgebra
from sage.sets.family import Family
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector

class FiniteDimensionalLieAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    Category of finite dimensional Lie algebras with a basis.

    .. TODO::

        Many of these tests should use non-abelian Lie algebras and need to
        be added after :trac:`16820`.
    """
    _base_category_class_and_axiom = [LieAlgebras.FiniteDimensional, "WithBasis"]

    def example(self, n=3):
        """
        Return an example of a finite dimensional Lie algebra with basis as per
        :meth:`Category.example <sage.categories.category.Category.example>`.

        EXAMPLES::

            sage: C = LieAlgebras(QQ).FiniteDimensional().WithBasis()
            sage: C.example()
            An example of a finite dimensional Lie algebra with basis:
             the 3-dimensional abelian Lie algebra over Rational Field

        Other dimensions can be specified as an optional argument::

            sage: C.example(5)
            An example of a finite dimensional Lie algebra with basis:
             the 5-dimensional abelian Lie algebra over Rational Field
        """
        from sage.categories.examples.finite_dimensional_lie_algebras_with_basis import Example
        return Example(self.base_ring(), n)

    class ParentMethods:
        @cached_method
        def _construct_UEA(self):
            """
            Construct the universal enveloping algebra of ``self``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: UEA = L._construct_UEA(); UEA
                Noncommutative Multivariate Polynomial Ring in b0, b1, b2
                 over Rational Field, nc-relations: {}
                sage: UEA.relations(add_commutative=True)
                {b1*b0: b0*b1, b2*b0: b0*b2, b2*b1: b1*b2}

            ::

                sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}, ('y','z'):{'x':1}, ('z','x'):{'y':1}})
                sage: UEA = L._construct_UEA(); UEA
                Noncommutative Multivariate Polynomial Ring in x, y, z over Rational Field,
                 nc-relations: {...}
                sage: sorted(UEA.relations().items(), key=str)
                [(y*x, x*y - z), (z*x, x*z + y), (z*y, y*z - x)]
            """
            # Create the UEA relations
            # We need to get names for the basis elements, not just the generators
            I = self._basis_ordering
            try:
                names = [str(x) for x in I]
                def names_map(x): return x
                F = FreeAlgebra(self.base_ring(), names)
            except ValueError:
                names = ['b{}'.format(i) for i in range(self.dimension())]
                self._UEA_names_map = {g: names[i] for i,g in enumerate(I)}
                names_map = self._UEA_names_map.__getitem__
                F = FreeAlgebra(self.base_ring(), names)
            # ``F`` is the free algebra over the basis of ``self``. The
            # universal enveloping algebra of ``self`` will be constructed
            # as a quotient of ``F``.
            d = F.gens_dict()
            rels = {}
            S = self.structure_coefficients(True)
            # Construct the map from indices to names of the UEA
            def get_var(g): return d[names_map(g)]
            # The function ``get_var`` sends an element of the basis of
            # ``self`` to the corresponding element of ``F``.
            for k in S.keys():
                g0 = get_var(k[0])
                g1 = get_var(k[1])
                if g0 < g1:
                    rels[g1*g0] = g0*g1 - F.sum(val*get_var(g) for g, val in S[k])
                else:
                    rels[g0*g1] = g1*g0 + F.sum(val*get_var(g) for g, val in S[k])
            return F.g_algebra(rels)

        @lazy_attribute
        def _basis_ordering(self):
            """
            Return the indices of the basis of ``self`` as a tuple in
            a fixed order.

            Override this attribute to get a specific ordering of the basis.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L._basis_ordering
                (0, 1, 2)
            """
            return tuple(self.basis().keys())

        def _dense_free_module(self, R=None):
            """
            Return a dense free module associated to ``self`` over ``R``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L._dense_free_module()
                Vector space of dimension 3 over Rational Field
            """
            if R is None:
                R = self.base_ring()
            from sage.modules.free_module import FreeModule
            return FreeModule(R, self.dimension())

        module = _dense_free_module

        def from_vector(self, v):
            """
            Return the element of ``self`` corresponding to the
            vector ``v`` in ``self.module()``.

            Implement this if you implement :meth:`module`; see the
            documentation of
            :meth:`sage.categories.lie_algebras.LieAlgebras.module`
            for how this is to be done.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: u = L.from_vector(vector(QQ, (1, 0, 0))); u
                (1, 0, 0)
                sage: parent(u) is L
                True
            """
            B = self.basis()
            return self.sum(v[i] * B[k] for i,k in enumerate(self._basis_ordering)
                            if v[i] != 0)

        def killing_matrix(self, x, y):
            r"""
            Return the Killing matrix of ``x`` and ``y``, where ``x``
            and ``y`` are two elements of ``self``.

            The Killing matrix is defined as the matrix corresponding
            to the action of
            `\operatorname{ad}_x \circ \operatorname{ad}_y` in the
            basis of ``self``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: a,b,c = L.lie_algebra_generators()
                sage: L.killing_matrix(a, b)
                [0 0 0]
                [0 0 0]
                [0 0 0]

            ::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}})
                sage: L.killing_matrix(x, y)
                [ 0  0]
                [-1  0]
            """
            return x.adjoint_matrix() * y.adjoint_matrix()

        def killing_form(self, x, y):
            r"""
            Return the Killing form on ``x`` and ``y``, where ``x``
            and ``y`` are two elements of ``self``.

            The Killing form is defined as

            .. MATH::

                \langle x \mid y \rangle
                = \operatorname{tr}\left( \operatorname{ad}_x
                \circ \operatorname{ad}_y \right).

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: a,b,c = L.lie_algebra_generators()
                sage: L.killing_form(a, b)
                0
            """
            return self.killing_matrix(x, y).trace()

        @cached_method
        def killing_form_matrix(self):
            """
            Return the matrix of the Killing form of ``self``.

            The rows and the columns of this matrix are indexed by the
            elements of the basis of ``self`` (in the order provided by
            :meth:`basis`).

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.killing_form_matrix()
                [0 0 0]
                [0 0 0]
                [0 0 0]

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example(0)
                sage: m = L.killing_form_matrix(); m
                []
                sage: parent(m)
                Full MatrixSpace of 0 by 0 dense matrices over Rational Field
            """
            B = self.basis()
            m = matrix(self.base_ring(),
                       [[self.killing_form(x, y) for x in B] for y in B])
            m.set_immutable()
            return m

        @cached_method
        def structure_coefficients(self, include_zeros=False):
            """
            Return the structure coefficients of ``self``.

            INPUT:

            - ``include_zeros`` -- (default: ``False``) if ``True``, then
              include the `[x, y] = 0` pairs in the output

            OUTPUT:

            A dictionary whose keys are pairs of basis indices `(i, j)`
            with `i < j`, and whose values are the corresponding
            *elements* `[b_i, b_j]` in the Lie algebra.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.structure_coefficients()
                Finite family {}
                sage: L.structure_coefficients(True)
                Finite family {(0, 1): (0, 0, 0), (1, 2): (0, 0, 0), (0, 2): (0, 0, 0)}

            ::

                sage: G = SymmetricGroup(3)
                sage: S = GroupAlgebra(G, QQ)
                sage: L = LieAlgebra(associative=S)
                sage: L.structure_coefficients()
                Finite family {((1,3,2), (1,3)): (2,3) - (1,2),
                               ((1,2), (1,2,3)): -(2,3) + (1,3),
                               ((1,2,3), (1,3)): -(2,3) + (1,2),
                               ((2,3), (1,3,2)): -(1,2) + (1,3),
                               ((2,3), (1,3)): -(1,2,3) + (1,3,2),
                               ((2,3), (1,2)): (1,2,3) - (1,3,2),
                               ((2,3), (1,2,3)): (1,2) - (1,3),
                               ((1,2), (1,3,2)): (2,3) - (1,3),
                               ((1,2), (1,3)): (1,2,3) - (1,3,2)}
            """
            d = {}
            B = self.basis()
            K = list(B.keys())
            zero = self.zero()
            for i, x in enumerate(K):
                for y in K[i + 1:]:
                    bx = B[x]
                    by = B[y]
                    val = self.bracket(bx, by)
                    if not include_zeros and val == zero:
                        continue
                    if self._basis_key(x) > self._basis_key(y):
                        d[y,x] = -val
                    else:
                        d[x,y] = val
            return Family(d)

        def centralizer_basis(self, S):
            """
            Return a basis of the centralizer of ``S`` in ``self``.

            INPUT:

            - ``S`` -- a subalgebra of ``self`` or a list of elements that
              represent generators for a subalgebra

            .. SEEALSO::

                :meth:`centralizer`

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: a,b,c = L.lie_algebra_generators()
                sage: L.centralizer_basis([a + b, 2*a + c])
                [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

                sage: H = lie_algebras.Heisenberg(QQ, 2)
                sage: H.centralizer_basis(H)
                [z]


                sage: D = DescentAlgebra(QQ, 4).D()
                sage: L = LieAlgebra(associative=D)
                sage: L.centralizer_basis(L)
                [D{},
                 D{1} + D{1, 2} + D{2, 3} + D{3},
                 D{1, 2, 3} + D{1, 3} + D{2}]
                sage: D.center_basis()
                (D{},
                 D{1} + D{1, 2} + D{2, 3} + D{3},
                 D{1, 2, 3} + D{1, 3} + D{2})
            """
            #from sage.algebras.lie_algebras.subalgebra import LieSubalgebra
            #if isinstance(S, LieSubalgebra) or S is self:
            if S is self:
                from sage.matrix.special import identity_matrix
                m = identity_matrix(self.base_ring(), self.dimension())
            elif isinstance(S, (list, tuple)):
                m = matrix([v.to_vector() for v in self.echelon_form(S)])
            else:
                m = self.subalgebra(S).basis_matrix()

            S = self.structure_coefficients()
            sc = {}
            for k in S.keys():
                v = S[k].to_vector()
                sc[k] = v
                sc[k[1],k[0]] = -v
            X = self.basis().keys()
            d = len(X)
            c_mat = matrix(self.base_ring(),
                           [[sum(m[i,j] * sc[x,xp][k] for j,xp in enumerate(X)
                                 if (x, xp) in sc)
                             for x in X]
                            for i in range(d) for k in range(d)])
            C = c_mat.right_kernel().basis_matrix()
            return [self.from_vector(v) for v in C]

        def centralizer(self, S):
            """
            Return the centralizer of ``S`` in ``self``.

            INPUT:

            - ``S`` -- a subalgebra of ``self`` or a list of elements that
              represent generators for a subalgebra

            .. SEEALSO::

                :meth:`centralizer_basis`

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: a,b,c = L.lie_algebra_generators()
                sage: S = L.centralizer([a + b, 2*a + c]); S
                An example of a finite dimensional Lie algebra with basis:
                 the 3-dimensional abelian Lie algebra over Rational Field
                sage: S.basis_matrix()
                [1 0 0]
                [0 1 0]
                [0 0 1]
            """
            return self.subalgebra(self.centralizer_basis(S))

        def center(self):
            """
            Return the center of ``self``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: Z = L.center(); Z
                An example of a finite dimensional Lie algebra with basis: the
                 3-dimensional abelian Lie algebra over Rational Field
                sage: Z.basis_matrix()
                [1 0 0]
                [0 1 0]
                [0 0 1]
            """
            return self.centralizer(self)

        @cached_method
        def is_ideal(self, A):
            """
            Return if ``self`` is an ideal of ``A``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: a, b, c = L.lie_algebra_generators()
                sage: I = L.ideal([2*a - c, b + c])
                sage: I.is_ideal(L)
                True

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}})
                sage: L.is_ideal(L)
                True

                sage: F = LieAlgebra(QQ, 'F', representation='polynomial')
                sage: L.is_ideal(F)
                Traceback (most recent call last):
                ...
                NotImplementedError: A must be a finite dimensional Lie algebra
                 with basis
            """
            if A == self:
                return True
            if A not in LieAlgebras(self.base_ring()).FiniteDimensional().WithBasis():
                raise NotImplementedError("A must be a finite dimensional"
                                          " Lie algebra with basis")
            B = self.basis()
            AB = A.basis()
            try:
                b_mat = matrix(A.base_ring(), [A.bracket(b, ab).to_vector()
                                               for b in B for ab in AB])
            except (ValueError, TypeError):
                return False
            return b_mat.row_space().is_submodule(self.module())

        def product_space(self, L, submodule=False):
            r"""
            Return the product space ``[self, L]``.

            INPUT:

            - ``L`` -- a Lie subalgebra of ``self``
            - ``submodule`` -- (default: ``False``) if ``True``, then the
              result is forced to be a submodule of ``self``

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: a,b,c = L.lie_algebra_generators()
                sage: X = L.subalgebra([a, b+c])
                sage: L.product_space(X)
                An example of a finite dimensional Lie algebra with basis:
                 the 0-dimensional abelian Lie algebra over Rational Field
                 with basis matrix:
                []
                sage: Y = L.subalgebra([a, 2*b-c])
                sage: X.product_space(Y)
                An example of a finite dimensional Lie algebra with basis:
                 the 0-dimensional abelian Lie algebra over Rational
                 Field with basis matrix:
                []

            ::

                sage: H = lie_algebras.Heisenberg(ZZ, 4)
                sage: Hp = H.product_space(H, submodule=True).basis()
                sage: [H.from_vector(v) for v in Hp]
                [z]

            ::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}})
                sage: Lp = L.product_space(L) # todo: not implemented - #17416
                sage: Lp # todo: not implemented - #17416
                Subalgebra generated of Lie algebra on 2 generators (x, y) over Rational Field with basis:
                (x,)
                sage: Lp.product_space(L) # todo: not implemented - #17416
                Subalgebra generated of Lie algebra on 2 generators (x, y) over Rational Field with basis:
                (x,)
                sage: L.product_space(Lp) # todo: not implemented - #17416
                Subalgebra generated of Lie algebra on 2 generators (x, y) over Rational Field with basis:
                (x,)
                sage: Lp.product_space(Lp) # todo: not implemented - #17416
                Subalgebra generated of Lie algebra on 2 generators (x, y) over Rational Field with basis:
                ()
            """
            # Make sure we lift everything to the ambient space
            try:
                A = self._ambient
            except AttributeError:
                try:
                    A = L._ambient
                except AttributeError:
                    A = self

            B = self.basis()
            LB = L.basis()
            b_mat = matrix(A.base_ring(), [A.bracket(b, lb).to_vector()
                                           for b in B for lb in LB])
            if submodule is True or not (self.is_ideal(A) and L.is_ideal(A)):
                return b_mat.row_space()
            # We echelonize the matrix here
            # TODO: Do we want to?
            b_mat.echelonize()
            r = b_mat.rank()
            gens = [A.from_vector(row) for row in b_mat.rows()[:r]]
            return A.ideal(gens)

        @cached_method
        def derived_subalgebra(self):
            """
            Return the derived subalgebra of ``self``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.derived_subalgebra()
                An example of a finite dimensional Lie algebra with basis:
                 the 0-dimensional abelian Lie algebra over Rational Field
                 with basis matrix:
                []
            """
            return self.product_space(self)

        @cached_method
        def derived_series(self):
            r"""
            Return the derived series `(\mathfrak{g}^{(i)})_i` of ``self``
            where the rightmost
            `\mathfrak{g}^{(k)} = \mathfrak{g}^{(k+1)} = \cdots`.

            We define the derived series of a Lie algebra `\mathfrak{g}`
            recursively by `\mathfrak{g}^{(0)} := \mathfrak{g}` and

            .. MATH::

                \mathfrak{g}^{(k+1)} =
                [\mathfrak{g}^{(k)}, \mathfrak{g}^{(k)}]

            and recall that
            `\mathfrak{g}^{(k)} \supseteq \mathfrak{g}^{(k+1)}`.
            Alternatively we can express this as

            .. MATH::

                \mathfrak{g} \supseteq [\mathfrak{g}, \mathfrak{g}] \supseteq
                \bigl[ [\mathfrak{g}, \mathfrak{g}], [\mathfrak{g},
                \mathfrak{g}] \bigr] \supseteq
                \biggl[ \bigl[ [\mathfrak{g}, \mathfrak{g}], [\mathfrak{g},
                \mathfrak{g}] \bigr], \bigl[ [\mathfrak{g}, \mathfrak{g}],
                [\mathfrak{g}, \mathfrak{g}] \bigr] \biggr] \supseteq \cdots.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.derived_series()
                (An example of a finite dimensional Lie algebra with basis:
                    the 3-dimensional abelian Lie algebra over Rational Field,
                 An example of a finite dimensional Lie algebra with basis:
                    the 0-dimensional abelian Lie algebra over Rational Field
                    with basis matrix:
                    [])

            ::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}})
                sage: L.derived_series() # todo: not implemented - #17416
                (Lie algebra on 2 generators (x, y) over Rational Field,
                 Subalgebra generated of Lie algebra on 2 generators (x, y) over Rational Field with basis:
                (x,),
                 Subalgebra generated of Lie algebra on 2 generators (x, y) over Rational Field with basis:
                ())
            """
            L = [self]
            while L[-1].dimension() > 0:
                p = L[-1].derived_subalgebra()
                if L[-1].dimension() == p.dimension():
                    break
                L.append(p)
            return tuple(L)

        @cached_method
        def lower_central_series(self):
            r"""
            Return the lower central series `(\mathfrak{g}_{i})_i`
            of ``self`` where the rightmost
            `\mathfrak{g}_k = \mathfrak{g}_{k+1} = \cdots`.

            We define the lower central series of a Lie algebra `\mathfrak{g}`
            recursively by `\mathfrak{g}_0 := \mathfrak{g}` and

            .. MATH::

                \mathfrak{g}_{k+1} = [\mathfrak{g}, \mathfrak{g}_{k}]

            and recall that `\mathfrak{g}_{k} \supseteq \mathfrak{g}_{k+1}`.
            Alternatively we can express this as

            .. MATH::

                \mathfrak{g} \supseteq [\mathfrak{g}, \mathfrak{g}] \supseteq
                \bigl[ [\mathfrak{g}, \mathfrak{g}], \mathfrak{g} \bigr]
                \supseteq\biggl[\bigl[ [\mathfrak{g}, \mathfrak{g}],
                \mathfrak{g} \bigr], \mathfrak{g}\biggr] \supseteq \cdots.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.derived_series()
                (An example of a finite dimensional Lie algebra with basis:
                    the 3-dimensional abelian Lie algebra over Rational Field,
                 An example of a finite dimensional Lie algebra with basis:
                    the 0-dimensional abelian Lie algebra over Rational Field
                    with basis matrix:
                    [])

            ::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}})
                sage: L.lower_central_series() # todo: not implemented - #17416
                (Lie algebra on 2 generators (x, y) over Rational Field,
                 Subalgebra generated of Lie algebra on 2 generators (x, y) over Rational Field with basis:
                (x,))
            """
            L = [self]
            while L[-1].dimension() > 0:
                s = self.product_space(L[-1])
                if L[-1].dimension() == s.dimension():
                    break
                L.append(s)
            return tuple(L)

        def is_abelian(self):
            """
            Return if ``self`` is an abelian Lie algebra.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.is_abelian()
                True

            ::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'): {'x':1}})
                sage: L.is_abelian()
                False
            """
            return len(self.structure_coefficients()) == 0
            # TODO: boolean handling of empty family
            #return not self.structure_coefficients()

        def is_solvable(self):
            r"""
            Return if ``self`` is a solvable Lie algebra.

            A Lie algebra is solvable if the derived series eventually
            becomes `0`.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.is_solvable()
                True

            ::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}})
                sage: L.is_solvable() # todo: not implemented - #17416
                False
            """
            return not self.derived_series()[-1].dimension()

        def is_nilpotent(self):
            r"""
            Return if ``self`` is a nilpotent Lie algebra.

            A Lie algebra is nilpotent if the lower central series eventually
            becomes `0`.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.is_nilpotent()
                True
            """
            return not self.lower_central_series()[-1].dimension()

        def is_semisimple(self):
            """
            Return if ``self`` if a semisimple Lie algebra.

            A Lie algebra is semisimple if the solvable radical is zero. In
            characteristic 0, this is equivalent to saying the Killing form
            is non-degenerate.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.is_semisimple()
                False
            """
            return not self.killing_form_matrix().is_singular()

        @cached_method(key=lambda self,M,d,s,n: (M,d,s))
        def chevalley_eilenberg_complex(self, M=None, dual=False, sparse=True, ncpus=None):
            r"""
            Return the Chevalley-Eilenberg complex of ``self``.

            Let `\mathfrak{g}` be a Lie algebra and `M` be a right
            `\mathfrak{g}`-module. The *Chevalley-Eilenberg complex*
            is the chain complex on

            .. MATH::

                C_{\bullet}(\mathfrak{g}, M) =
                M \otimes \bigwedge\nolimits^{\bullet} \mathfrak{g},

            where the differential is given by

            .. MATH::

                d(m \otimes g_1 \wedge \cdots \wedge g_p) =
                \sum_{i=1}^p (-1)^{i+1}
                  (m g_i) \otimes g_1 \wedge \cdots \wedge
                  \hat{g}_i \wedge \cdots \wedge g_p +
                \sum_{1 \leq i < j \leq p} (-1)^{i+j}
                  m \otimes [g_i, g_j] \wedge
                  g_1 \wedge \cdots \wedge \hat{g}_i
                  \wedge \cdots \wedge \hat{g}_j
                  \wedge \cdots \wedge g_p.

            INPUT:

            - ``M`` -- (default: the trivial 1-dimensional module)
              the module `M`
            - ``dual`` -- (default: ``False``) if ``True``, causes
              the dual of the complex to be computed
            - ``sparse`` -- (default: ``True``) whether to use sparse
              or dense matrices
            - ``ncpus`` -- (optional) how many cpus to use

            EXAMPLES::

                sage: L = lie_algebras.sl(ZZ, 2)
                sage: C = L.chevalley_eilenberg_complex(); C
                Chain complex with at most 4 nonzero terms over Integer Ring
                sage: ascii_art(C)
                                          [ 2  0  0]       [0]
                                          [ 0 -1  0]       [0]
                            [0 0 0]       [ 0  0  2]       [0]
                 0 <-- C_0 <-------- C_1 <----------- C_2 <---- C_3 <-- 0

                sage: L = LieAlgebra(QQ, cartan_type=['C',2])
                sage: C = L.chevalley_eilenberg_complex()  # long time
                sage: [C.free_module_rank(i) for i in range(11)]  # long time
                [1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1]

            REFERENCES:

            - :wikipedia:`Lie_algebra_cohomology#Chevalley-Eilenberg_complex`
            - [Wei1994]_ Chapter 7

            .. TODO::

                Currently this is only implemented for coefficients
                given by the trivial module `R`, where `R` is the
                base ring and `g R = 0` for all `g \in \mathfrak{g}`.
                Allow generic coefficient modules `M`.
            """
            if dual:
                return self.chevalley_eilenberg_complex(M, dual=False,
                                                        sparse=sparse,
                                                        ncpus=ncpus).dual()

            if M is not None:
                raise NotImplementedError("only implemented for the default"
                                          " (the trivial module)")

            from itertools import combinations
            from sage.functions.other import binomial
            from sage.matrix.matrix_space import MatrixSpace
            R = self.base_ring()
            zero = R.zero()
            mone = -R.one()
            if M is not None:
                raise NotImplementedError("coefficient module M cannot be passed")

            # Make sure we specify the ordering of the basis
            B = self.basis()
            K = list(B.keys())
            B = [B[k] for k in K]
            Ind = list(range(len(K)))

            def sgn(k, X):
                """
                Insert a new entry ``k`` into a strictly increasing
                list ``X`` in such a way that the resulting list is
                still strictly increasing.
                The return value is the pair ``(s, Y)``, where ``Y``
                is the resulting list (as tuple) and ``s`` is the
                Koszul sign incurred by the insertion (with the
                understanding that ``k`` originally stood to the
                left of the list).
                If ``k`` is already in ``X``, then the return value
                is ``(zero, None)``.
                """
                Y = list(X)
                for i in range(len(X)-1, -1, -1):
                    val = X[i]
                    if val == k:
                        return zero, None
                    if k > val:
                        Y.insert(i+1, k)
                        return mone**(i+1), tuple(Y)
                Y.insert(0, k)
                return R.one(), tuple(Y)

            from sage.parallel.decorate import parallel
            @parallel(ncpus=ncpus)
            def compute_diff(k):
                """
                Build the ``k``-th differential (in parallel).
                """
                indices = {tuple(X): i for i,X in enumerate(combinations(Ind, k-1))}
                if sparse:
                    data = {}
                    row = 0
                else:
                    data = []
                for X in combinations(Ind, k):
                    if not sparse:
                        ret = [zero] * len(indices)
                    for i in range(k):
                        Y = list(X)
                        Y.pop(i)
                        # We do mone**i because we are 0-based
                        # This is where we would do the action on
                        #   the coefficients module
                        #ret[indices[tuple(Y)]] += mone**i * zero
                        for j in range(i+1,k):
                            # We shift j by 1 because we already removed
                            #   an earlier element from X.
                            Z = tuple(Y[:j-1] + Y[j:])
                            elt = mone**(i+j) * B[X[i]].bracket(B[X[j]])
                            for key, coeff in elt.to_vector().iteritems():
                                s, A = sgn(key, Z)
                                if A is None:
                                    continue
                                if sparse:
                                    coords = (row,indices[A])
                                    if coords in data:
                                        data[coords] += s * coeff
                                    else:
                                        data[coords] = s * coeff
                                else:
                                    ret[indices[A]] += s * coeff
                    if sparse:
                        row += 1
                    else:
                        data.append(ret)
                nrows = binomial(len(Ind), k)
                ncols = binomial(len(Ind), k-1)
                MS = MatrixSpace(R, nrows, ncols, sparse=sparse)
                ret = MS(data).transpose()
                ret.set_immutable()
                return ret

            chain_data = {X[0][0]: M for X, M in compute_diff(list( range(1,len(Ind)+1) ))}

            from sage.homology.chain_complex import ChainComplex
            try:
                return ChainComplex(chain_data, degree_of_differential=-1)
            except TypeError:
                return chain_data

        def homology(self, deg=None, M=None, sparse=True, ncpus=None):
            r"""
            Return the Lie algebra homology of ``self``.

            The Lie algebra homology is the homology of the
            Chevalley-Eilenberg chain complex.

            INPUT:

            - ``deg`` -- the degree of the homology (optional)
            - ``M`` -- (default: the trivial module) a right module
              of ``self``
            - ``sparse`` -- (default: ``True``) whether to use sparse
              matrices for the Chevalley-Eilenberg chain complex
            - ``ncpus`` -- (optional) how many cpus to use when
              computing the Chevalley-Eilenberg chain complex

            EXAMPLES::

                sage: L = lie_algebras.cross_product(QQ)
                sage: L.homology()
                {0: Vector space of dimension 1 over Rational Field,
                 1: Vector space of dimension 0 over Rational Field,
                 2: Vector space of dimension 0 over Rational Field,
                 3: Vector space of dimension 1 over Rational Field}

                sage: L = lie_algebras.pwitt(GF(5), 5)
                sage: L.homology()
                {0: Vector space of dimension 1 over Finite Field of size 5,
                 1: Vector space of dimension 0 over Finite Field of size 5,
                 2: Vector space of dimension 1 over Finite Field of size 5,
                 3: Vector space of dimension 1 over Finite Field of size 5,
                 4: Vector space of dimension 0 over Finite Field of size 5,
                 5: Vector space of dimension 1 over Finite Field of size 5}

                sage: d = {('x', 'y'): {'y': 2}}
                sage: L.<x,y> = LieAlgebra(ZZ, d)
                sage: L.homology()
                {0: Z, 1: Z x C2, 2: 0}

            .. SEEALSO::

                :meth:`chevalley_eilenberg_complex`
            """
            C = self.chevalley_eilenberg_complex(M=M, sparse=sparse,
                                                 ncpus=ncpus)
            return C.homology(deg=deg)

        def cohomology(self, deg=None, M=None, sparse=True, ncpus=None):
            r"""
            Return the Lie algebra cohomology of ``self``.

            The Lie algebra cohomology is the cohomology of the
            Chevalley-Eilenberg cochain complex (which is the dual
            of the Chevalley-Eilenberg chain complex).

            Let `\mathfrak{g}` be a Lie algebra and `M` a left
            `\mathfrak{g}`-module. It is known that `H^0(\mathfrak{g}; M)`
            is the subspace of `\mathfrak{g}`-invariants of `M`:

            .. MATH::

                H^0(\mathfrak{g}; M) = M^{\mathfrak{g}}
                = \{ m \in M \mid g m = 0
                    \text{ for all } g \in \mathfrak{g} \}.

            Additionally, `H^1(\mathfrak{g}; M)` is the space of
            derivations `\mathfrak{g} \to M`
            modulo the space of inner derivations, and
            `H^2(\mathfrak{g}; M)` is the space of equivalence classes
            of Lie algebra extensions of `\mathfrak{g}` by `M`.

            INPUT:

            - ``deg`` -- the degree of the homology (optional)
            - ``M`` -- (default: the trivial module) a right module
              of ``self``
            - ``sparse`` -- (default: ``True``) whether to use sparse
              matrices for the Chevalley-Eilenberg chain complex
            - ``ncpus`` -- (optional) how many cpus to use when
              computing the Chevalley-Eilenberg chain complex

            EXAMPLES::

                sage: L = lie_algebras.so(QQ, 4)
                sage: L.cohomology()
                {0: Vector space of dimension 1 over Rational Field,
                 1: Vector space of dimension 0 over Rational Field,
                 2: Vector space of dimension 0 over Rational Field,
                 3: Vector space of dimension 2 over Rational Field,
                 4: Vector space of dimension 0 over Rational Field,
                 5: Vector space of dimension 0 over Rational Field,
                 6: Vector space of dimension 1 over Rational Field}

                sage: L = lie_algebras.Heisenberg(QQ, 2)
                sage: L.cohomology()
                {0: Vector space of dimension 1 over Rational Field,
                 1: Vector space of dimension 4 over Rational Field,
                 2: Vector space of dimension 5 over Rational Field,
                 3: Vector space of dimension 5 over Rational Field,
                 4: Vector space of dimension 4 over Rational Field,
                 5: Vector space of dimension 1 over Rational Field}

                sage: d = {('x', 'y'): {'y': 2}}
                sage: L.<x,y> = LieAlgebra(ZZ, d)
                sage: L.cohomology()
                {0: Z, 1: Z, 2: C2}

            .. SEEALSO::

                :meth:`chevalley_eilenberg_complex`

            REFERENCES:

            - :wikipedia:`Lie_algebra_cohomology`
            """
            C = self.chevalley_eilenberg_complex(M=M, dual=True, sparse=sparse,
                                                 ncpus=ncpus)
            return C.homology(deg=deg)

        def as_finite_dimensional_algebra(self):
            """
            Return ``self`` as a :class:`FiniteDimensionalAlgebra`.

            EXAMPLES::

                sage: L = lie_algebras.cross_product(QQ)
                sage: x,y,z = L.basis()
                sage: F = L.as_finite_dimensional_algebra()
                sage: X,Y,Z = F.basis()
                sage: x.bracket(y)
                Z
                sage: X * Y
                Z
            """
            K = self._basis_ordering
            B = self.basis()
            mats = []
            R = self.base_ring()
            S = dict(self.structure_coefficients())
            V = self._dense_free_module()
            zero_vec = V.zero()
            for k in K:
                M = []
                for kp in K:
                    if (k, kp) in S:
                        M.append( -S[k,kp].to_vector() )
                    elif (kp, k) in S:
                        M.append( S[kp,k].to_vector() )
                    else:
                        M.append( zero_vec )
                mats.append(matrix(R, M))
            from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra import FiniteDimensionalAlgebra
            return FiniteDimensionalAlgebra(R, mats, names=self._names)

    class ElementMethods:
        def adjoint_matrix(self): # In #11111 (more or less) by using matrix of a morphism
            """
            Return the matrix of the adjoint action of ``self``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.an_element().adjoint_matrix()
                [0 0 0]
                [0 0 0]
                [0 0 0]

            ::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}})
                sage: x.adjoint_matrix()
                [0 0]
                [1 0]
                sage: y.adjoint_matrix()
                [-1  0]
                [ 0  0]
            """
            P = self.parent()
            basis = P.basis()
            return matrix(self.base_ring(),
                          [P.bracket(self, b).to_vector() for b in basis])

        def to_vector(self):
            """
            Return the vector in ``g.module()`` corresponding to the
            element ``self`` of ``g`` (where ``g`` is the parent of
            ``self``).

            Implement this if you implement ``g.module()``.
            See :meth:`sage.categories.lie_algebras.LieAlgebras.module`
            for how this is to be done.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.an_element().to_vector()
                (0, 0, 0)

                sage: D = DescentAlgebra(QQ, 4).D()
                sage: L = LieAlgebra(associative=D)
                sage: L.an_element().to_vector()
                (1, 1, 1, 1, 1, 1, 1, 1)

            TESTS:

            Check that the error raised agrees with the one
            from ``monomial_coefficients()`` (see :trac:`25007`)::

                sage: L = lie_algebras.sp(QQ, 4, representation='matrix')
                sage: x = L.an_element()
                sage: x.monomial_coefficients()
                Traceback (most recent call last):
                ...
                NotImplementedError: the basis is not defined
                sage: x.to_vector()
                Traceback (most recent call last):
                ...
                NotImplementedError: the basis is not defined
            """
            mc = self.monomial_coefficients(copy=False)
            M = self.parent().module()
            B = M.basis()
            return M.sum(mc[k] * B[i] for i,k in enumerate(self.parent()._basis_ordering)
                         if k in mc)

        _vector_ = to_vector

    class Subobjects(SubobjectsCategory):
        """
        A category for subalgebras of a finite dimensional Lie algebra
        with basis.
        """
        class ParentMethods:
            @abstract_method
            def ambient(self):
                """
                Return the ambient Lie algebra of ``self``.

                EXAMPLES::

                    sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                    sage: a, b, c = L.lie_algebra_generators()
                    sage: S = L.subalgebra([2*a+b, b + c])
                    sage: S.ambient() == L
                    True
                """

            @abstract_method
            def basis_matrix(self):
                """
                Return the basis matrix of ``self``.

                EXAMPLES::

                    sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                    sage: a, b, c = L.lie_algebra_generators()
                    sage: S = L.subalgebra([2*a+b, b + c])
                    sage: S.basis_matrix()
                    [   1    0 -1/2]
                    [   0    1    1]
                """

