r"""
Finite Dimensional Lie Algebras With Basis

AUTHORS:

- Travis Scrimshaw (07-15-2013): Initial implementation
"""
#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.lie_algebras import LieAlgebras
from sage.categories.subobjects import SubobjectsCategory
from sage.rings.all import ZZ
from sage.algebras.free_algebra import FreeAlgebra
from sage.sets.family import Family
from sage.matrix.constructor import matrix

class FiniteDimensionalLieAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    Category of finite dimensional Lie algebras with a basis.

    .. TODO::

        Many of these tests should use non-abelian Lie algebras and need to
        be added after :trac:`16820`.
    """
    _base_category_class_and_axiom = [LieAlgebras.FiniteDimensional, "WithBasis"]

    def example(self, names=('a', 'b', 'c')):
        """
        Return an example of a finite dimensional Lie algebra with basis as per
        :meth:`Category.example <sage.categories.category.Category.example>`.

        EXAMPLES::

            sage: C = LieAlgebras(QQ).FiniteDimensional().WithBasis()
            sage: C.example()
            An example of a finite dimensional Lie algebra with basis:
             the abelian Lie algebra with generators ('a', 'b', 'c')
             over Rational Field

        Other names of generators can be specified as an optional argument::

            sage: C.example(('x','y','z'))
            An example of a finite dimensional Lie algebra with basis:
             the abelian Lie algebra with generators ('x', 'y', 'z')
             over Rational Field
        """
        from sage.categories.examples.finite_dimensional_lie_algebras_with_basis import Example
        return Example(self.base_ring(), names)

    class ParentMethods:
        @cached_method
        def _construct_UEA(self):
            """
            Construct the universal enveloping algebra of ``self``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: UEA = L._construct_UEA(); UEA
                Noncommutative Multivariate Polynomial Ring in a, b, c
                 over Rational Field, nc-relations: {}
                sage: UEA.relations(add_commutative=True)
                {b*a: a*b, c*a: a*c, c*b: b*c}

            ::

                sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}, ('y','z'):{'x':1}, ('z','x'):{'y':1}}) # todo: not implemented - #16820
                sage: L._construct_UEA() # todo: not implemented - #16820
            """
            # Create the UEA relations
            # We need to get names for the basis elements, not just the generators
            try:
                names = self.variable_names()
            except ValueError:
                names = tuple('b{}'.format(i) for i in range(self.dimension()))
            I = self._ordered_indices
            F = FreeAlgebra(self.base_ring(), names)
            #gens = F.gens()
            d = F.gens_dict()
            rels = {}
            S = self.structure_coefficients(True)
            for k in S.keys():
                g0 = d[names[I.index(k[0])]]
                g1 = d[names[I.index(k[1])]]
                if g0 < g1:
                    rels[g1*g0] = g0*g1 - sum(val*d[names[I.index(g)]] for g, val in S[k])
                else:
                    rels[g0*g1] = g1*g0 + sum(val*d[names[I.index(g)]] for g, val in S[k])
            return F.g_algebra(rels)

        def killing_matrix(self, x, y):
            r"""
            Return the Killing matrix of ``x`` and ``y``.

            The Killing matrix is defined as the matrix corresponding
            to the action of `\mathrm{ad}_x \circ \mathrm{ad}_y` in
            the basis of ``self``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: a,b,c = L.lie_algebra_generators()
                sage: L.killing_matrix(a, b)
                [0 0 0]
                [0 0 0]
                [0 0 0]

            ::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}}) # todo: not implemented - #16820
                sage: L.killing_matrix(x, y) # todo: not implemented - #16820
                [ 0  0]
                [-1  0]
            """
            return x.adjoint_matrix() * y.adjoint_matrix()

        def killing_form(self, x, y):
            r"""
            Return the Killing form on ``x`` and ``y``.

            The Killing form is defined as

            .. MATH::

                \langle x \mid y \rangle = \mathrm{tr}\left( \mathrm{ad}_x
                \circ \mathrm{ad}_y \right).

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

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.killing_form_matrix()
                [0 0 0]
                [0 0 0]
                [0 0 0]
            """
            B = self.basis()
            m = matrix([[self.killing_form(x, y) for x in B] for y in B])
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

            A dictionary whose keys are pairs of basis indices `(i, j)` and
            whose values are the corresponding *element* of `[b_i, b_j]`
            in the Lie algebra.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.structure_coefficients()
                Finite family {}
                sage: L.structure_coefficients(True)
                Finite family {('b', 'c'): (0, 0, 0), ('a', 'b'): (0, 0, 0),
                               ('a', 'c'): (0, 0, 0)}
            """
            d = {}
            B = self.basis()
            K = self.basis().keys()
            zero = self.zero()
            one = self.base_ring().one()
            for i,x in enumerate(K):
                for y in K[i+1:]:
                    bx = B[x]
                    by = B[x]
                    val = self.bracket(bx, by)
                    if not include_zeros and val == zero:
                        continue
                    if self._basis_cmp(x, y) > 0:
                        d[(y, x)] = -val
                    else:
                        d[(x, y)] = val
            return Family(d)

        def centralizer(self, S):
            """
            Return the centralizer of ``S`` in ``self``.

            .. TODO::

                What is S?

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: a,b,c = L.lie_algebra_generators()
                sage: S = L.centralizer([a + b, 2*a + c]); S
                An example of a finite dimensional Lie algebra with basis:
                 the abelian Lie algebra with generators ('x0', 'x1', 'x2')
                 over Rational Field
                sage: S.basis_matrix()
                [1 0 0]
                [0 1 0]
                [0 0 1]
            """
            #from sage.algebras.lie_algebras.subalgebra import LieSubalgebra
            #if isinstance(S, LieSubalgebra) or S is self:
            if S is self:
                K = S
            else:
                K = self.subalgebra(S)

            m = K.basis_matrix()
            S = self.structure_coefficients()
            sc = {k: S[k].to_vector() for k in S.keys()}
            X = self.basis().keys()
            d = self.dimension()
            K = sc.keys()
            c_mat = matrix([[sum(r[j]*sc[x,X[j]][k] for j in range(d) if (x, X[j]) in K)
                             for x in X]
                            for r in m for k in range(d)])
            C = c_mat.right_kernel().basis_matrix()
            return self.subalgebra(map(self, C))

        def center(self):
            """
            Return the center of ``self``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: Z = L.center(); Z
                An example of a finite dimensional Lie algebra with basis:
                 the abelian Lie algebra with generators ('x0', 'x1', 'x2')
                 over Rational Field
                sage: Z.basis_matrix()
                [1 0 0]
                [0 1 0]
                [0 0 1]
            """
            return self.centralizer(self)

        def product_space(self, L):
            r"""
            Return the product space ``[self, L]``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: a,b,c = L.lie_algebra_generators()
                sage: X = L.subalgebra([a, b+c])
                sage: L.product_space(X)
                An example of a finite dimensional Lie algebra with basis:
                 the abelian Lie algebra with generators ()
                 over Rational Field with basis matrix:
                 []
                sage: Y = L.subalgebra([a, 2*b-c])
                sage: X.product_space(Y)
                An example of a finite dimensional Lie algebra with basis:
                 the abelian Lie algebra with generators ()
                 over Rational Field with basis matrix:
                 []

            ::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}}) # todo: not implemented - #16820
                sage: Lp = L.product_space(L) # todo: not implemented - #16820
                sage: Lp # todo: not implemented - #16820
                Subalgebra generated of Lie algebra on 2 generators (x, y) over Rational Field with basis:
                (x,)
                sage: Lp.product_space(L) # todo: not implemented - #16820
                Subalgebra generated of Lie algebra on 2 generators (x, y) over Rational Field with basis:
                (x,)
                sage: L.product_space(Lp) # todo: not implemented - #16820
                Subalgebra generated of Lie algebra on 2 generators (x, y) over Rational Field with basis:
                (x,)
                sage: Lp.product_space(Lp) # todo: not implemented - #16820
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
            K = B.keys()
            LK = LB.keys()
            # We echelonize the matrix here
            b_mat = matrix(A.base_ring(), [A.bracket(B[a], LB[b]).to_vector()
                                           for a in K for b in LK])
            b_mat.echelonize()
            r = b_mat.rank()
            I = A._ordered_indices
            gens = [A.element_class(A, {I[i]: v for i,v in row.iteritems()})
                    for row in b_mat.rows()[:r]]
            return A.subalgebra(gens)

        @cached_method
        def derived_subalgebra(self):
            """
            Return the derived subalgebra of ``self``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.derived_subalgebra()
                An example of a finite dimensional Lie algebra with basis:
                 the abelian Lie algebra with generators () over
                 Rational Field with basis matrix:
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
                    the abelian Lie algebra with generators ('a', 'b', 'c')
                    over Rational Field,
                 An example of a finite dimensional Lie algebra with basis:
                    the abelian Lie algebra with generators () over
                    Rational Field with basis matrix:
                    [])

            ::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}}) # todo: not implemented - #16820
                sage: L.derived_series() # todo: not implemented - #16820
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

            and recall that `\mathfrak{g}_k} \supseteq \mathfrak{g}_{k+1}`.
            Alternatively we can express this as

            .. MATH::

                \mathfrak{g} \supseteq [\mathfrak{g}, \mathfrak{g}] \supseteq
                \bigl[ [\mathfrak{g}, \mathfrak{g}], \mathfrak{g} \bigr]
                \supseteq\biggl[\bigl[ [\mathfrak{g}, \mathfrak{g}],
                \mathfrak{g} \bigr], \mathfrak{g}\biggr] \supseteq \cdots.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.lower_central_series()
                (An example of a finite dimensional Lie algebra with basis:
                    the abelian Lie algebra with generators ('a', 'b', 'c')
                    over Rational Field,
                 An example of a finite dimensional Lie algebra with basis:
                    the abelian Lie algebra with generators () over
                    Rational Field with basis matrix:
                    [])

            ::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}}) # todo: not implemented - #16820
                sage: L.lower_central_series() # todo: not implemented - #16820
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

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}}) # todo: not implemented - #16820
                sage: L.is_abelian() # todo: not implemented - #16820
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

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}}) # todo: not implemented - #16820
                sage: L.is_solvable() # todo: not implemented - #16820
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

        def dimension(self):
            """
            Return the dimension of ``self``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.dimension()
                3

            ::

                sage: L = LieAlgebra(QQ, 'x,y', {('x','y'): {'x':1}}) # todo: not implemented - #16820
                sage: L.dimension() # todo: not implemented - #16820
                2
            """
            return self.basis().cardinality()

    class ElementMethods:
        def to_vector(self):
            """
            Return ``self`` as a vector.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.an_element().to_vector()
                (0, 0, 0)
            """
            M = self.parent().free_module()
            if not self:
                return M.zero()
            B = M.basis()
            return M.sum(B[k]*self[k] for k in self.parent()._ordered_indices)

        def adjoint_matrix(self): # In #11111 (more or less) by using matrix of a mophism
            """
            Return the matrix of the adjoint action of ``self``.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                sage: L.an_element().adjoint_matrix()
                [0 0 0]
                [0 0 0]
                [0 0 0]

            ::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}}) # todo: not implemented - #16820
                sage: x.adjoint_matrix() # todo: not implemented - #16820
                [1 0]
                [0 0]
                sage: y.adjoint_matrix() # todo: not implemented - #16820
                [ 0  0]
                [-1  0]
            """
            P = self.parent()
            basis = P.basis()
            return matrix([P.bracket(self, b).to_vector() for b in basis])

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
                    sage: L.inject_variables()
                    Defining a, b, c
                    sage: S = L.subalgebra([2*a+b, b + c], 'x,y')
                    sage: S.ambient() == L
                    True
                """

            @abstract_method
            def basis_matrix(self):
                """
                Return the basis matrix of ``self``.

                EXAMPLES::

                    sage: L = LieAlgebras(QQ).FiniteDimensional().WithBasis().example()
                    sage: L.inject_variables()
                    Defining a, b, c
                    sage: S = L.subalgebra([2*a+b, b + c], 'x,y')
                    sage: S.basis_matrix()
                    [   1    0 -1/2]
                    [   0    1    1]
                """

