r"""
Finite dimensional algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2011-2014 Nicolas M. Thiery <nthiery at users.sf.net>
#                2011-2014 Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

import operator
from sage.misc.cachefunc import cached_method, cached_function
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.algebras import Algebras
from sage.categories.associative_algebras import AssociativeAlgebras
from sage.categories.semisimple_algebras import SemisimpleAlgebras
from sage.matrix.constructor import Matrix
from sage.functions.other import sqrt


class FiniteDimensionalAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    r"""
    The category of finite dimensional algebras with a distinguished basis.

    EXAMPLES::

        sage: C = FiniteDimensionalAlgebrasWithBasis(QQ); C
        Category of finite dimensional algebras with basis over Rational Field
        sage: C.super_categories()
        [Category of algebras with basis over Rational Field,
         Category of finite dimensional modules with basis over Rational Field]
        sage: C.example()
        An example of a finite dimensional algebra with basis:
        the path algebra of the Kronecker quiver
        (containing the arrows a:x->y and b:x->y) over Rational Field

    TESTS::

        sage: TestSuite(C).run()
        sage: C is Algebras(QQ).FiniteDimensional().WithBasis()
        True
        sage: C is Algebras(QQ).WithBasis().FiniteDimensional()
        True
    """

    class ParentMethods:

        @cached_method
        def radical_basis(self, cache_products=True):
            r"""
            Return a basis of the Jacobson radical of this algebra.

            .. NOTE::

               This implementation also works for algebras over fields of
               finite characteristic `p` in which we can compute `x^{1/p}`.

            INPUT:

            - ``cache_products`` -- boolean (default: ``True``); if ``True``
              then all products computed in this method are cached.

            OUTPUT:

            - ``list`` of elements of ``self``

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: A.radical_basis()
                [a, b]

            We construct the group algebra of the Klein Four-Group over the rationals::

                sage: A = KleinFourGroup().algebra(QQ)

            This algebra belongs to the category of finite dimensional
            algebras over the rationals::

                sage: A in FiniteDimensionalAlgebrasWithBasis(QQ)
                True

            Since the field has characteristic `0`, Maschke's Theorem
            tells us that the group algebra is semisimple. So its
            radical is the zero ideal::

                sage: A.radical_basis()
                []

            Let's work instead over a field of characteristic `2`::

                sage: A = KleinFourGroup().algebra(GF(2))
                sage: A.radical_basis()
                [B[()] + B[(1,2)(3,4)], B[(3,4)] + B[(1,2)(3,4)], B[(1,2)] + B[(1,2)(3,4)]]

            TESTS::

                sage: A = KleinFourGroup().algebra(GF(2))
                sage: A.radical_basis(cache_products=True)
                [B[()] + B[(1,2)(3,4)], B[(3,4)] + B[(1,2)(3,4)], B[(1,2)] + B[(1,2)(3,4)]]
                sage: A.radical_basis(cache_products=False)
                [B[()] + B[(1,2)(3,4)], B[(3,4)] + B[(1,2)(3,4)], B[(1,2)] + B[(1,2)(3,4)]]

            ::

                sage: A = KleinFourGroup().algebra(QQ)
                sage: A.radical_basis(cache_products=True)
                []
                sage: A.radical_basis(cache_products=False)
                []

            .. TODO:: explain and check this example

            ::

                sage: class AnAlgebra(CombinatorialFreeModule):
                ...       def __init__(self, F):
                ...           R.<x> = PolynomialRing(F)
                ...           I = R.ideal(x**F.characteristic()-F.one())
                ...           self._xbar = R.quotient(I).gen()
                ...           basis_keys = [self._xbar**i for i in range(F.characteristic())]
                ...           CombinatorialFreeModule.__init__(self, F, basis_keys,
                ...                   category = FiniteDimensionalAlgebrasWithBasis(F))
                ...       def one(self):
                ...           return self.basis()[self.base_ring().one()]
                ...       def product_on_basis(self, w1, w2):
                ...           return self.from_vector(vector(w1*w2))
                sage: AnAlgebra(GF(3)).radical_basis()
                [B[1] + 2*B[xbar^2], B[xbar] + 2*B[xbar^2]]
                sage: AnAlgebra(GF(16,'a')).radical_basis()
                [B[1] + B[xbar]]
                sage: AnAlgebra(GF(49,'a')).radical_basis()
                [B[1] + 6*B[xbar^6], B[xbar] + 6*B[xbar^6], B[xbar^2] + 6*B[xbar^6], B[xbar^3] + 6*B[xbar^6], B[xbar^4] + 6*B[xbar^6], B[xbar^5] + 6*B[xbar^6]]

            .. SEEALSO:: :meth:`radical`, :class:`SemisimpleAlgebras`

            AUTHORS: TODO: polish this!

            - Franco Saliola
            """
            F = self.base_ring()
            if not F.is_field():
                raise NotImplementedError, "the base ring must be a field"
            p = F.characteristic()
            from sage.matrix.constructor import matrix
            from sage.modules.free_module_element import vector

            if self in SemisimpleAlgebras(self.base_ring()):
                return []

            if cache_products is True:
                product_on_basis = cached_function(self.product_on_basis)
            else:
                product_on_basis = self.product_on_basis

            if p == 0:
                keys = self.basis().keys()
                mat = matrix(self.base_ring(), [
                        [sum(product_on_basis(x,j).coefficient(i) * product_on_basis(y,i).coefficient(j)
                            for i in keys for j in keys) for x in keys] for y in keys
                        ])
                rad_basis = mat.kernel().basis()
            else:
                # TODO: some finite field elements in Sage have both an
                # ``nth_root`` method and a ``pth_root`` method (such as ``GF(9,'a')``),
                # some only have a ``nth_root`` element such as ``GF(2)``
                # I imagine that ``pth_root`` would be fastest, but it is not
                # always available....
                if hasattr(self.base_ring().one(), 'nth_root'):
                    root_fcn = lambda s, x : x.nth_root(s)
                else:
                    root_fcn = lambda s, x : x**(1/s)

                s, n = 1, self.dimension()
                B = [b.on_left_matrix() for b in self.basis()]
                I = B[0].parent().one()
                while s <= n:
                    BB = B + [I]
                    G = matrix([ [(-1)**s * (b*bb).characteristic_polynomial()[n-s]
                                    for bb in BB] for b in B])
                    C = G.left_kernel().basis()
                    if 1 < s < F.order():
                        C = [vector(F, [root_fcn(s, ci) for ci in c]) for c in C]
                    B = [ sum(ci*b for (ci,b) in zip(c,B)) for c in C ]
                    s = p * s
                e = vector(self.one())
                rad_basis = [b*e for b in B]

            return [self.from_vector(vec) for vec in rad_basis]

        @cached_method
        def radical(self):
            r"""
            Return the Jacobson radical of ``self``.

            .. SEEALSO:: :meth:`radical_basis`, :meth:`semisimple_quotient`

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: radical = A.radical(); radical
                Radical of An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: from sage.categories.associative_algebras import AssociativeAlgebras
                sage: radical in AssociativeAlgebras(QQ).WithBasis().FiniteDimensional()
                True
                sage: radical.dimension()
                2
                sage: radical.basis()
                Finite family {0: B[0], 1: B[1]}
                sage: radical.ambient() is A
                True
                sage: [c.lift() for c in radical.basis()]
                [a, b]

            .. TODO::

                - This is in fact an ideal.
                - Add references
                - Pickling by construction, as ``A.center()``
                - Lazy evaluation of ``_repr_``

            TESTS::

                sage: TestSuite(radical).run()
            """
            category = AssociativeAlgebras(self.base_ring()).WithBasis().FiniteDimensional().Subobjects()
            radical = self.submodule(self.radical_basis(),
                                     category = category,
                                     already_echelonized = True)
            radical.rename("Radical of {}".format(self))
            return radical

        @cached_method
        def semisimple_quotient(self):
            """
            Return the semisimple quotient of ``self``.

            This is the quotient of ``self`` by its radical.

            .. SEEALSO:: :meth:`radical`

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example(); A
                An example of a finite dimensional algebra with basis: the path algebra of the Kronecker quiver (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: a,b,x,y = sorted(A.basis())
                sage: S = A.semisimple_quotient(); S
                Semisimple quotient of An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: S in SemisimpleAlgebras
                True
                sage: S.basis()
                Finite family {'y': B['y'], 'x': B['x']}
                sage: xs,ys = sorted(S.basis())
                sage: (xs + ys) * xs
                B['x']

            .. TODO::

               - This example is not very interesting because the
                 semisimple quotient is actually a subalgebra.
               - Pickling by construction, as ``A.center()``
               - Lazy evaluation of ``_repr_``

            TESTS::

                sage: TestSuite(S).run()
            """
            ring = self.base_ring()
            category=Algebras(ring).WithBasis().FiniteDimensional().Quotients() \
                      & SemisimpleAlgebras(ring)
            result = self.quotient(self.radical(), category=category)
            result.rename("Semisimple quotient of {}".format(self))
            return result


        @cached_method
        def center_basis(self):
            r"""
            Return a basis of the center of ``self``.

            OUTPUT:

            - ``list`` of elements of ``self``

            .. SEEALSO:: :meth:`center`

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example(); A
                An example of a finite dimensional algebra with basis: the path algebra of the Kronecker quiver (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: A.center_basis()
                [x + y]
            """
            return self.annihilator_basis(self.algebra_generators(), self.bracket)

        @cached_method
        def center(self):
            r"""
            Return the center of ``self``.

            .. SEEALSO:: :meth:`center_basis`

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: center = A.center(); center
                Center of An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: center in Algebras(QQ).WithBasis().FiniteDimensional().Commutative()
                True
                sage: center.dimension()
                1
                sage: center.basis()
                Finite family {0: B[0]}
                sage: center.ambient() is A
                True
                sage: [c.lift() for c in center.basis()]
                [x + y]
                sage: DihedralGroup(6).algebra(QQ).center() in SemisimpleAlgebras
                True

            .. TODO::

                - Pickling by construction, as ``A.center()``
                - Lazy evaluation of ``_repr_``

            TESTS::

                sage: TestSuite(center).run()
            """
            if self in SemisimpleAlgebras(self.base_ring()):
                category = SemisimpleAlgebras(self.base_ring()).FiniteDimensional().WithBasis().Subobjects().Commutative()
            else:
                category = Algebras(self.base_ring()).FiniteDimensional().WithBasis().Subobjects().Commutative()
            center = self.submodule(self.center_basis(),
                                    category = category,
                                    already_echelonized = True)
            center.rename("Center of {}".format(self))
            return center

        def principal_ideal(self, a, side='left'):
            r"""
            Construct the ``side`` A-module generated by ``a``.

            EXAMPLE::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: A.principal_ideal(A.an_element())
                Free module generated by {0, 1, 2, 3} over Rational Field

            """
            B = self.basis()
            if side == 'left':
                phi = self.module_morphism(on_basis=lambda i: a*B[i],
                        codomain=self,
                        triangular=True)
            elif side == 'right':
                phi = self.module_morphism(on_basis=lambda i: B[i]*a,
                        codomain=self,
                        triangluar=True)
            else:
                raise Exception("Side must be ``left`` or ``right``")
            ideal = phi.matrix().image()
            #now let's make it a submodule of A
            return self.submodule([self.from_vector(v) for v in ideal.basis()],
                    already_echelonized=True)

        def orthogonal_idempotents(self):
            r"""
            Return a maximal family of orthogonal idempotents of ``self``.

            INPUT:

            - ``self`` -- a finite dimensional algebra

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: A
                An example of a finite dimensional algebra with basis: the path
                algebra of the Kronecker quiver (containing the arrows a:x->y
                and b:x->y) over Rational Field
                sage: sorted(A.orthogonal_idempotents(), key=str)
                [x, y]
                sage: Monoids().Finite().example()
                An example of a finite multiplicative monoid: the integers
                modulo 12
                sage: Z12 = Monoids().Finite().example()
                sage: A = Z12.algebra(QQ)
                sage: sorted(A.orthogonal_idempotents(), key=str)
                [-1/2*B[3] + 1/2*B[9], -1/2*B[8] + 1/2*B[4], -B[0] + 1/2*B[3] +
                1/2*B[9], 1/2*B[8] + 1/2*B[4] - B[0], 1/4*B[1] + 1/2*B[3] +
                1/4*B[5] - 1/4*B[7] - 1/2*B[9] - 1/4*B[11], 1/4*B[1] + 1/4*B[11]
                - 1/4*B[5] - 1/4*B[7], 1/4*B[1] - 1/2*B[4] - 1/4*B[5] + 1/4*B[7]
                + 1/2*B[8] - 1/4*B[11], B[0], B[0] + 1/4*B[1] - 1/2*B[3] -
                1/2*B[4] + 1/4*B[5] + 1/4*B[7] - 1/2*B[8] - 1/2*B[9] +
                1/4*B[11]]


            """
            Aquo = self.semisimple_quotient()
            orth_quo = Aquo.orthogonal_idempotents()
            return [self._lift_idempotent(x) for x in orth_quo]

        def _lift_idempotent(self, x):
            r"""
            Lift an idempotent of the semisimple quotient of ``self`` into an
            idempotent of ``self``.

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: Aquo = A.semisimple_quotient()
                sage: orth = Aquo.orthogonal_idempotents()
                sage: A._lift_idempotent(orth[1])
                y
            """
            idempOld = None
            assert x in self.semisimple_quotient()
            idemp = x.lift()
            p = idemp.parent()
            while idemp <> idempOld:
                tmp = idemp
                idemp = (p.one() - (p.one() - idemp**2)**2)
                idempOld = tmp
            return idemp

        @cached_method
        def cartan_invariant_matrix(self, side='left'):
            r"""
            Return the Cartan invariant matrix of the algebra ``self``.

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: A.cartan_invariant_matrix()
                [1 0]
                [2 1]
                sage: A3 = SymmetricGroup(3).algebra(QQ)
                sage: A3.cartan_invariant_matrix()
                [1 0 0]
                [0 1 0]
                [0 0 1]
                sage: Z12 = Monoids().Finite().example()
                sage: A = Z12.algebra(QQ)
                sage: A.cartan_invariant_matrix()
                [1 0 0 0 0 0 0 0 0]
                [0 1 0 0 0 0 0 0 0]
                [0 0 2 0 0 0 0 0 0]
                [0 0 0 1 0 0 0 0 0]
                [0 0 0 0 1 0 0 0 0]
                [0 0 0 0 0 2 0 0 0]
                [0 0 0 0 0 0 1 0 0]
                [0 0 0 0 0 0 0 1 0]
                [0 0 0 0 0 0 0 0 2]
            """
            Aquo = self.semisimple_quotient()
            orth_quo = Aquo.orthogonal_idempotents()
            # Dimension of simple modules
            dimSimples = [sqrt(Aquo.principal_ideal(e, side).dimension()) for e in
                    orth_quo]
            orth = [x._lift_idempotent() for x in orth_quo]
            return Matrix(self.base_ring(),
                    len(orth),
                    lambda i,j: self._cartan_matrix_coef(orth[i], orth[j])/(dimSimples[i]*dimSimples[j]))

        def projective_decomposition(self, side='left'):
            r"""
            Return the list of indecomposable projective ``side``-sided
            ``self``-modules.

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: projs = A.projective_decomposition()
                sage: sorted(projs, key=str)
                [Free module generated by {0, 1, 2} over Rational Field,
                 Free module generated by {0} over Rational Field]

            We check that the sum of the dimensions of the indecomposable
            projective module is the dimension of ``self``:

                sage: sum([P.dimension() for P in projs]) == A.dimension()
                True

            """
            return [self.principal_ideal(e, side) for e in self.orthogonal_idempotents()]


        def _cartan_matrix_coef(self, ei, ej):
            B = self.basis()
            phi = self.module_morphism(on_basis=lambda k: ei * B[k] * ej,
                    codomain=self,
                    triangular=True)
            return phi.matrix().rank()

    class ElementMethods:

        def to_matrix(self, base_ring=None, action=operator.mul, side='left'):
            """
            Return the matrix of the action of ``self`` on the algebra.

            INPUT::

            - ``base_ring`` -- the base ring for the matrix to be constructed
            - ``action`` -- a function (default: :func:`operator.mul`)
            - ``side`` -- 'left' or 'right' (default: 'right')

            EXAMPLES::

                sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
                sage: a = QS3([2,1,3])
                sage: a.to_matrix(side='left')
                [0 0 1 0 0 0]
                [0 0 0 0 1 0]
                [1 0 0 0 0 0]
                [0 0 0 0 0 1]
                [0 1 0 0 0 0]
                [0 0 0 1 0 0]
                sage: a.to_matrix(base_ring=RDF, side="left")
                [0.0 0.0 1.0 0.0 0.0 0.0]
                [0.0 0.0 0.0 0.0 1.0 0.0]
                [1.0 0.0 0.0 0.0 0.0 0.0]
                [0.0 0.0 0.0 0.0 0.0 1.0]
                [0.0 1.0 0.0 0.0 0.0 0.0]
                [0.0 0.0 0.0 1.0 0.0 0.0]

            AUTHORS: Mike Hansen # TODO: polish this!
            """
            basis = self.parent().basis()
            action_left = action
            if side == 'right':
                action = lambda x: action_left(basis[x], self)
            else:
                action = lambda x: action_left(self, basis[x])
            endo = self.parent().module_morphism(on_basis=action, codomain=self.parent())
            return endo.matrix(base_ring=base_ring)

        _matrix_ = to_matrix  # For temporary backward compatibility
        on_left_matrix = to_matrix
