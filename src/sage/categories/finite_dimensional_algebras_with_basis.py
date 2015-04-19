# -*- coding: utf-8 -*-
r"""
Finite dimensional algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2011-2015 Nicolas M. Thiéry <nthiery at users.sf.net>
#                2011-2015 Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

import operator
from sage.misc.cachefunc import cached_method, cached_function
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.algebras import Algebras
from sage.categories.associative_algebras import AssociativeAlgebras
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

               This implementation handles algebras over fields of
               characteristic zero (using Dixon's lemma) or fields of
               characteristic `p` in which we can compute `x^{1/p}`
               [FR85], [Eb89].

            REFERENCES:

            [Eb89] Eberly, Wayne. "Computations for algebras and group
            representations." Ph.D. Thesis, University of Toronto, 1989.

            [FR85] Friedl, Katalin, and Lajos Rónyai. "Polynomial time
            solutions of some problems of computational algebra." Proceedings
            of the seventeenth annual ACM symposium on Theory of computing.
            ACM, 1985.

            INPUT:

            - ``cache_products`` -- a boolean (default: ``True``); if
              ``True`` then all products computed in this method are cached.

            OUTPUT:

            - a list of elements of ``self``.

            .. SEEALSO:: :meth:`radical`, :class:`Algebras.Semisimple`

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: A.radical_basis()
                [a, b]

            We construct the group algebra of the Klein Four-Group
            over the rationals::

                sage: A = KleinFourGroup().algebra(QQ)

            This algebra belongs to the category of finite dimensional
            algebras over the rationals::

                sage: A in FiniteDimensionalAlgebrasWithBasis(QQ)
                True

            Since the field has characteristic `0`, Maschke's Theorem
            tells us that the group algebra is semisimple. So its
            radical is the zero ideal::

                sage: A in Algebras(QQ).Semisimple()
                True
                sage: A.radical_basis()
                []

            Let's work instead over a field of characteristic `2`::

                sage: A = KleinFourGroup().algebra(GF(2))
                sage: A in Algebras(GF(2)).Semisimple()
                False
                sage: A.radical_basis()
                [B[()] + B[(1,2)(3,4)], B[(3,4)] + B[(1,2)(3,4)], B[(1,2)] + B[(1,2)(3,4)]]

            We now implement the algebra `A = K[x] / x^p-1`, where `K`
            is a finite field of characteristic `p`, and check its
            radical; alas, we currently need to wrap `A` to make it a
            proper :class:`ModulesWithBasis`::

                sage: class AnAlgebra(CombinatorialFreeModule):
                ....:     def __init__(self, F):
                ....:         R.<x> = PolynomialRing(F)
                ....:         I = R.ideal(x**F.characteristic()-F.one())
                ....:         self._xbar = R.quotient(I).gen()
                ....:         basis_keys = [self._xbar**i for i in range(F.characteristic())]
                ....:         CombinatorialFreeModule.__init__(self, F, basis_keys,
                ....:                 category=FiniteDimensionalAlgebrasWithBasis(F))
                ....:     def one(self):
                ....:         return self.basis()[self.base_ring().one()]
                ....:     def product_on_basis(self, w1, w2):
                ....:         return self.from_vector(vector(w1*w2))
                sage: AnAlgebra(GF(3)).radical_basis()
                [B[1] + 2*B[xbar^2], B[xbar] + 2*B[xbar^2]]
                sage: AnAlgebra(GF(16,'a')).radical_basis()
                [B[1] + B[xbar]]
                sage: AnAlgebra(GF(49,'a')).radical_basis()
                [B[1] + 6*B[xbar^6], B[xbar] + 6*B[xbar^6], B[xbar^2] + 6*B[xbar^6],
                 B[xbar^3] + 6*B[xbar^6], B[xbar^4] + 6*B[xbar^6], B[xbar^5] + 6*B[xbar^6]]

            TESTS::

                sage: A = KleinFourGroup().algebra(GF(2))
                sage: A.radical_basis(cache_products=True)
                [B[()] + B[(1,2)(3,4)], B[(3,4)] + B[(1,2)(3,4)], B[(1,2)] + B[(1,2)(3,4)]]
                sage: A.radical_basis(cache_products=False)
                [B[()] + B[(1,2)(3,4)], B[(3,4)] + B[(1,2)(3,4)], B[(1,2)] + B[(1,2)(3,4)]]

                sage: A = KleinFourGroup().algebra(QQ, category=Monoids())
                sage: A.radical_basis.__module__
                'sage.categories.finite_dimensional_algebras_with_basis'
                sage: A.radical_basis(cache_products=True)
                []
                sage: A.radical_basis(cache_products=False)
                []
            """
            F = self.base_ring()
            if not F.is_field():
                raise NotImplementedError("the base ring must be a field")
            p = F.characteristic()
            from sage.matrix.constructor import matrix
            from sage.modules.free_module_element import vector

            product_on_basis = self.product_on_basis
            if cache_products is True:
                product_on_basis = cached_function(product_on_basis)

            if p == 0:
                keys = list(self.basis().keys())
                cache = [{(i,j): c
                            for i in keys
                            for j,c in product_on_basis(y,i)}
                         for y in keys]
                mat = [ [ sum(x.get((j, i), 0) * c for (i,j),c in y.items())
                            for x in cache]
                        for y in cache]

                mat = matrix(self.base_ring(), mat)
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

            This uses :meth:`radical_basis`, whose default
            implementation handles algebras over fields of
            characteristic zero or fields of characteristic `p` in
            which we can compute `x^{1/p}`.

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

            The radical is an ideal of `A`, and thus a finite
            dimensional non unital associative algebra::

                sage: from sage.categories.associative_algebras import AssociativeAlgebras
                sage: radical in AssociativeAlgebras(QQ).WithBasis().FiniteDimensional()
                True
                sage: radical in Algebras(QQ)
                False

                sage: radical.dimension()
                2
                sage: radical.basis()
                Finite family {0: B[0], 1: B[1]}
                sage: radical.ambient() is A
                True
                sage: [c.lift() for c in radical.basis()]
                [a, b]

            .. TODO::

                - Tell Sage that the radical is in fact an ideal;
                - Pickling by construction, as ``A.center()``;
                - Lazy evaluation of ``_repr_``.

            TESTS::

                sage: TestSuite(radical).run()
            """
            category = AssociativeAlgebras(self.base_ring()).WithBasis().FiniteDimensional().Subobjects()
            radical = self.submodule(self.radical_basis(),
                                     category=category,
                                     already_echelonized=True)
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
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: a,b,x,y = sorted(A.basis())
                sage: S = A.semisimple_quotient(); S
                Semisimple quotient of An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: S in Algebras(QQ).Semisimple()
                True
                sage: S.basis()
                Finite family {'y': B['y'], 'x': B['x']}
                sage: xs,ys = sorted(S.basis())
                sage: (xs + ys) * xs
                B['x']

            Sanity check: the semisimple quotient of the `n`-th
            descent algebra of the symmetric group is of dimension the
            number of partitions of `n`::

                sage: [ DescentAlgebra(QQ,n).B().semisimple_quotient().dimension()
                ....:   for n in range(6) ]
                [1, 1, 2, 3, 5, 7]
                sage: [Partitions(n).cardinality() for n in range(10)]
                [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]

            .. TODO::

               - Pickling by construction, as ``A.semisimple_quotient()``?
               - Lazy evaluation of ``_repr_``

            TESTS::

                sage: TestSuite(S).run()
            """
            ring = self.base_ring()
            category = Algebras(ring).WithBasis().FiniteDimensional().Quotients().Semisimple()
            result = self.quotient_module(self.radical(), category=category)
            result.rename("Semisimple quotient of {}".format(self))
            return result


        @cached_method
        def center_basis(self):
            r"""
            Return a basis of the center of ``self``.

            OUTPUT:

            - a list of elements of ``self``.

            .. SEEALSO:: :meth:`center`

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
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

            The center of a semisimple algebra is semisimple::

                sage: DihedralGroup(6).algebra(QQ).center() in Algebras(QQ).Semisimple()
                True

            .. TODO::

                - Pickling by construction, as ``A.center()``?
                - Lazy evaluation of ``_repr_``

            TESTS::

                sage: TestSuite(center).run()
            """
            category = Algebras(self.base_ring()).FiniteDimensional().Subobjects().Commutative().WithBasis()
            if self in Algebras.Semisimple:
                category = category.Semisimple()
            center = self.submodule(self.center_basis(),
                                    category=category,
                                    already_echelonized=True)
            center.rename("Center of {}".format(self))
            return center

        def principal_ideal(self, a, side='left'):
            r"""
            Construct the ``side`` `A`-module generated by ``a``.

            EXAMPLE::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: A.principal_ideal(A.an_element())
                Free module generated by {0, 1, 2, 3} over Rational Field

            """
            B = self.basis()
            if side == 'left':
                phi = self.module_morphism(on_basis=lambda i: a*B[i],
                        codomain=self,
                        triangular='lower')
            elif side == 'right':
                phi = self.module_morphism(on_basis=lambda i: B[i]*a,
                        codomain=self,
                        triangular='lower')
            else:
                raise Exception("Side must be ``left`` or ``right``")
            ideal = phi.matrix().image()
            #now let's make it a submodule of A
            return self.submodule([self.from_vector(v) for v in ideal.basis()],
                                  already_echelonized=True)

        def central_orthogonal_idempotents(self):
            r"""
            Return a maximal family of central orthogonal idempotents of
            ``self``.

            INPUT:

            - ``self`` -- a finite dimensional algebra

            EXAMPLES::

                sage: import itertools
                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example(); A
                An example of a finite dimensional algebra with basis: the path
                algebra of the Kronecker quiver (containing the arrows a:x->y
                and b:x->y) over Rational Field
                sage: A.central_orthogonal_idempotents()
                [y, x]
                sage: Z12 = Monoids().Finite().example(); Z12
                An example of a finite multiplicative monoid: the integers
                modulo 12
                sage: A = Z12.algebra(QQ)
                sage: orth = A.central_orthogonal_idempotents(); orth
                [-1/2*B[8] + 1/2*B[4], 1/4*B[1] - 1/2*B[4] - 1/4*B[5] +
                1/4*B[7] + 1/2*B[8] - 1/4*B[11], 1/4*B[1] + 1/2*B[3] + 1/4*B[5]
                - 1/4*B[7] - 1/2*B[9] - 1/4*B[11], 1/2*B[9] - 1/2*B[3], B[0],
                -B[0] + 1/2*B[8] + 1/2*B[4], 1/4*B[1] - 1/4*B[5] - 1/4*B[7] +
                1/4*B[11], -B[0] + 1/2*B[9] + 1/2*B[3], B[0] + 1/4*B[1] -
                1/2*B[3] - 1/2*B[4] + 1/4*B[5] + 1/4*B[7] - 1/2*B[8] - 1/2*B[9]
                + 1/4*B[11]]
                sage: all(e*f == f*e for e,f in itertools.product(orth, orth))
                True

            We construct the minimum orthogonal idempotents of the `0`-Hecke
            monoid algebra::

                sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
                sage: W = WeylGroup(['A', 3]); W.rename("W")
                sage: ambient_monoid = FiniteSetMaps(W, action="right")
                sage: pi = W.simple_projections(length_increasing=True).map(ambient_monoid)
                sage: M = AutomaticSemigroup(pi, one=ambient_monoid.one()); M
                A submonoid of (Maps from W to itself) with 3 generators
                sage: A = M.algebra(QQ)
                sage: orth = A.central_orthogonal_idempotents()
                sage: all(e*f == f*e for e,f in itertools.product(orth, orth))
                True

            TODO: More documentation, explanation algorithm, merge with
            _lift_idempotent ?
            """
            Aquo = self.semisimple_quotient()
            orths = [self._lift_idempotent(x) for x in Aquo.central_orthogonal_idempotents()]
            # Construction of the orthogonal idempotents
            idems = []
            one = self.one()
            for g in orths:
                f = sum(idems) # 0 if idems is empty
                idems.append((one - f) * g * (one - f))
            return idems

        def _lift_idempotent(self, x):
            r"""
            Lift an idempotent of the semisimple quotient of ``self`` into an
            idempotent of ``self``.

            ALGORITHM:

            Let '\overline{e}' be a central orthogonal idempotent of the
            semisimple quotient `\overline{A}` of a finite dimensional algebra
            `A`. Let `e \in A` such that the projection of `e` in
            `\overline{A}` is `\overline{e}`. We iterate the formula
            `1 - (1 - e^2)^2` until having an idempotent. See [CR62] for
            correctness and termination proofs.

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: Aquo = A.semisimple_quotient()
                sage: orth = Aquo.central_orthogonal_idempotents()
                sage: A._lift_idempotent(orth[1])
                x

            REFERENCES:

            ..  [CR62] Curtis, Charles W.; Reiner, Irving
                "Representation theory of finite groups and associative
                algebras."
                Pure and Applied Mathematics, Vol. XI Interscience Publishers, a
                division of John Wiley & Sons, New York-London 1962
                pp 545--546

            TODO: merge in central_orthogonal_idempotents?
            """
            idempOld = None
            assert x in self.semisimple_quotient()
            assert (x == x*x), "Element is not an idempotent"
            idemp = x.lift()
            one = idemp.parent().one()
            while idemp != idempOld:
                tmp = idemp
                idemp = (one - (one - idemp**2)**2)
                idempOld = tmp
            return idemp

        @cached_method
        def cartan_invariants_matrix(self, side='left'):
            r"""
            Return the Cartan invariants matrix of the algebra ``self``.

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: A.cartan_invariants_matrix()
                [1 0]
                [2 1]
                sage: A3 = SymmetricGroup(3).algebra(QQ)
                sage: A3.cartan_invariants_matrix()
                [1 0 0]
                [0 1 0]
                [0 0 1]
                sage: Z12 = Monoids().Finite().example()
                sage: A = Z12.algebra(QQ)
                sage: A.cartan_invariants_matrix()
                [1 0 0 0 0 0 0 0 0]
                [0 2 0 0 0 0 0 0 0]
                [0 0 1 0 0 0 0 0 0]
                [0 0 0 1 0 0 0 0 0]
                [0 0 0 0 1 0 0 0 0]
                [0 0 0 0 0 1 0 0 0]
                [0 0 0 0 0 0 1 0 0]
                [0 0 0 0 0 0 0 2 0]
                [0 0 0 0 0 0 0 0 2]

            With the 0-Hecke monoid algebra::

                sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
                sage: W = WeylGroup(['A', 3]); W.rename("W")
                sage: ambient_monoid = FiniteSetMaps(W, action="right")
                sage: pi = W.simple_projections(length_increasing=True).map(ambient_monoid)
                sage: M = AutomaticSemigroup(pi, one=ambient_monoid.one()); M
                A submonoid of (Maps from W to itself) with 3 generators
                sage: A = M.algebra(QQ)
                sage: A.cartan_invariants_matrix()
                [1 0 0 0 0 0 0 0]
                [0 1 1 1 0 0 0 0]
                [0 1 1 1 0 0 0 0]
                [0 1 1 2 0 1 0 0]
                [0 0 0 0 1 1 1 0]
                [0 0 0 1 1 2 1 0]
                [0 0 0 0 1 1 1 0]
                [0 0 0 0 0 0 0 1]
            """
            Aquo = self.semisimple_quotient()
            orth_quo = Aquo.central_orthogonal_idempotents()
            # Dimension of simple modules
            dimSimples = [sqrt(Aquo.principal_ideal(e, side).dimension()) for e in
                    orth_quo]
            orth = [self._lift_idempotent(x) for x in orth_quo]
            return Matrix(self.base_ring(),
                    len(orth),
                    lambda i,j: self.pierce_decomposition_component(orth[i],
                        orth[j]).dimension()/(dimSimples[i]*dimSimples[j]))

        def projective_indecomposables(self, side='left'):
            r"""
            Return the list of indecomposable projective ``side``
            ``self``-modules.

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: projs = A.projective_indecomposables(); projs
                [Free module generated by {0} over Rational Field,
                 Free module generated by {0, 1, 2} over Rational Field]

            We check that the sum of the dimensions of the indecomposable
            projective module is the dimension of ``self``:

                sage: sum([P.dimension() for P in projs]) == A.dimension()
                True

            """
            return [self.principal_ideal(e, side) for e in self.central_orthogonal_idempotents()]

        def pierce_decomposition_component(self, ei, ej):
            r"""
            Return the Pierce decomposition component of ``self``.

            INPUT:

            - ``self`` -- an algebra `A`
            - `e_i` -- an orthogonal idempotent of ``self``
            - `e_j` -- an orthogonal idempotent of ``self``

            OUTPUT:

            - `e_i A e_j` as a subalgebra of ``self``

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: idemp = A.central_orthogonal_idempotents()
                sage: A.pierce_decomposition_component(idemp[1], idemp[0])
                Free module generated by {0, 1} over Rational Field
                sage: A.pierce_decomposition_component(idemp[0], idemp[1])
                Free module generated by {} over Rational Field


            We recover the unique 2 dimensional representations of S4::

                sage: A4 = SymmetricGroup(4).algebra(QQ)
                sage: e = A4.central_orthogonal_idempotents()[2]
                sage: A4.pierce_decomposition_component(e, e)
                Free module generated by {0, 1, 2, 3} over Rational Field
            """
            B = self.basis()
            phi = self.module_morphism(on_basis=lambda k: ei * B[k] * ej,
                                       codomain=self, triangular='lower')
            ideal = phi.matrix().image()
            return self.submodule([self.from_vector(v) for v in ideal.basis()],
                                  already_echelonized=True)

        def pierce_decomposition(self):
            r"""
            Return the Pierce decomposition of ``self``

            OUTPUT:

            - list of subalgebras `e_i A e_j` for `e_i, e_j` orthogonal
              idempotents

              EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: A.pierce_decomposition()
                [Free module generated by {0} over Rational Field, Free module
                generated by {} over Rational Field, Free module generated by
                {0, 1} over Rational Field, Free module generated by {0} over
                Rational Field]

            """
            import itertools
            decomp = []
            for (ei, ej) in itertools.product(self.central_orthogonal_idempotents(),
                    repeat=2):
                decomp.append(self.pierce_decomposition_component(ei, ej))
            return decomp


    class ElementMethods:

        def to_matrix(self, base_ring=None, action=operator.mul, side='left'):
            """
            Return the matrix of the action of ``self`` on the algebra.

            INPUT:

            - ``base_ring`` -- the base ring for the matrix to be constructed
            - ``action`` -- a bivariate function (default: :func:`operator.mul`)
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
                sage: a.to_matrix(side='right')
                [0 0 1 0 0 0]
                [0 0 0 1 0 0]
                [1 0 0 0 0 0]
                [0 1 0 0 0 0]
                [0 0 0 0 0 1]
                [0 0 0 0 1 0]
                sage: a.to_matrix(base_ring=RDF, side="left")
                [0.0 0.0 1.0 0.0 0.0 0.0]
                [0.0 0.0 0.0 0.0 1.0 0.0]
                [1.0 0.0 0.0 0.0 0.0 0.0]
                [0.0 0.0 0.0 0.0 0.0 1.0]
                [0.0 1.0 0.0 0.0 0.0 0.0]
                [0.0 0.0 0.0 1.0 0.0 0.0]

            AUTHORS: Mike Hansen, ...
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
