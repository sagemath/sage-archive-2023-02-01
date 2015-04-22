# -*- coding: utf-8 -*-
r"""
Finite dimensional algebras with basis


REFERENCES:

..  [CR62] Curtis, Charles W.; Reiner, Irving
    "Representation theory of finite groups and associative
    algebras."
    Pure and Applied Mathematics, Vol. XI Interscience Publishers, a
    division of John Wiley & Sons, New York-London 1962
    pp 545--547
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
from sage.misc.cachefunc import cached_method
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
        def radical_basis(self):
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
                sage: A.radical_basis()
                [B[()] + B[(1,2)(3,4)], B[(3,4)] + B[(1,2)(3,4)], B[(1,2)] + B[(1,2)(3,4)]]

                sage: A = KleinFourGroup().algebra(QQ, category=Monoids())
                sage: A.radical_basis.__module__
                'sage.categories.finite_dimensional_algebras_with_basis'
                sage: A.radical_basis()
                []
            """
            F = self.base_ring()
            if not F.is_field():
                raise NotImplementedError("the base ring must be a field")
            p = F.characteristic()
            from sage.matrix.constructor import matrix
            from sage.modules.free_module_element import vector

            product_on_basis = self.product_on_basis

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
            Construct the ``side`` principal ideal generated by ``a``.

            EXAMPLES:

            In order to highlight the difference between left and
            right principal ideals, our first example deals with a non
            commutative algebra::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: x, y, a, b = A.basis()

            In this algebra, multiplication on the right by `x`
            annihilates all basis elements but `x`::

                sage: x*x, y*x, a*x, b*x
                (x, 0, 0, 0)

            so the left ideal generated by `x` is one-dimensional::

                sage: Ax = A.principal_ideal(x, side='left'); Ax
                Free module generated by {0} over Rational Field
                sage: [B.lift() for B in Ax.basis()]
                [x]

            Multiplication on the left by `x` annihilates
            only `x` and fixes the other basis elements::

                sage: x*x, x*y, x*a, x*b
                (x, 0, a, b)

            so the right ideal generated by `x` is 3-dimensional::

                sage: xA = A.principal_ideal(x, side='right'); xA
                Free module generated by {0, 1, 2} over Rational Field
                sage: [B.lift() for B in xA.basis()]
                [x, a, b]
            """
            return self.submodule([(a * b if side=='right' else b * a)
                                    for b in self.basis()])

        @cached_method
        def orthogonal_idempotents_central_mod_rad(self):
            r"""
            Return a family of orthogonal idempotents of ``self`` that project
            on the central orthogonal idempotents of the semisimple quotient.

            INPUT:

            - ``self`` -- a finite dimensional algebra

            OUTPUT:

            - a list of orthogonal idempotents obtained by lifting the central
              orthogonal idempotents of the semisimple quotient.

            ALGORITHM:

            The orthogonal idempotents of `A` are obtained by lifting the
            central orthogonal idempotents of the semisimple quotient
            `\overline{A}`.

            Let `(\overline{e_i})` be a set of central orthogonal idempotents
            of the semisimple quotient of `A` and `(e_i)` their lift in `A`.
            We recursively construct orthogonal idempotents of `A`: if
            `(f_i)_{i < n}` is a set of already constructed orthogonal
            idempotent, we then construct `f_n` by lifting the element
            `(1 - \sum_{i < n} f_i) e_n (1 - \sum_{i < n} f_i)`.

            See [CR62] for correctness and termination proofs.

            .. SEEALSO:: :meth:`idempotent_lift`


            EXAMPLES::

                sage: import itertools
                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example(); A
                An example of a finite dimensional algebra with basis:
                the path algebra of the Kronecker quiver
                (containing the arrows a:x->y and b:x->y) over Rational Field
                sage: A.orthogonal_idempotents_central_mod_rad()
                [y, x]

            ::

                sage: Z12 = Monoids().Finite().example(); Z12
                An example of a finite multiplicative monoid: the integers modulo 12
                sage: A = Z12.algebra(QQ)
                sage: orth = A.orthogonal_idempotents_central_mod_rad(); orth
                [-1/2*B[8] + 1/2*B[4],
                 1/4*B[1] - 1/2*B[4] - 1/4*B[5] + 1/4*B[7] + 1/2*B[8] - 1/4*B[11],
                 1/4*B[1] + 1/2*B[3] + 1/4*B[5] - 1/4*B[7] - 1/2*B[9] - 1/4*B[11],
                 -1/2*B[3] + 1/2*B[9],
                 B[0],
                 1/2*B[8] + 1/2*B[4] - B[0],
                 1/4*B[1] + 1/4*B[11] - 1/4*B[5] - 1/4*B[7],
                 -B[0] + 1/2*B[3] + 1/2*B[9],
                 B[0] + 1/4*B[1] - 1/2*B[3] - 1/2*B[4] + 1/4*B[5] + 1/4*B[7] - 1/2*B[8] - 1/2*B[9] + 1/4*B[11]]
                sage: all(e*e == e for e in orth)
                True
                sage: all(e*f == f*e and e*f == 0 for e,f in itertools.product(orth, orth) if e!= f)
                True

            We construct the minimal orthogonal idempotents of the `0`-Hecke
            monoid algebra::

                sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
                sage: W = WeylGroup(['A', 3]); W.rename("W")
                sage: ambient_monoid = FiniteSetMaps(W, action="right")
                sage: pi = W.simple_projections(length_increasing=True).map(ambient_monoid)
                sage: M = AutomaticSemigroup(pi, one=ambient_monoid.one()); M
                A submonoid of (Maps from W to itself) with 3 generators
                sage: A = M.algebra(QQ)
                sage: orth = A.orthogonal_idempotents_central_mod_rad()
                sage: all(e*e == e for e in orth)
                True
                sage: all(e*f == f*e and e*f == 0 for e,f in itertools.product(orth, orth) if e!= f)
                True
            """
            Aquo = self.semisimple_quotient()
            one = self.one()
            # Construction of the orthogonal idempotents
            idems = []
            f = 0
            for g in Aquo.central_orthogonal_idempotents():
                idems.append(self.idempotent_lift((one - f) * g.lift() * (one - f)))
                f = f + idems[-1]
            return idems

        def idempotent_lift(self, x):
            r"""
            Return an idempotent of ``self`` which projection on the semisimple
            quotient is the same as `x`.

            Let `\pi` be the projection `A \rightarrow \overline{A}` on the
            quotient by the radical.

            INPUT:

            1. either `x` -- an idempotent of the semisimple quotient of `A`
            2. or `x` -- an element of an algebra `A` which project on an
               idempotent element of the semisimple quotient of `A`

            OUTPUT:

            1. if `x` is in the semisimple quotient of `A`, return an idempotent
               `e` such that `\pi(e) = x`.
            2. if `x` is in `A`, return the unique idempotent `e` of `A` such
               that `\pi(e) = \pi(x)` and `e` is polynomial in `x`.

            ALGORITHM:

            For `\overline{e}` an idempotent of `\overline{A}`, we construct `e
            \in A` idempotent such that `\pi(e) = \overline{e}`. Namely, we
            find an element of `A` for which the projection on `\overline{A}`
            is `\overline{e}` and then iterate the formula `1 - (1 - e^2)^2`
            until having an idempotent.

            See [CR62] for correctness and termination proofs.

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: S = A.semisimple_quotient()
                sage: A.idempotent_lift(S.basis()['x'])
                x
                sage: A.idempotent_lift(A.basis()['y'])
                y
            """
            if not self.is_parent_of(x):
                x = x.lift()
            p = self.semisimple_quotient().retract(x)
            if p*p != p:
                raise ValueError("%s does not retract to an idempotent."%p)
            xOld = None
            one = self.one()
            while x != xOld:
                tmp = x
                x = (one - (one - x**2)**2)
                xOld = tmp
            return x

        @cached_method
        def cartan_invariants_matrix(self):
            r"""
            Return the Cartan invariants matrix of the algebra.

            OUTPUT:

            - The Cartan invariants matrix of the algebra.

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
            dimSimples = [sqrt(Aquo.principal_ideal(e).dimension()) for e in
                    orth_quo]
            # Orthogonal idempotents
            orth = self.orthogonal_idempotents_central_mod_rad()
            return Matrix(self.base_ring(),
                    len(orth),
                    lambda i,j: self.peirce_summand(orth[i],
                        orth[j]).dimension()/(dimSimples[i]*dimSimples[j]))

        def projective_isotypics(self, side='left'):
            r"""
            Return the list of isotypic projective ``side`` ``self``-modules.

            Let `(e_i)` be the orthogonal idempotents lifted from the central
            orthogonal idempotents of the semisimple algebra. The regular
            representation `A` can be decomposed as a direct sum
            `\bigoplus_i P_i` with `P_i = e_i A`.
            The projective modules `P_i` are isotypic in the sense that all
            their direct summands are isomorphic.

            .. NOTE::

                The number of summands of `P_i` is the dimension of the associated
                simple module.

            OUTPUT:

            - return a list of modules `e_i A` for `(e_i)` a set of orthogonal
              idempotents lifted from central orthogonal idempotents of the
              semisimple quotient.

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: projs = A.projective_isotypics(side="left"); projs
                [Free module generated by {0, 1, 2} over Rational Field,
                 Free module generated by {0} over Rational Field]

            We check that the sum of the dimensions of the indecomposable
            projective module is the dimension of ``self``:

                sage: sum([P.dimension() for P in projs]) == A.dimension()
                True

            """
            return [self.principal_ideal(e, side) for e in
                    self.orthogonal_idempotents_central_mod_rad()]

        @cached_method
        def peirce_summand(self, ei, ej):
            r"""
            Return the Peirce decomposition summand `e_i A e_j` as a
            subspace of the algebra.

            INPUT:

            - ``self`` -- an algebra `A`
            - `e_i` -- an idempotent of `A`
            - `e_j` -- an idempotent of `A`

            OUTPUT:

            - `e_i A e_j` as a subspace of `A`

            .. SEEALSO:: :meth:`peirce_decomposition`

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: idemp = A.orthogonal_idempotents_central_mod_rad()
                sage: A.peirce_summand(idemp[1], idemp[0])
                Free module generated by {0, 1} over Rational Field
                sage: A.peirce_summand(idemp[0], idemp[1])
                Free module generated by {} over Rational Field

            We recover the unique 2 dimensional representation of S4::

                sage: A4 = SymmetricGroup(4).algebra(QQ)
                sage: e = A4.central_orthogonal_idempotents()[2]
                sage: A4.peirce_summand(e, e)
                Free module generated by {0, 1, 2, 3} over Rational Field
            """
            B = self.basis()
            phi = self.module_morphism(on_basis=lambda k: ei * B[k] * ej,
                                       codomain=self, triangular='lower')
            ideal = phi.matrix().image()
            return self.submodule([self.from_vector(v) for v in ideal.basis()],
                                  already_echelonized=True)


        def peirce_decomposition(self, list_idempotents=None, check=True):
            r"""
            Return the Peirce decomposition of ``self``.

            Let `(e_i)` be a set of orthogonal idempotents of `A` with sum `1`,
            then `A` is the direct sum of the subspaces `e_i A e_j`.

            INPUT:

            - ``list_idempotents`` -- (default: None) A list of orthogonal
              idempotents of the algebra. By default, the list of orthogonal
              idempotents will be lifted from the central orthogonal
              idempotents of the semisimple quotient.
            - ``check`` -- (default: True) If set to True ``list_idempotents``
              is checked to be orthogonal.

            OUTPUT:

            - The list of subspaces `e_i A e_j`, where `e_i` and `e_j` run
            through ``list_idempotents``.

            .. SEEALSO:: :meth:`orthogonal_idempotents_central_mod_rad`

            EXAMPLES::

                sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example()
                sage: A.peirce_decomposition()
                [Free module generated by {0} over Rational Field, Free module
                generated by {} over Rational Field, Free module generated by
                {0, 1} over Rational Field, Free module generated by {0} over
                Rational Field]
            """
            import itertools
            if list_idempotents is None:
                list_idempotents = self.orthogonal_idempotents_central_mod_rad()
            if check:
                if not (all(e*e == e for e in list_idempotents) and
                        all(e*f == f*e and e*f == 0 for e, f in
                            itertools.product(list_idempotents, repeat=2) if e != f)):
                    raise ValueError("Not an orthogonal list of idempotents.")
            return [self.peirce_summand(ei, ej) for ei, ej in
                    itertools.product(list_idempotents, repeat=2)]


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
