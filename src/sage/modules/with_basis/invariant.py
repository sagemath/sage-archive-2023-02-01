r"""
Invariant modules
"""

# ****************************************************************************
#       Copyright (C) 2021 Trevor K. Karn <karnx018 at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import operator
from sage.modules.with_basis.subquotient import SubmoduleWithBasis
from sage.categories.finitely_generated_semigroups import FinitelyGeneratedSemigroups
from sage.categories.finite_dimensional_modules_with_basis import FiniteDimensionalModulesWithBasis
from sage.categories.groups import Groups
from sage.sets.family import Family

class FiniteDimensionalInvariantModule(SubmoduleWithBasis):
    r"""
    The invariant submodule under a semigroup action.

    When a semigroup `S` acts on a module `M`, the invariant module is the
    set of elements `m \in M` such that `s \cdot m = m` for all `s \in S`:

    .. MATH::

        M^S := \{m \in M : s \cdot m = m,\, \forall s \in S \}.

    INPUT:

    - ``R`` -- an instance of a ``Representation`` of a semigroup `S`
      acting on the module `M`

    EXAMPLES:

    First, we create the invariant defined by the cyclic group action on the
    free module with basis `\{1,2,3\}`::

        sage: G = CyclicPermutationGroup(3)
        sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
        sage: action = lambda g, m: M.monomial(g(m))  # cyclically permute coordinates

    In order to give the module an action of ``G``, we create a
    :class:`~sage.modules.with_basis.representation.Representation`::

        sage: from sage.modules.with_basis.representation import Representation
        sage: R = Representation(G, M, action)
        sage: I = R.invariant_module()

    Then we can lift the basis from the invariant to the original module::

        sage: [I.lift(b) for b in I.basis()]
        [M[1] + M[2] + M[3]]
    
    The we could also have the action be a right-action, instead of the
    default left-action::

        sage: def rt_action(g, m): return M.monomial(g(m))  # cyclically permute coordinates
        sage: R = Representation(G, M, rt_action, side='right')  # same as last but on right
        sage: g = G.an_element(); g
        (1,2,3)
        sage: r = R.an_element(); r
        2*M[1] + 2*M[2] + 3*M[3]
        sage: R.side()
        'right'

    So now we can see that multiplication with ``g`` on the right sends
    ``M[1]`` to ``M[2]`` and so on::

        sage: r * g
        3*M[1] + 2*M[2] + 2*M[3]
        sage: I = R.invariant_module()
        sage: [I.lift(b) for b in I.basis()]
        [M[1] + M[2] + M[3]]

    Now we will take the regular representation of the symmetric group on
    three elements to be the module, and compute its invariant submodule::

        sage: G = SymmetricGroup(3)
        sage: R = G.regular_representation(QQ)
        sage: I = R.invariant_module()
        sage: [I.lift(b).to_vector() for b in I.basis()]
        [(1, 1, 1, 1, 1, 1)]

    We can also check the scalar multiplication by elements of the base ring
    (for this example, the rational field)::

        sage: [I.lift(3*b).to_vector() for b in I.basis()]
        [(3, 3, 3, 3, 3, 3)]

    A more subtle example is the invariant submodule of a skew-commutative
    module, for example the exterior algebra `E[x_0,x_1,x_2]` generated
    by three elements::

        sage: G = CyclicPermutationGroup(3)
        sage: M = algebras.Exterior(QQ, 'x', 3)
        sage: def cyclic_ext_action(g, m):
        ....:     # cyclically permute generators
        ....:     return M.prod([M.monomial((g(j+1)-1,)) for j in m])

    If you care about being able to exploit the algebra structure of the
    exterior algebra (i.e. if you want to multiply elements together), you
    should make sure the representation knows it is also an algebra with
    the semigroup action being by algebra endomorphisms::

        sage: cat = Algebras(QQ).WithBasis().FiniteDimensional()
        sage: R = Representation(G, M, cyclic_ext_action, category=cat)
        sage: I = R.invariant_module()

    We can express the basis in the ambient algebra (`E[x_0,x_1,x_2]`)::

        sage: [I.lift(b) for b in I.basis()]
        [1, x0 + x1 + x2, x0*x1 - x0*x2 + x1*x2, x0*x1*x2]

    or we can express the basis intrinsicallly to the invariant ``I``::

        sage: B = I.basis()
        sage: m = 3*B[0] + 2*B[1] + 7*B[3]

    This lifts to the exterior algebra::

        sage: I.lift(m)
        3 + 2*x0 + 7*x0*x1*x2 + 2*x1 + 2*x2

    We can also check using the invariant element ``m`` that arithmetic works::

        sage: m^2
        9*B[0] + 12*B[1] + 42*B[3]
        sage: m+m
        6*B[0] + 4*B[1] + 14*B[3]

    To see the actual elements expressed in the exterior algebra, we lift them
    again::

        sage: I.lift(m+m)
        6 + 4*x0 + 14*x0*x1*x2 + 4*x1 + 4*x2
        sage: 7*m
        21*B[0] + 14*B[1] + 49*B[3]
        sage: I.lift(7*m)
        21 + 14*x0 + 49*x0*x1*x2 + 14*x1 + 14*x2

    The classic example of an invariant module is the module of symmetric
    functions, which is the invariant module of polynomials whose variables
    are acted upon by permutation. We can create a module isomorphic to the
    homogeneous component of a a polynomial ring in `n` variable of a fixed
    degree `d` by looking at  weak  compositions of `d` of length `n`, which
    we consider as the exponent vector. For example, `x^2yz \in \QQ[x,y,z]`
    would have the exponent vector `(2,1,1)`. The vector `(2,1,1)` is a
    weak composition of `4`, with length `3`, and so we can think of it as
    being in the degree-`4` homogeneous component of a polynomial ring
    in three variables::

        sage: C = IntegerVectors(4, length=3, min_part=0)  # representing degree-4 monomials
        sage: M = CombinatorialFreeModule(QQ, C)  # isomorphic to deg-4 homog. polynomials
        sage: G = SymmetricGroup(3)
        sage: def perm_action(g,x): return M.monomial(C(g(list(x))))
        sage: perm_action(G((1,2,3)), C([4,3,2]))
        B[[3, 2, 4]]
        sage: R = Representation(G, M, perm_action)
        sage: I = R.invariant_module()
        sage: [I.lift(b) for b in I.basis()]
        [B[[0, 0, 4]] + B[[0, 4, 0]] + B[[4, 0, 0]],
         B[[0, 1, 3]] + B[[0, 3, 1]] + B[[1, 0, 3]]
         + B[[1, 3, 0]] + B[[3, 0, 1]] + B[[3, 1, 0]],
         B[[0, 2, 2]] + B[[2, 0, 2]] + B[[2, 2, 0]],
         B[[1, 1, 2]] + B[[1, 2, 1]] + B[[2, 1, 1]]]

    These are the monomial symmetric functions, which are a well-known
    basis for the symmetric functions. For comparison::

        sage: Sym = SymmetricFunctions(QQ)
        sage: m = Sym.monomial()
        sage: [m[mu].expand(3) for mu in Partitions(4)]
        [x0^4 + x1^4 + x2^4,
         x0^3*x1 + x0*x1^3 + x0^3*x2 + x1^3*x2 + x0*x2^3 + x1*x2^3,
         x0^2*x1^2 + x0^2*x2^2 + x1^2*x2^2,
         x0^2*x1*x2 + x0*x1^2*x2 + x0*x1*x2^2,
         0]

    .. NOTE::

        The current implementation works when `S` is a finitely-generated
        semigroup, and when `M` is a finite-dimensional free module with
        a distinguished basis.

    .. TODO::

        Extend this to have multiple actions, including actions on both sides.

    .. TODO::

        Extend when `M` does not have a basis and `S` is a permutation
        group using:

        - https://arxiv.org/abs/0812.3082
        - https://www.dmtcs.org/pdfpapers/dmAA0123.pdf
    """
    def __init__(self, M, S, action=operator.mul, side='left', *args, **kwargs):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(3)
            sage: R = G.regular_representation()
            sage: I = R.invariant_module()
            sage: TestSuite(I).run()

        TESTS::

            sage: G = GroupExp()(QQ) # a group that is not finitely generated
            sage: M = CombinatorialFreeModule(QQ, [1,2,3])
            sage: def on_basis(g,m): return M.monomial(m)  # trivial rep'n
            sage: from sage.modules.with_basis.representation import Representation
            sage: R = Representation(G, M, on_basis)
            sage: R.invariant_module()
            Traceback (most recent call last):
            ...
            ValueError: Multiplicative form of Rational Field is not finitely generated
        """
        if S not in FinitelyGeneratedSemigroups():
            raise ValueError(f"{S} is not finitely generated")
        if M not in FiniteDimensionalModulesWithBasis:
            raise ValueError(f"{M} is not a finite dimensional module with a distinguished basis")

        if side == "left":
            def _invariant_map(g, x):
                return action(g, x) - x
        elif side == "right":
            def _invariant_map(g, x):
                return action(x, g) - x
        else:
            raise ValueError("side must either be 'left' or 'right'")

        self._side = side
        self._action = action
        self._semigroup = S

        category = kwargs.pop("category", M.category().Subobjects())

        # Give the intersection of kernels of the map `s*x-x` to determine when
        # `s*x = x` for all generators `s` of `S`
        basis = M.annihilator_basis(S.gens(), action=_invariant_map, side="left")

        super().__init__(Family(basis),
                         support_order=M._compute_support_order(basis),
                         ambient=M,
                         unitriangular=False,
                         category=category,
                         *args, **kwargs)

    def _repr_(self):
        r"""
        Return a string representaion of ``self``.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(3)
            sage: R = G.trivial_representation()
            sage: R.invariant_module()
            (Cyclic group of order 3 as a permutation group)-invariant submodule of
             Trivial representation of Cyclic group of order 3 as a permutation group over Integer Ring

            sage: G = CyclicPermutationGroup(3)
            sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
            sage: action = lambda g, m: M.monomial(g(m))  # cyclically permute coordinates
            sage: M.invariant_module(G, action_on_basis=action)
            (Cyclic group of order 3 as a permutation group)-invariant submodule of
             Free module generated by {1, 2, 3} over Rational Field
        """
        M = self._ambient
        from sage.modules.with_basis.representation import Representation
        if isinstance(self._ambient, Representation):
            M = M._module
        return f"({self._semigroup})-invariant submodule of {M}"

    def _latex_(self):
        r"""
        Return a latex representaion of ``self``.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(3)
            sage: R = G.algebra(QQ)
            sage: latex(R.invariant_module(G))
            \left( \Bold{Q}[\langle (1,2,3) \rangle] \right)^{\langle (1,2,3) \rangle}
        """
        M = self._ambient
        from sage.modules.with_basis.representation import Representation
        if isinstance(self._ambient, Representation):
            M = M._module
        from sage.misc.latex import latex
        return "\\left( {} \\right)^{{{}}}".format(latex(M), latex(self._semigroup))

    def _test_invariant(self, **options): ## Lift to representation and check that the element is invariant
        """
        Check (on some elements) that ``self`` is invariant.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
            sage: def action(g, x): return M.monomial(g(x))
            sage: I = M.invariant_module(G, action_on_basis=action)
            sage: I._test_invariant()

            sage: G = SymmetricGroup(10)
            sage: M = CombinatorialFreeModule(QQ, list(range(1,11)), prefix='M')
            sage: def action(g, x): return M.monomial(g(x))
            sage: I = M.invariant_module(G, action_on_basis=action)
            sage: I._test_invariant(max_runs=10)
        """
        tester = self._tester(**options)
        X = tester.some_elements()
        L = []
        max_len = tester._max_runs

        # FIXME: This is max_len * dim number of runs!!!
        for i, x in enumerate(self._semigroup):
            L.append(x)
            if i >= max_len:
                break

        for x in L:
            for elt in X:
                lifted = self.lift(elt)
                if self._side == 'left':
                    tester.assertEqual(self._action(x, lifted), lifted)
                else:
                    tester.assertEqual(self._action(lifted, x), lifted)

    def semigroup(self):
        r"""
        Return the semigroup `S` whose action ``self`` is invariant under.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
            sage: def action(g,x): return M.monomial(g(x))
            sage: I = M.invariant_module(G, action_on_basis=action)
            sage: I.semigroup()
            Symmetric group of order 3! as a permutation group
        """
        return self._semigroup

    semigroup_representation = SubmoduleWithBasis.ambient

    class Element(SubmoduleWithBasis.Element):
        def _mul_(self, other):
            r"""
            Multiply ``self`` and ``other``.

            EXAMPLES:

            In general, there is not a well defined multiplication between
            two elements of a given module, but there is a multiplication
            with scalars::

                sage: M = CombinatorialFreeModule(QQ,[1,2,3],prefix='M');
                sage: G = CyclicPermutationGroup(3); G.rename('G')
                sage: g = G.an_element(); g
                (1,2,3)
                sage: from sage.modules.with_basis.representation import Representation
                sage: R = Representation(G,M,lambda g,x:M.monomial(g(x))); R.rename('R')
                sage: I = R.invariant_module()
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [M[1] + M[2] + M[3]]
                sage: v = B[0]
                sage: v*v
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: 'R' and 'R'
                sage: (1/2) * v
                1/2*B[0]
                sage: v * (1/2)
                1/2*B[0]
                sage: R.rename()  # reset name

            Sometimes, the module is also a ring. To ensure the multiplication
            works as desired, we should be sure to pass the correct category to
            the :class:`~sage.modules.with_basis.representation.Representation`.
            In the following example, we use the exterior algebra over `\QQ`
            with three generators, which is in the category of finite
            dimensional `\QQ`-algebras with a basis::

                sage: G = CyclicPermutationGroup(3); G.rename('G')
                sage: M = algebras.Exterior(QQ, 'x', 3)
                sage: def on_basis(g,m): return M.prod([M.monomial((g(j+1)-1,)) for j in m])  # cyclically permute generators
                sage: R = Representation(G, M, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional(), side='right')
                sage: I = R.invariant_module(); I.rename('I')
                sage: B = I.basis()
                sage: v = B[0] + 2*B[1]; I.lift(v)
                1 + 2*x0 + 2*x1 + 2*x2
                sage: w = B[2]; I.lift(w)
                x0*x1 - x0*x2 + x1*x2
                sage: v * w
                B[2] + 6*B[3]
                sage: I.lift(v*w)
                x0*x1 + 6*x0*x1*x2 - x0*x2 + x1*x2
                sage: w * v
                B[2] + 6*B[3]
                sage: (1/2) * v
                1/2*B[0] + B[1]
                sage: w * (1/2)
                1/2*B[2]
                sage: g = G((1,3,2))
                sage: v * g
                B[0] + 2*B[1]
                sage: w * g
                B[2]
                sage: g * v
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: 'G' and 'I'
                sage: I.rename()  # reset name

                sage: R = Representation(G, M, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional())
                sage: I = R.invariant_module(); I.rename('I')
                sage: B = I.basis()
                sage: v = B[0] + 2*B[1]; I.lift(v)
                1 + 2*x0 + 2*x1 + 2*x2
                sage: w = B[2]; I.lift(w)
                x0*x1 - x0*x2 + x1*x2
                sage: v * w
                B[2] + 6*B[3]
                sage: I.lift(v*w)
                x0*x1 + 6*x0*x1*x2 - x0*x2 + x1*x2
                sage: w * v
                B[2] + 6*B[3]
                sage: (1/2) * v
                1/2*B[0] + B[1]
                sage: w * (1/2)
                1/2*B[2]
                sage: g = G((1,3,2))
                sage: v * v
                B[0] + 4*B[1]
                sage: g * w
                B[2]
                sage: v * g
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: 'I' and 'G'
                sage: G.rename(); I.rename()  # reset names
            """
            P = self.parent()
            return P.retract(P.lift(self) * P.lift(other))

        def _acted_upon_(self, scalar, self_on_left=False):
            """
            EXAMPLES::

                sage: G = CyclicPermutationGroup(3)
                sage: g = G.an_element(); g
                (1,2,3)
                sage: M = CombinatorialFreeModule(QQ, [1,2,3])
                sage: E = algebras.Exterior(QQ, 'x', 3)
                sage: from sage.modules.with_basis.representation import Representation
                sage: R = Representation(G, M, lambda g,x: M.monomial(g(x)))
                sage: I = R.invariant_module()
                sage: [b._acted_upon_(G((1,3,2))) for b in I.basis()]
                [B[0]]
                sage: v = I.an_element(); v
                2*B[0]
                sage: g * v
                2*B[0]
                sage: [g * v for g in G.list()]
                [2*B[0], 2*B[0], 2*B[0]]


                sage: def on_basis(g,m): return E.prod([E.monomial((g(j+1)-1,)) for j in m])  # cyclically permute generators
                sage: R = Representation(G, E, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional())
                sage: I = R.invariant_module()
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [1, x0 + x1 + x2, x0*x1 - x0*x2 + x1*x2, x0*x1*x2]
                sage: [[g*b for g in G] for b in B]
                [[B[0], B[0], B[0]],
                 [B[1], B[1], B[1]],
                 [B[2], B[2], B[2]],
                 [B[3], B[3], B[3]]]
                sage: 3 * I.basis()[0]
                3*B[0]
                sage: 3*B[0] + B[1]*2
                3*B[0] + 2*B[1]

                sage: R = G.regular_representation(QQ)
                sage: I = R.invariant_module()
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [() + (1,2,3) + (1,3,2)]
                sage: B[0]._acted_upon_(G((1,3,2)))
                B[0]
                sage: B[0]._acted_upon_(G((1,3,2)), self_on_left=True) is None
                True

                sage: R = G.regular_representation(QQ, side='right')
                sage: I = R.invariant_module()
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [() + (1,2,3) + (1,3,2)]
                sage: g = G((1,3,2))
                sage: B[0]._acted_upon_(g, self_on_left=True)
                B[0]
                sage: B[0]._acted_upon_(g, self_on_left=False) is None
                True

                sage: R = Representation(G, M, lambda g,x: M.monomial(g(x)), side='right')
                sage: I = R.invariant_module()
                sage: v = I.an_element(); v
                2*B[0]
                sage: v * g
                2*B[0]
                sage: [v * g for g in G.list()]
                [2*B[0], 2*B[0], 2*B[0]]
                sage: [b._acted_upon_(G((1,3,2)), self_on_left=True) for b in I.basis()]
                [B[0]]

                sage: def on_basis(g,m): return E.prod([E.monomial((g(j+1)-1,)) for j in m])  # cyclically permute generators
                sage: R = Representation(G, E, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional(), side='right')
                sage: I = R.invariant_module()
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [1, x0 + x1 + x2, x0*x1 - x0*x2 + x1*x2, x0*x1*x2]
                sage: [[b * g for g in G] for b in B]
                [[B[0], B[0], B[0]],
                 [B[1], B[1], B[1]],
                 [B[2], B[2], B[2]],
                 [B[3], B[3], B[3]]]
                sage: 3 * B[0] + B[1] * 2
                3*B[0] + 2*B[1]
            """
            if scalar in self.parent()._semigroup and self_on_left == (self.parent()._side == 'right'):
                return self
            return super()._acted_upon_(scalar, self_on_left)

