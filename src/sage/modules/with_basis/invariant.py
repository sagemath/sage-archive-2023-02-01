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

## TODO: COMB THORUGH TO MAKE SURE ALL STUFF IS HAPPENING IN REPN NOT IN MODULE

from sage.modules.with_basis.subquotient import SubmoduleWithBasis
from sage.categories.finitely_generated_semigroups import FinitelyGeneratedSemigroups
from sage.categories.finite_dimensional_modules_with_basis import FiniteDimensionalModulesWithBasis
from sage.categories.groups import Groups
from sage.sets.family import Family

class FiniteDimensionalInvariantModule(SubmoduleWithBasis):
    r"""
    Construct the `S`-invariant submodule of `M`. When a semigroup `S` acts on a module
    `M`, the invariant module is the collection of elements `m` in `M` such that
    `s \cdot m = m` for all `s \in S.

    .. MATH::

        M^S = \{m \in M : s\cdot m = m,\, \forall s \in S \}

    INPUTS:

    - ``R`` -- an instance of a ``Representation`` of a semigroup `S`
               acting on the module `M`.

    OUTPUTS:

    - ``MS`` -- the invariant algebra of the semigroup action of `S` on `M`, or
                equivalently, the isotypic component of the representation of
                `S` carried by `M` corresponding to the trivial character.

    EXAMPLES::

        sage: G = CyclicPermutationGroup(3)
        sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
        sage: from sage.modules.with_basis.representation import Representation
        sage: action = lambda g, m: M.term(g(m)) #cyclically permute coordinates
        sage: R = Representation(G, M, action)
        sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
        sage: I = FiniteDimensionalInvariantModule(R)
        sage: [I.lift(b) for b in I.basis()]
        [M[1] + M[2] + M[3]]

        sage: G = CyclicPermutationGroup(3)
        sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
        sage: from sage.modules.with_basis.representation import Representation
        sage: action = lambda g, m: M.term(g(m)) #cyclically permute coordinates
        sage: R = Representation(G, M, action, side='right') #same as last but on right
        sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
        sage: g = G.an_element(); g
        (1,2,3)
        sage: r = R.an_element(); r
        2*M[1] + 2*M[2] + 3*M[3]
        sage: R.side()
        'right'
        sage: r*g
        3*M[1] + 2*M[2] + 2*M[3]
        sage: I = FiniteDimensionalInvariantModule(R)
        sage: [I.lift(b) for b in I.basis()]
        [M[1] + M[2] + M[3]]

        sage: G = SymmetricGroup(3)
        sage: R = G.regular_representation(QQ)
        sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
        sage: I = FiniteDimensionalInvariantModule(R)
        sage: [I.lift(b).to_vector() for b in I.basis()]
        [(1, 1, 1, 1, 1, 1)]
        sage: [I.lift(3*b).to_vector() for b in I.basis()]
        [(3, 3, 3, 3, 3, 3)]

        sage: G = CyclicPermutationGroup(3)
        sage: M = algebras.Exterior(QQ, 'x', 3)
        sage: from sage.modules.with_basis.representation import Representation
        sage: on_basis = lambda g,m: M.prod([M.monomial((g(j+1)-1,)) for j in m]) #cyclically permute generators
        sage: from sage.categories.algebras import Algebras
        sage: R = Representation(G, M, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional())
        sage: I = FiniteDimensionalInvariantModule(R)
        sage: [I.lift(b) for b in I.basis()]
        [1, x0 + x1 + x2, x0*x1 - x0*x2 + x1*x2, x0*x1*x2]
        sage: B = I.basis()
        sage: m = 3*B[0] + 2*B[1] + 7*B[3]
        sage: I.lift(m)
        3 + 2*x0 + 7*x0*x1*x2 + 2*x1 + 2*x2
        sage: m^2
        9*B[0] + 12*B[1] + 42*B[3]
        sage: m+m
        6*B[0] + 4*B[1] + 14*B[3]
        sage: I.lift(m+m)
        6 + 4*x0 + 14*x0*x1*x2 + 4*x1 + 4*x2
        sage: 7*m
        21*B[0] + 14*B[1] + 49*B[3]
        sage: I.lift(7*m)
        21 + 14*x0 + 49*x0*x1*x2 + 14*x1 + 14*x2

    .. NOTE:: 

        The current implementation works when `S` is a finitely-generated semigroup,
        and when `M` is a finite-dimensional free module with a distinguished basis.

    .. TODO::

        Extend when `M` does not have a basis and `S` is a permutation group using:
        - https://arxiv.org/abs/0812.3082
        - https://www.dmtcs.org/pdfpapers/dmAA0123.pdf

    """

    def __init__(self, R, *args, **kwargs):
        """
        TESTS::
            
            sage: G = GroupExp()(QQ) # a group that is not finitely generated
            sage: M = CombinatorialFreeModule(QQ, [1,2,3])
            sage: on_basis = lambda g,m: M.term(m) # trivial rep'n
            sage: from sage.modules.with_basis.representation import Representation
            sage: R = Representation(G, M, on_basis)
            sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
            sage: I = FiniteDimensionalInvariantModule(R)
            Traceback (most recent call last):
            ...
            ValueError: Multiplicative form of Rational Field is not finitely generated

        """

        self._semigroup_representation = R
        self._semigroup = R.semigroup()

        # A check for self._semigroup_representation._module not in FiniteDimensionalModulesWithBasis
        # is not required, because a ``Representation`` cannot be built without a basis
        if self._semigroup not in FinitelyGeneratedSemigroups:
            raise ValueError(f'{self._semigroup} is not finitely generated')

        # The left/right multiplication is taken care of
        # by self._semigroup_representation, so here
        # we can just pass the left multiplication.
        # This means that the side argument of annihilator_basis
        # (see below) will always be side = 'left'
        if self._semigroup_representation.side() == 'left':
            def _invariant_map(g, x):
                return g*x - x
        elif self._semigroup_representation.side() == 'right':
            def _invariant_map(g, x):
                return x*g - x
        
        self._invariant_map = _invariant_map
        
        category = kwargs.pop('category', R.category().Subobjects())

        # Give the intersection of kernels of the map `s*x-x` to determine when
        # `s*x = x` for all generators `s` of `S`
        basis = self._semigroup_representation.annihilator_basis(
                    self._semigroup.gens(),
                    action = self._invariant_map,
                    side = 'left')

        super().__init__(Family(basis),
                        support_order = self._semigroup_representation._compute_support_order(basis),
                        ambient = self._semigroup_representation,
                        unitriangular = False,#is this right?
                        category = category,
                        *args, **kwargs)

    def _repr_(self):
        """
        EXAMPLES::

            sage: M = CombinatorialFreeModule(QQ,[1,2,3])
            sage: G = CyclicPermutationGroup(3)
            sage: from sage.modules.with_basis.representation import Representation
            sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
            sage: R = Representation(G,M,lambda g,x: M.monomial(g(x)))
            sage: FiniteDimensionalInvariantModule(R)
            (Cyclic group of order 3 as a permutation group)-invariant submodule of
             Free module generated by {1, 2, 3} over Rational Field

        """

        return f"({self._semigroup})-invariant submodule of {self._semigroup_representation._module}"


    def semigroup(self):
        """
        Return the semigroup `S` whose action ``self`` is invariant under.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
            sage: action = lambda g,x: M.term(g(x))
            sage: from sage.modules.with_basis.representation import Representation
            sage: R = Representation(G, M, action)
            sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
            sage: I = FiniteDimensionalInvariantModule(R)
            sage: I.semigroup()
            Symmetric group of order 3! as a permutation group

        """

        return self._semigroup

    def semigroup_representation(self):
        """
        Return the underlying representation of the invariant module.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: M = CombinatorialFreeModule(QQ, [1,2,3], prefix='M')
            sage: action = lambda g,x: M.term(g(x))
            sage: from sage.modules.with_basis.representation import Representation
            sage: R = Representation(G, M, action); R
            Representation of Symmetric group of order 3! as a permutation group indexed by {1, 2, 3} over Rational Field
            sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
            sage: I = FiniteDimensionalInvariantModule(R)
            sage: I.semigroup_representation()
            Representation of Symmetric group of order 3! as a permutation group indexed by {1, 2, 3} over Rational Field

        """

        return self._semigroup_representation

    #def _test_

    class Element(SubmoduleWithBasis.Element):

        def _mul_(self, other):
            """
            EXAMPLES::

                sage: M = CombinatorialFreeModule(QQ,[1,2,3],prefix='M');
                sage: G = CyclicPermutationGroup(3); G.rename('G')
                sage: g = G.an_element(); g
                (1,2,3)
                sage: from sage.modules.with_basis.representation import Representation
                sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
                sage: R = Representation(G,M,lambda g,x:M.monomial(g(x))); R.rename('R')
                sage: I = FiniteDimensionalInvariantModule(R)
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [M[1] + M[2] + M[3]]
                sage: v = B[0]
                sage: v*v
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: 'R' and 'R'
                sage: (1/2)*v
                1/2*B[0]
                sage: v*(1/2)
                1/2*B[0]

                sage: G = CyclicPermutationGroup(3); G.rename('G')
                sage: M = algebras.Exterior(QQ, 'x', 3)
                sage: from sage.modules.with_basis.representation import Representation
                sage: on_basis = lambda g,m: M.prod([M.monomial((g(j+1)-1,)) for j in m]) #cyclically permute generators
                sage: from sage.categories.algebras import Algebras
                sage: R = Representation(G, M, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional(), side = 'right')
                sage: I = FiniteDimensionalInvariantModule(R); I.rename('I')
                sage: B = I.basis()
                sage: v = B[0] + 2*B[1]; I.lift(v)
                1 + 2*x0 + 2*x1 + 2*x2
                sage: w = B[2]; I.lift(w)
                x0*x1 - x0*x2 + x1*x2
                sage: v*w
                B[2] + 6*B[3]
                sage: I.lift(v*w)
                x0*x1 + 6*x0*x1*x2 - x0*x2 + x1*x2 
                sage: w*v
                B[2] + 6*B[3]
                sage: (1/2)*v
                1/2*B[0] + B[1]
                sage: w*(1/2)
                1/2*B[2]
                sage: g = G((1,3,2))
                sage: v*g
                B[0] + 2*B[1]
                sage: w*g
                B[2]
                sage: g*v
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: 'G' and 'I'

                sage: R = Representation(G, M, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional())
                sage: I = FiniteDimensionalInvariantModule(R); I.rename('I')
                sage: B = I.basis()
                sage: v = B[0] + 2*B[1]; I.lift(v)
                1 + 2*x0 + 2*x1 + 2*x2
                sage: w = B[2]; I.lift(w)
                x0*x1 - x0*x2 + x1*x2
                sage: v*w
                B[2] + 6*B[3]
                sage: I.lift(v*w)
                x0*x1 + 6*x0*x1*x2 - x0*x2 + x1*x2 
                sage: w*v
                B[2] + 6*B[3]
                sage: (1/2)*v
                1/2*B[0] + B[1]
                sage: w*(1/2)
                1/2*B[2]
                sage: g = G((1,3,2))
                sage: v*v
                B[0] + 4*B[1]
                sage: g*w
                B[2]
                sage: v*g
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *: 'I' and 'G'

            """
            P = self.parent()
            try:
                return P.retract(P.lift(self) * P.lift(other))
            except:
                return P.retract(P.lift(self)*other)

        def _lmul_(self, right):
            """
            Give the product of ``self*right``

            EXAMPLES::

                sage: M = CombinatorialFreeModule(QQ,[1,2,3])
                sage: G = CyclicPermutationGroup(3)
                sage: g = G.an_element(); g
                (1,2,3)
                sage: from sage.modules.with_basis.representation import Representation
                sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
                sage: R = Representation(G,M,lambda g,x: M.monomial(g(x)), side = 'right')
                sage: I = FiniteDimensionalInvariantModule(R)
                sage: v = I.an_element(); v
                2*B[0]
                sage: v*g
                2*B[0]
                sage: [v*g for g in G.list()]
                [2*B[0], 2*B[0], 2*B[0]]

                sage: G = CyclicPermutationGroup(3)
                sage: M = algebras.Exterior(QQ, 'x', 3)
                sage: from sage.modules.with_basis.representation import Representation
                sage: on_basis = lambda g,m: M.prod([M.monomial((g(j+1)-1,)) for j in m]) #cyclically permute generators
                sage: from sage.categories.algebras import Algebras
                sage: R = Representation(G, M, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional(), side = 'right')
                sage: I = FiniteDimensionalInvariantModule(R)
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [1, x0 + x1 + x2, x0*x1 - x0*x2 + x1*x2, x0*x1*x2]
                sage: [[b*g for g in G] for b in B]
                [[B[0], B[0], B[0]],
                 [B[1], B[1], B[1]],
                 [B[2], B[2], B[2]],
                 [B[3], B[3], B[3]]]
                sage: 3*I.basis()[0]
                3*B[0]
                sage: 3*B[0] + B[1]*2
                3*B[0] + 2*B[1]

            """

            if right in self.parent()._semigroup and self.parent()._semigroup_representation.side() == 'right':
                return self

            elif right in self.parent()._semigroup_representation._module.base_ring():
                # This preserves the structure of the invariant as a 
                # ``.base_ring()``-module
                return self._mul_(right)

            return super()._lmul_(right)

        def _rmul_(self, left):
            """
            Give the product of ``left * self``

            EXAMPLES::

                sage: M = CombinatorialFreeModule(QQ,[1,2,3])
                sage: G = CyclicPermutationGroup(3)
                sage: g = G.an_element(); g
                (1,2,3)
                sage: from sage.modules.with_basis.representation import Representation
                sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
                sage: R = Representation(G,M,lambda g,x: M.monomial(g(x)))
                sage: I = FiniteDimensionalInvariantModule(R)
                sage: v = I.an_element(); v
                2*B[0]
                sage: g*v
                2*B[0]
                sage: [g*v for g in G.list()]
                [2*B[0], 2*B[0], 2*B[0]]

                sage: G = CyclicPermutationGroup(3)
                sage: M = algebras.Exterior(QQ, 'x', 3)
                sage: on_basis = lambda g,m: M.prod([M.monomial((g(j+1)-1,)) for j in m]) #cyclically permute generators
                sage: from sage.categories.algebras import Algebras
                sage: R = Representation(G, M, on_basis, category=Algebras(QQ).WithBasis().FiniteDimensional())
                sage: I = FiniteDimensionalInvariantModule(R)
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [1, x0 + x1 + x2, x0*x1 - x0*x2 + x1*x2, x0*x1*x2]
                sage: [[g*b for g in G] for b in B]
                [[B[0], B[0], B[0]], 
                 [B[1], B[1], B[1]], 
                 [B[2], B[2], B[2]], 
                 [B[3], B[3], B[3]]]
                sage: 3*I.basis()[0]
                3*B[0]
                sage: 3*B[0] + B[1]*2
                3*B[0] + 2*B[1]
            """
            if left in self.parent()._semigroup and self.parent()._semigroup_representation.side() == 'left':
                return self

            elif left in self.parent()._semigroup_representation._module.base_ring():
                return self._mul_(left)

            return super()._rmul_(left)

        def _acted_upon_(self, scalar, self_on_left = False):
            """
            EXAMPLES::
                
                sage: G = CyclicPermutationGroup(3)
                sage: M = CombinatorialFreeModule(QQ,[1,2,3])
                sage: from sage.modules.with_basis.representation import Representation
                sage: R = Representation(G, M, lambda g,x: M.monomial(g(x)))
                sage: from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
                sage: I = FiniteDimensionalInvariantModule(R)
                sage: B = I.basis()
                sage: [b._acted_upon_(G((1,3,2))) for b in B]
                [B[0]]

                sage: R = Representation(G, M, lambda g,x: M.monomial(g(x)), side = 'right')
                sage: I = FiniteDimensionalInvariantModule(R)
                sage: B = I.basis()
                sage: [b._acted_upon_(G((1,3,2)), self_on_left = True) for b in B]
                [B[0]]

                sage: R = G.regular_representation(QQ)
                sage: I = FiniteDimensionalInvariantModule(R)
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [() + (1,2,3) + (1,3,2)]
                sage: B[0]._acted_upon_(G((1,3,2)))
                B[0]
                sage: B[0]._acted_upon_(G((1,3,2)), self_on_left=True) == None
                True

                sage: R = G.regular_representation(QQ, side = 'right')
                sage: I = FiniteDimensionalInvariantModule(R)
                sage: B = I.basis()
                sage: [I.lift(b) for b in B]
                [() + (1,2,3) + (1,3,2)]
                sage: g = G((1,3,2))
                sage: B[0]._acted_upon_(g, self_on_left = True)
                B[0]
                sage: B[0]._acted_upon_(g, self_on_left = False) == None
                True

            """

            if scalar in self.parent()._semigroup and self_on_left == (self.parent()._semigroup_representation.side() == 'right'):

                return self

            return None