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

from sage.misc.cachefunc import cached_method
from sage.modules.with_basis.subquotient import SubmoduleWithBasis
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finitely_generated_semigroups import FinitelyGeneratedSemigroups
from sage.categories.groups import Groups

class FiniteDimensionalInvariantModule(SubmoduleWithBasis):
    r"""
    Construct the `S`-invariant submodule of `M`. When a semigroup `S` acts on a module
    `M`, the invariant module is the collection of elements `m` in `M` such that
    `s \cdot m = m` for all `s \in S.

    ..MATH::

        M^S = \{m \in M : s\cdot m = m,\, \forall s \in S \}

    NOTE: The current implementation works when `S` is a finitely-generated semigroup,
    and when `M` is a finite-dimensional free module with a distinguished basis.

    TODO::
        Extend when `M` does not have a basis and `S` is a permutation group using:
        - https://arxiv.org/abs/0812.3082
        - https://www.dmtcs.org/pdfpapers/dmAA0123.pdf

    """

    def __init__(self, R, *args, **kwargs):
        """
        INPUTS::

        - ``R`` -- an instance of a ``Representation`` of a semigroup `S`
                   acting on the module `M`.

        OUTPUTS::

        - ``MS`` -- the invariant algebra of the semigroup action of `S` on `M`, or
                    equivalently, the isotypic component of the representation of
                    `S` carried by `M` corresponding to the trivial character.

        EXAMPLES::
        
            

        TESTS::

        """
     
        self._semigroup = R.semigroup()
        
        if self._semigroup not in FinitelyGeneratedSemigroups:
            raise ValueError(f'{self._semigroup} is not finitely generated')

        self._ambient_module = R._module

        if self._ambient_module not in FiniteDimensionalModulesWithBasis:
            raise ValueError(f'{self._ambient_module} is not finite dimensional')

        self._on_ambient_basis = R._on_basis

        self._side = kwargs.pop('side', R.side())

        unitriangular = kwargs.pop('unitriangular', False)
        category = kwargs.pop('category', R.category().Subobjects())

        if self._side == 'left':            
            self._invariant_map = lambda g,x: self._on_ambient_basis(g,x.support()[0]) - x
        else:
            self._invariant_map = lambda x,g: self._on_ambient_basis(x.support()[0],g) - x

        
        basis = self._ambient_module.annihilator_basis(
                    self._semigroup.gens(),
                    action = self._invariant_map,
                    side = self._side)


        super().__init__(Family(basis), support_order = self._ambient_module._compute_support_order(basis),
                                                ambient = self._ambient_module, # bug - should lift to the representation
                                                unitriangular = unitriangular,#is this right?
                                                category = category, 
                                                *args, **kwargs)

    def _repr_(self):
        """
        EXAMPLES::

            sage: V = VectorSpace(QQ,3)
            sage: G = CyclicPermutationGroup(3)
            sage: FiniteDimensionalInvariantAlgebra(V, G)
            (Cyclic group of order 3 as a permutation group)-invariant subalgebra of Vector space of dimension 3 over Rational Field

        """

        return f"({self._semigroup})-invariant submodule of {self._ambient_module}"


    def semigroup(self):
        """
        Return the semigroup `S` whose action ``self`` is invariant under.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: A.<x,y,z> = PolynomialRing(QQ)
            sage: import operator
            sage: AG = FiniteDimensionalInvariantAlgebra(G,A,operator.mul)
            sage: AG.group()
            Symmetric group of order 3! as a permutation group

        """

        return self._semigroup

    def action_on_basis(self,*args):

        if len(args) == 0:
            return self._action_on_basis

        return self._action_on_basis(args)

    class Element(SubmoduleWithBasis.Element):
        
        def _mul_(self, other):
            P = self.parent()
            return P.retract(P.lift(self) * P.lift(other))
        
        # lmul -- self*right
        def _lmul_(self, right):
            """
            EXAMPLES::

                sage:        #### create invariant module
                sage:        #### multiply group element
            """

            if right in self.parent()._semigroup and self.parent()._side == 'right':
                return self

            return super()._lmul_(right)

        # rmul -- left * self
        def _rmul_(self, left):
            """
            EXAMPLES::

                sage: ### create invariant module
                sage: ### multiply group element
            """
            if left in self.parent()._semigroup and self.parent()._side == 'left':
                return self
            
            return super()._rmul_(left)

        def _acted_upon_(self, scalar, self_on_left = False):
            """
            EXAMPLES::

                sage: 
            """

            if scalar in self.parent()._semigroup and self_on_left == (self.parent()._side == 'right'):

                return self

            return None
