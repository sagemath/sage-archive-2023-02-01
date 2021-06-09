r"""
Invariant algebras
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

class FiniteDimensionalInvariantWithBasis(CombinatorialFreeModule):
    r"""
    When a group `G` acts on a module `M`, the invariant algebra is the collection 
    of elements `m` in `M` such that `g*m = m`.

    ..MATH::

        M^G = \{m \in M : g\cdot m = m}

    The current implementation works when `G` is a finitely-generated group, and
    when `M` is a finite-dimensional free module.
    """

    def __init__(self, M, G, action_on_basis = None, **kwds):
        """
        EXAMPLES::

        sage: V = VectorSpace(QQ,3)
        sage: G = CyclicPermutationGroup(3)
        sage: action = lambda x,g: x
        sage: VG = FiniteDimensionalInvariantWithBasis(V,G,action)
        Cyclic group of order 3 as a permutation group-invariant submodule of Vector space of dimension 3 over Rational Field
        sage: TestSuite(VG).run()

        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: G = CyclicPermutationGroup(3)
        sage: RG = FiniteDimensionalInvariantWithBasis(V,G,action)
        sage: TestSuite(RG).run()
        
        """
        # TODO: check assumtions on M
        self._ambient_module = M
        self._ambient_module_basis = M.basis()
        self._base_ring = M.base_ring()

        # TODO: add check for G to be finitely generated
        self._group = G

        # TODO: Determine if there are any checks to do here.
        self._action_on_basis = action_on_basis if action_on_basis else lambda x,g: x*g

        self._basis = M.annihilator_basis(G.gens(), lambda x,g: self._action_on_basis(x,g) - x)

        cat = kwds.pop("category", None)

        if cat is None:        
            cat = Modules(self._base_ring).FiniteDimensional().WithBasis()

        CombinatorialFreeModule.__init__(self, self._base_ring, self._basis, 
            category = cat)
        
    def _repr_(self):
        """
        EXAMPLES::
        
            sage: V = VectorSpace(QQ,3)
            sage: G = CyclicPermutationGroup(3)
            sage: action = lambda g, x: x
            sage: FiniteDimensionalInvariantWithBasis(V, G, action)
            (Cyclic group of order 3 as a permutation group)-invariant submodule of Vector space of dimension 3 over Rational Field

        """

        return f"({self._group})-invariant submodule of {self._ambient_module}"

    def base_ring(self):
        """
        EXAMPLES::

            sage: V = VectorSpace(QQ,3)
            sage: G = CyclicPermutationGroup(3)
            sage: VG = FiniteDimensionalInvariantWithBasis(V,G,lambda x,g: x)
            sage: VG.base_ring()
            Rational Field
        """

        return self._base_ring

