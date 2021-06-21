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

from sage.modules.with_basis.invariant import FiniteDimensionalInvariantModule
from sage.categories.algebras import Algebras

class FiniteDimensionalInvariantAlgebra(FiniteDimensionalInvariantModule):

    def __init__(self, R, *args, **kwargs):
        """
        INPUTS::

        - ``R`` -- an instance of a ``Representation`` of a semigroup `S`
                   acting on the algebra `A`.

        OUTPUTS::

        - ``AS`` -- the invariant algebra of the semigroup action of `S` on `A`, or
                    equivalently, the isotypic component of the representation of
                    `S` carried by `A` corresponding to the trivial character.

        EXAMPLES::

            sage: S = Set({1,2,3})
            sage: L = S.subsets_lattice()
            sage: M = L.moebius_algebra(QQ)
            sage: A = M.a_realization()
            sage: G = SymmetricGroup(3)
            sage: on_basis = lambda g,x: A.monomial(Set(g(s) for s in S))
            sage: from sage.modules.with_basis.representation import Representation
            sage: R = Representation(G, A, on_basis)
            sage: from sage.algebras.finite_dimensional_algebras.finite_dimensional_invariant_algebras import FiniteDimensionalInvariantAlgebra
            sage: I = FiniteDimensionalInvariantAlgebra(R)
            sage: I.basis()
            (E[{}], E[{1}] + E[{2}] + E[{3}], E[{1,2}] + E[{1,3}] + E[{2,3}], E[{1,2,3}])
            sage: e = I.an_element(); e

            sage: G = SymmetricGroup(4)
            sage: table = [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])]
            sage: A = FiniteDimensionalAlgebra(GF(3), table, names = 'e')
            sage: R = Representation(G, A, on_basis)
            sage: I = FiniteDimensionalInvariantAlgebra(R)

            sage: G = CyclicPermutationGroup(3)
            sage: M = algebras.Exterior(QQ, 'x', 3)
            sage: from sage.modules.with_basis.representation import Representation
            sage: on_basis = lambda g,m: M.prod([M.monomial(tuple([g(j+1)-1])) for j in m]) #cyclically permute generators
            sage: from sage.categories.algebras import Algebras
            sage: R = Representation(G, M, on_basis, category = Algebras(QQ).WithBasis().FiniteDimensional())
            sage: I = FiniteDimensionalInvariantAlgebra(R); I
            (Cyclic group of order 3 as a permutation group)-invariant subalgebra of
             The exterior algebra of rank 3 over Rational Field
            sage: I.basis()
            (1, x0 + x1 + x2, x0*x1*x2)

        TESTS::

            sage: G = CyclicPermutationGroup(3)
            sage: M = algebras.Exterior(QQ, 'x', 3)
            sage: from sage.modules.with_basis.representation import Representation
            sage: on_basis = lambda g,m: M.prod([M.monomial(tuple([g(j+1)-1])) for j in m]) #cyclically permute generators
            sage: from sage.categories.algebras import Algebras
            sage: R = Representation(G, M, on_basis, category = Algebras(QQ))
            sage: from sage.algebras.finite_dimensional_algebras.finite_dimensional_invariant_algebras \
            ....: import FiniteDimensionalInvariantAlgebra
            sage: I = FiniteDimensionalInvariantAlgebra(R)
            Traceback (most recent call last):
            ...
            ValueError: 'category' keyword argument must be FiniteDimensional

            sage: R = Representation(G, M, on_basis, category = Algebras(QQ).FiniteDimensional())
            sage: I = FiniteDimensionalInvariantAlgebra(R)
            Traceback (most recent call last):
            ...
            ValueError: 'category' keyword argument must be WithBasis

        """

        # Check that the representation indeed comes from an algebra
        if R._module not in Algebras():
            raise ValueError(f'{R._module} is not an algebra')

        if R._module not in Algebras().FiniteDimensional().WithBasis():
            raise NotImplementedError(f'{R._module} must be finite-dimensional with a basis')

        if 'category' not in kwargs:
            category = Algebras().FiniteDimensional().WithBasis()
            kwargs['category'] = category
        else:
            if 'FiniteDimensional' not in kwargs['category'].axioms():
                raise ValueError("'category' keyword argument must be FiniteDimensional")
            if 'WithBasis' not in kwargs['category'].axioms():
                raise ValueError("'category' keyword argument must be WithBasis")

        super().__init__(R, *args, **kwargs)

    def _repr_(self):
        """
        EXAMPLES::

        """

        return f"({self._semigroup})-invariant subalgebra of {self._ambient_module}"
