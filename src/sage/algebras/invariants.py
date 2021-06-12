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

from sage.misc.cachefunc import cached_method
from sage.modules.with_basis.subquotient import SubmoduleWithBasis
from sage.structure.unique_representation import UniqueRepresentation

class FiniteDimensionalInvariantAlgebra(UniqueRepresentation, SubmoduleWithBasis):
    r"""
    Construct the `G`-invariant subalgebra of `A`. When a group `G` acts on an algebra 
    `A`, the invariant algebra is the collection of elements `a` in `A` such that 
    `g*a = a`.

    ..MATH::

        A^G = \{a \in A : g\cdot a = a}

    NOTE: The current implementation works when `G` is a finitely-generated group, and
    when `A` is a finite-dimensional free module.

    TODO::
        Extend when `A` does not have a basis and `G` is a permutation group using:
        - https://arxiv.org/abs/0812.3082
        - https://www.dmtcs.org/pdfpapers/dmAA0123.pdf

    """

    def __init__(self, A, G, action_on_basis = lambda x, g: x, **kwargs):
        """

        TESTS::

            sage: from sage.algebras.invariants import FiniteDimensionalInvariantAlgebra
            sage: A.<x,y,z> = PolynomialRing(QQ)
            sage: G = SymmetricGroup(3); G.rename('S3')
            sage: import operator
            sage: AG = FiniteDimensionalInvariantAlgebra(A, G, operator.mul); AG
            (S3)-invariant subalgebra of Polynomial Ring over Rational Field

        """

        self._group = G
        self._action_on_basis = action_on_basis

        side = kwargs.pop('side','right')
        support_order = kwargs.pop('support_order', None)
        unitriangular = kwargs.pop('unitriangular', False)
        category = kwargs.pop('category', None)

        basis = A.annihilator_basis(G.gens(), action = lambda x,g: action_on_basis(x,g) - x, side=side)

        super(SubmoduleWithBasis,self).__init__(basis.keys(),
                                                support_order = support_order,
                                                ambient = A,
                                                unitriangular = unitriangular,
                                                category = category,
                                                **kwargs)

    def _repr_(self):
        """
        EXAMPLES::

            sage: V = VectorSpace(QQ,3)
            sage: G = CyclicPermutationGroup(3)
            sage: FiniteDimensionalInvariantAlgebra(V, G)
            (Cyclic group of order 3 as a permutation group)-invariant subalgebra of Vector space of dimension 3 over Rational Field

        """

        return f"({self._group})-invariant subalgebra of {self._ambient_module}"


    def group(self):
        """
        Return the group `G` whose action ``self`` is invariant under.

        EXAMPLES::

            sage: G = SymmetricGroup(3)
            sage: A.<x,y,z> = PolynomialRing(QQ)
            sage: import operator
            sage: AG = FiniteDimensionalInvariantAlgebra(G,A,operator.mul)
            sage: AG.group()
            Symmetric group of order 3! as a permutation group

        """

        return self._group