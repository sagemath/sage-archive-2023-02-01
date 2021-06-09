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
from sage.categories.finite_dimensional_modules_with_basis import FiniteDimensionalModuleWithBasis


class InvariantModule(Parent):

    def __init__(self):
        pass


class FiniteDimensionalInvariantWithBasis(Parent, CombinatorialFreeModule):
    r"""
    When a group ``G`` acts on a module ``M``, the invariant algebra is the collection 
    of elements ``m`` in ``M`` such that ``g*m = m``.

    """

    def __init__(self, M, G, action_on_basis = None, R = ZZ):
        """
        EXAMPLES::
            
        """

        if not isinstance(G, ):

        if not isinstance(M, ):
        
        if not isinstance(action_on_basis, function):

        self._group = G
        self._action_on_basis = action_on_basis if action_on_basis else lambda x: x
        self._ambient_module = M

        cat = Modules(R).FiniteDimensional().WithBasis()

        if M.is_graded():
            cat = cat.Graded()


        
    def _repr_(self):
        """
        EXAMPLES::
        
            sage: V = VectorSpace(QQ,3)
            sage: G = CyclicPermutationGroup(3)
            sage: action = lambda g,x: x
            Cyclic group of order 3 as a permutation group-invariant submodule of Vector space of dimension 3 over Rational Field

        """

        return f"{self._group}-invariant submodule of {self._ambient_module}"

 