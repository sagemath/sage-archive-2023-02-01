r"""
Finitely generated commutative graded algebras with finite degree

AUTHORS:

- Michael Jung (2021): initial version

"""

#*****************************************************************************
#       Copyright (C) 2021 Michael Jung <m.jung at vu.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.all import Algebras
from sage.misc.cachefunc import cached_method
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
from sage.rings.ring import Algebra

class FiniteCommutativeGradedAlgebra(CombinatorialFreeModule, Algebra):
    r"""

    """

    Element = CombinatorialFreeModule.Element

    def __init__(self, base_ring, degrees, max_deg, names=None):
        r"""

        """
        from sage.arith.misc import gcd

        self._names = names
        self.__ngens = len(self._names)
        self._degrees = degrees
        self._max_deg = max_deg
        self._weighted_vectors = WeightedIntegerVectors(degrees)
        step = gcd(degrees)
        indices = [w for k in range(0, self._max_deg + 1, step)
                   for w in WeightedIntegerVectors(k, degrees)]
        sorting_key = self._weighted_vectors.grading
        cat = Algebras(base_ring).WithBasis().Graded().FiniteDimensional()
        CombinatorialFreeModule.__init__(self, base_ring, indices,
                                         sorting_key=sorting_key,
                                         category=cat)

    def _repr_(self):
        """

        """
        desc = f'Graded commutative algebra with generators {self._names} in '
        desc += f'degrees {self._degrees} with maximal finite '
        desc += f'degree {self._max_deg}'
        return desc

    def quotient(self, *args, **kwargs):
        r"""

        """
        raise NotImplementedError('no quotient implemented for {}'.format(self))

    @cached_method
    def product_on_basis(self, w1, w2):
        r"""
        Return the product of two basis vectors.

        """
        # TODO: Add graded-commutativity (-1)^(d_1 d_2)
        grading = self._weighted_vectors.grading
        deg_left = grading(w1)
        deg_right = grading(w2)
        deg_tot = deg_left + deg_right
        if deg_tot > self._max_deg:
            return self.zero()
        w_tot = self._weighted_vectors([sum(w) for w in zip(w1, w2)])
        return self.basis()[w_tot]

    def degree_on_basis(self, i):
        r"""
        Return the degree of the a homogeneous element with index `i`.

        """
        return self._weighted_vectors.grading(i)

    def _repr_term(self, w):
        r"""

        """
        # Trivial case:
        if sum(w) == 0:
            return '1'
        # Non-trivial case:
        res = ''
        for i in reversed(range(len(w))):
            if w[i] == 0:
                continue
            elif w[i] == 1:
                res += self._names[i] + '*'
            else:
                res += self._names[i] + '^{{{}}}'.format(w[i]) + '*'
        return res[:-1]

    def _latex_term(self, w):
        r"""

        """
        # Trivial case:
        if sum(w) == 0:
            return '1'
        # Non-trivial case:
        res = ''
        for i in reversed(range(len(w))):
            if w[i] == 0:
                continue
            elif w[i] == 1:
                res += self._names[i]
            else:
                res += self._names[i] + '^{{{}}}'.format(w[i])
        return res

    def algebra_generators(self):
        r"""
        Return the generators of ``self`` as a
        :class:`sage.sets.family.Family`.

        """
        from sage.sets.family import Family

        return Family(self.gens())

    @cached_method
    def one_basis(self):
        r"""
        Return the index of the one element of ``self``.

        """
        n = len(self._degrees)
        return self._weighted_vectors([0 for _ in range(n)])

    def gens(self):
        r"""
        Return the generators of ``self`` as a list.

        """
        n = len(self._degrees)
        zero = [0 for _ in range(n)]
        indices = []
        for k in range(n):
            ind = list(zero)
            ind[k] = 1
            indices.append(self._weighted_vectors(ind))
        return [self.monomial(ind) for ind in indices]

    @cached_method
    def gen(self, i):
        r"""
        Return the `i`-th generator of ``self``.

        """
        return self.gens()[i]
