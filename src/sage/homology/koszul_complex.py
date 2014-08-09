"""
Koszul Complexes
"""

########################################################################
#       Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
###################################################.#####################

from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.choose_nk import rank
from sage.rings.arith import binomial
from sage.rings.all import ZZ
from sage.matrix.constructor import matrix
from sage.homology.chain_complex import ChainComplex_class

import itertools

class KoszulComplex(ChainComplex_class, UniqueRepresentation):
    r"""
    A Koszul complex.

    EXAMPLES::
    """
    @staticmethod
    def __classcall_private__(cls, R, elements=None):
        """
        Normalize input to ensure a unique representation.

        TESTS::
        """
        if elements is None:
            elements = R
            if not elements:
                R = ZZ # default to ZZ as the base ring if no elements are given
            else:
                R = elements[0].parent()
        return super(KoszulComplex, cls).__classcall__(cls, R, tuple(elements))

    def __init__(self, R, elements):
        """
        Initialize ``self``.

        EXAMPLES::
        """
        # Generate the differentials
        self._elements = elements
        n = len(elements)
        I = range(n)
        diff = {}
        zero = R.zero()
        for i in I:
            M = matrix(R, binomial(n,i+1), binomial(n,i), zero)
            j = 0
            for comb in itertools.combinations(I, i+1):
                for k,val in enumerate(comb):
                    r = rank(comb[:j] + comb[j+1:], n, False)
                    M[j,r] = (-1)**j * elements[val]
                j += 1
            M.set_immutable()
            diff[i+1] = M
        diff[0] = matrix(R, 0, 1, zero)
        diff[0].set_immutable()
        diff[n+1] = matrix(R, 1, 0, zero)
        diff[n+1].set_immutable()
        #print diff
        #from sage.homology.chain_complex import ChainComplex
        #print ChainComplex(diff, degree_of_differential=-1)
        #print diff
        ChainComplex_class.__init__(self, ZZ, ZZ(-1), R, diff)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::
        """
        return "Koszul complex defined by {} over {}".format(elements, self.base_ring())

