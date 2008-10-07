"""
Sum species

"""
#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from species import GenericCombinatorialSpecies
from structure import SpeciesStructureWrapper
import __builtin__
from sage.misc.cachefunc import cached_function


class SumSpeciesStructure(SpeciesStructureWrapper):
    pass

@cached_function
def SumSpecies(*args, **kwds):
    """
    Returns the sum of two species.

    EXAMPLES:
        sage: S = species.PermutationSpecies()
        sage: A = S+S
        sage: A.generating_series().coefficients(5)
        [2, 2, 2, 2, 2]

    TESTS:
        sage: A = species.SingletonSpecies() + species.SingletonSpecies()
        sage: B = species.SingletonSpecies() + species.SingletonSpecies()
        sage: C = species.SingletonSpecies() + species.SingletonSpecies(min=2)
        sage: A is B
        True
        sage: (A is C) or (A == C)
        False

    """
    return SumSpecies_class(*args, **kwds)

class SumSpecies_class(GenericCombinatorialSpecies):
    def __init__(self, F, G, min=None, max=None, weight=None):
        """
        EXAMPLES:
            sage: P = species.PermutationSpecies()
            sage: F = P + P
            sage: F._check()
            True
            sage: F == loads(dumps(F))
            True
        """
        self._F = F
        self._G = G

        self._state_info = [F, G]

        GenericCombinatorialSpecies.__init__(self, min=None, max=None, weight=None)

    _default_structure_class = SumSpeciesStructure

    _cached_constructor = staticmethod(SumSpecies)

    def _name(self):
        """
        Note that we use a function to return the name of this species
        because we can't do it in the __init__ method due to it requiring
        that self._F and self._G already be unpickled.

        EXAMPLES:
            sage: P = species.PermutationSpecies()
            sage: F = P + P
            sage: F._name()
            'Sum of (Permutation species) and (Permutation species)'
        """
        return "Sum of (%s) and (%s)"%(self._F, self._G)

    def _structures(self, structure_class, labels):
        """
        EXAMPLES:
            sage: P = species.PermutationSpecies()
            sage: F = P + P
            sage: F.structures([1,2]).list()
            [[1, 2], [2, 1], [1, 2], [2, 1]]
        """
        for res in self._F.structures(labels):
            yield structure_class(self, res, tag="left")

        for res in self._G.structures(labels):
            yield structure_class(self, res, tag="right")

    def _isotypes(self, structure_class, labels):
        """
        EXAMPLES:
            sage: P = species.PermutationSpecies()
            sage: F = P + P
            sage: F.isotypes([1,2]).list()
            [[2, 1], [1, 2], [2, 1], [1, 2]]

        """
        for res in self._F.isotypes(labels):
            yield structure_class(self, res, tag="left")

        for res in self._G.isotypes(labels):
            yield structure_class(self, res, tag="right")

    def _gs(self, series_ring, base_ring):
        """
        Returns the cycle index series of this species.

        EXAMPLES:
            sage: P = species.PermutationSpecies()
            sage: F = P + P
            sage: F.generating_series().coefficients(5)
            [2, 2, 2, 2, 2]
        """
        return self._F.generating_series(base_ring) + self._G.generating_series(base_ring)


    def _itgs(self, series_ring, base_ring):
        """
        Returns the isomorphism type generating series of this species.

        EXAMPLES:
            sage: P = species.PermutationSpecies()
            sage: F = P + P
            sage: F.isotype_generating_series().coefficients(5)
            [2, 2, 4, 6, 10]
        """
        return (self._F.isotype_generating_series(base_ring) +
                self._G.isotype_generating_series(base_ring))

    def _cis(self, series_ring, base_ring):
        """
        Returns the generating series of this species.

        EXAMPLES:
            sage: P = species.PermutationSpecies()
            sage: F = P + P
            sage: F.cycle_index_series().coefficients(5)
            [2*p[],
             2*p[1],
             2*p[1, 1] + 2*p[2],
             2*p[1, 1, 1] + 2*p[2, 1] + 2*p[3],
             2*p[1, 1, 1, 1] + 2*p[2, 1, 1] + 2*p[2, 2] + 2*p[3, 1] + 2*p[4]]

        """
        return self._F.cycle_index_series(base_ring) + self._G.cycle_index_series(base_ring)

    def weight_ring(self):
        """
        Returns the weight ring for this species.  This is determined by asking
        Sage's coercion model what the result is when you add elements of the
        weight rings for each of the operands.

        EXAMPLES:
            sage: S = species.SetSpecies()
            sage: C = S+S
            sage: C.weight_ring()
            Rational Field

            sage: S = species.SetSpecies(weight=QQ['t'].gen())
            sage: C = S + S
            sage: C.weight_ring()
            Univariate Polynomial Ring in t over Rational Field
        """
        return self._common_parent([self._F.weight_ring(), self._G.weight_ring()])

    def _equation(self, var_mapping):
        """
        Returns the right hand side of an algebraic equation satisfied by
        this species.  This is a utility function called by the
        algebraic_equation_system method.

        EXAMPLES:
            sage: X = species.SingletonSpecies()
            sage: S = X + X
            sage: S.algebraic_equation_system()
            [node1 - 2*z]

        """
        return sum(var_mapping[operand] for operand in self._state_info)

