#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from group import Group
from sage.rings.all import IntegerRing, is_Ring, Integer


class LinearGroup_generic(Group):
    def __init__(self, n, R):
        self.__n = Integer(n)
        self.__R = R
        if not is_Ring(R):
            raise TypeError, "R (=%s) must be a ring"%R

    def degree(self):
        return self.__n

    def base_ring(self):
        return self.__R

    def is_finite(self):
        """
        EXAMPLES:
            sage: G = GL(2,GF(3))
            sage: G.is_finite()
            True
        """
        return self.__R.is_finite()

class LinearGroup_finite_field(LinearGroup_generic):
    def conjugacy_class_representatives(self):
        """
        Return a complete set of representatives of the conjugacy
        classes of the group.

        EXAMPLES:
            sage: G = GL(2,GF(3))
            sage: C = G.conjugacy_class_representatives()
            sage: len(C)
            8
            sage: C[0]
            [1 0]
            [0 1]
            sage: C
            [[1 0]
            [0 1], [0 1]
            [2 1], [2 0]
            [0 2], [0 1]
            [2 2], [0 1]
            [2 0], [0 1]
            [1 2], [0 1]
            [1 1], [2 0]
            [0 1]]
            sage: G = GL(2,GF(4))
            sage: C = G.conjugacy_class_representatives()
            sage: [list(g) for g in C]      # prints more nicely
            [[(1, 0), (0, 1)], [(0, 1), (1, 0)], [(a, 0), (0, a)], [(0, 1), (a + 1, 0)],
            [(a + 1, 0), (0, a + 1)], [(0, 1), (a, 0)], [(0, 1), (1, a)],
            [(0, 1), (1, a + 1)], [(0, 1), (a, 1)], [(0, 1), (a, a)],
            [(0, 1), (a + 1, 1)], [(0, 1), (a + 1, a + 1)], [(1, 0), (0, a)],
            [(1, 0), (0, a + 1)], [(a, 0), (0, a + 1)]]
        """
        G = self._gap_()
        C = G.ConjugacyClasses()
        gap = G.parent()
        reps = gap.List(C, 'x -> Representative(x)')
        K = self.base_ring()
        self.__reps = [x._matrix_(K) for x in reps]
        return self.__reps


