"""
Abelian groups

AUTHOR:
    - William Stein: initial version
    - David Joyner (2006-01-01): added examples,
    - William Stein (2006-01-03): added examples and _latex_ method
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.misc.misc
import group

class AbelianGroup(group.AbelianGroup):
    r"""
    A finitely generated abstract abelian group given
    by integer invariants.

    Use $0$ for each infinite cyclic factor.

    EXAMPLES:
    We create the group $\Z/2\Z \times \Z/3\Z$.

        sage: AbelianGroup([2,3])
        Abelian Group with invariants [2, 3]
        sage: G = AbelianGroup([2,3]); G
        Abelian Group with invariants [2, 3]
        sage: G.invariants()
        [2, 3]
        sage: G.is_abelian()
        True
        sage: G.exponent()
        6

    We create $\Z/3\Z \times \Z/3\Z$.
        sage: G = AbelianGroup([3,3]); G
        Abelian Group with invariants [3, 3]
        sage: G.order()
        9
        sage: G.exponent()
        3
    """
    def __init__(self, invariants):
        invariants.sort()
        invariants = [x for x in invariants if x != 1]
        self.__invariants = invariants

    def _latex_(self):
        """
        Return latex representation of self.

        EXAMPLES:
            sage: G = AbelianGroup([2,3,7]); G
            Abelian Group with invariants [2, 3, 7]
            sage: latex(G)
            '\\Z/2\\Z \\times \\Z/3\\Z \\times \\Z/7\\Z'
        """
        return ' \\times '.join(['\\Z/%s\Z'%a for a in self.__invariants])

    def _repr_(self):
        """
        Return representation of self.

        EXAMPLES:
            sage: G = AbelianGroup([2,2])
            sage: G
            Abelian Group with invariants [2, 2]
            sage: G.rename('Klein 4 group')
            sage: G
            Klein 4 group
        """
        return "Abelian Group with invariants %s"%self.__invariants

    def exponent(self):
        """
        Return the exponent of this group.

        EXAMPLES:
            sage: G = AbelianGroup([2,3,7]); G
            Abelian Group with invariants [2, 3, 7]
            sage: G.exponent()
            42
        """
        import sage.rings.all
        return sage.rings.all.LCM(self.__invariants)

    def invariants(self):
        """
        Return the invariants of this group.

        EXAMPLES:
            sage: G = AbelianGroup([2,3,7]); G
            Abelian Group with invariants [2, 3, 7]
            sage: G.invariants()
            [2, 3, 7]
        """
        return self.__invariants

    def is_commutative(self):
        """
        Return True, since abelian groups are commutative.

        EXAMPLES:
            sage: G = AbelianGroup([2,2,3])
            sage: G.is_commutative()
            True
        """
        return True

    def is_abelian(self):
        """
        Return True, since abelian groups are abelian.

        EXAMPLES:
            sage: G = AbelianGroup([2,2,3])
            sage: G.is_commutative()
            True
        """
        return True

    def order(self):
        """
        Return the order of this group.

        EXAMPLES:
            sage: G = AbelianGroup([2,3])
            sage: G.order()
            6
            sage: G = AbelianGroup([2,3,0])
            sage: G.order()
            Infinity
        """
        import sage.rings.all
        try:
            return self.__len
        except AttributeError:
            self.__len = sage.misc.misc.mul(self.__invariants)
            if self.__len == 0:
                self.__len = sage.rings.all.infinity
        return self.__len

    def permutation_group(self):
        r"""
        Return the permutation group isomorphic to this abelian group.

        If the invariants are $q_1, \ldots, q_n$ then the generators
        of the permutation will be of order $q_1, \ldots, q_n$,
        respectively.

        EXAMPLES:
            sage: G = AbelianGroup([2,3]); G
            Abelian Group with invariants [2, 3]
            sage: G.permutation_group()
            Permutation Group with generators [(1,4)(2,5)(3,6), (1,2,3)(4,5,6)]
        """
        from permgroup import PermutationGroup
        from sage.interfaces.all import gap
        invs = self.__invariants
        s = 'Image(IsomorphismPermGroup(AbelianGroup(%s)))'%invs
        return PermutationGroup(gap(s), from_group=True)

