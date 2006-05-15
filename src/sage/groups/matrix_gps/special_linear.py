"""
Special Linear Groups

AUTHOR:
    -- William Stein: initial version
    -- David Joyner (2006-05) - added examples, _latex_, __str__, gens,
                                      as_matrix_group

EXAMPLES:
        sage: SL(2, ZZ)
        Modular Group SL(2,Z)
        sage: G = SL(2,GF(3)); G
        Special Linear Group of degree 2 over Finite Field of size 3
        sage: G.is_finite()
        True
        sage: G.conjugacy_class_representatives()
        [[1 0]
        [0 1], [0 2]
        [1 1], [0 1]
        [2 1], [2 0]
        [0 2], [0 2]
        [1 2], [0 1]
        [2 2], [0 2]
        [1 0]]
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.modular.all import SL2Z
from sage.rings.all import IntegerRing, is_FiniteField
from linear import LinearGroup_generic, LinearGroup_finite_field

def SL(n, R):
    if n == 2 and R == IntegerRing():
        return SL2Z()
    if is_FiniteField(R):
        return SpecialLinearGroup_finite_field(n, R)
    else:
        return SpecialLinearGroup_generic(n, R)

class SpecialLinearGroup_generic(LinearGroup_generic):
    def _gap_init_(self):
        """
        EXAMPLES:
            sage: G = SL(6,GF(5))
            sage: print G
            SL(6, GF(5))

        """
        return "SL(%s, GF(%s))"%(self.degree(), self.base_ring().order())

    def _latex_(self):
        """
        EXAMPLES:
            sage: G = SL(6,GF(5))
            sage: G._latex_()
            'SL$(6, GF(5))$'

        """
        return "SL$(%s, GF(%s))$"%(self.degree(), self.base_ring().order())

    def __str__(self):
        """
        EXAMPLES:
            sage: G = SL(6,GF(5))
            sage: print G
            SL(6, GF(5))

        """
        return "SL(%s, GF(%s))"%(self.degree(), self.base_ring().order())

    def __repr__(self):
        return "Special Linear Group of degree %s over %s"%(self.degree(), self.base_ring())

    def gens(self):
        """
        EXAMPLES:
            sage: G = SL(6,GF(5))
            sage: G.gens()
            [[2 0 0 0 0 0]
            [0 3 0 0 0 0]
            [0 0 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 0]
            [0 0 0 0 0 1],
            [4 0 0 0 0 1]
            [4 0 0 0 0 0]
            [0 4 0 0 0 0]
            [0 0 4 0 0 0]
            [0 0 0 4 0 0]
            [0 0 0 0 4 0]]
        """
        from sage.interfaces.all import gap
        F = self.base_ring()
        G = self._gap_init_()
        n = eval(gap.eval("Length(GeneratorsOfGroup(%s))"%G))
        gens = [gap("GeneratorsOfGroup(%s)[%s]"%(G,i))._matrix_(F) for i in range(1,n+1)]
        return gens

    def as_matrix_group(self):
        """
        EXAMPLES:
            sage: G = SL(6,GF(5))
            sage: G.as_matrix_group()
            Matrix group over Finite Field of size 5 with 2 generators:
            [[[2, 0, 0, 0, 0, 0], [0, 3, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]], [[4, 0, 0, 0, 0, 1], [4, 0, 0, 0, 0, 0], [0, 4, 0, 0, 0, 0], [0, 0, 4, 0, 0, 0], [0, 0, 0, 4, 0, 0], [0, 0, 0, 0, 4, 0]]]
        """
        from sage.groups.matrix_gps.matrix_group import MatrixGroup
        gns = self.gens()
        G = MatrixGroup(gns)
        return G

class SpecialLinearGroup_finite_field(SpecialLinearGroup_generic, LinearGroup_finite_field):
    pass

