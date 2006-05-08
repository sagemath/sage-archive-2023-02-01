"""
Special Linear Groups

AUTHOR:
    -- William Stein: initial version

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
        return "SL(%s, %s)"%(self.degree(), self.base_ring().order())

    def _repr_(self):
        return "Special Linear Group of degree %s over %s"%(self.degree(), self.base_ring())

class SpecialLinearGroup_finite_field(SpecialLinearGroup_generic, LinearGroup_finite_field):
    pass

