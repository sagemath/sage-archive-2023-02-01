"""
Symplectic Linear Groups

AUTHOR:
    -- David Joyner: initial version (2006-3), modified from
                     special_linear (by W. Stein)

EXAMPLES:
    sage: G = Sp(4,GF(7))
    sage: G._gap_init_()
    'Sp(4, 7)'
    sage: G
    Symplectic Group of rank 2 over Finite Field of size 7
    sage: G.random()
    [5 5 5 1]
    [0 2 6 3]
    [5 0 1 0]
    [4 6 3 4]
    sage: G.order()
    276595200

"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.modular.all import SL2Z
from sage.rings.all import IntegerRing, is_FiniteField
from linear import LinearGroup_generic, LinearGroup_finite_field

def Sp(n, R):
    if n%2!=0:
        raise ValueError, "\n n must be even.\n"
    if n == 2 and R == IntegerRing():
        return SL2Z()
    if is_FiniteField(R):
        return SymplecticGroup_finite_field(n, R)
    else:
        return SymplecticGroup_generic(n, R)

class SymplecticGroup_generic(LinearGroup_generic):
    def _gap_init_(self):
        return "Sp(%s, %s)"%(self.degree(), self.base_ring().order())

    def _repr_(self):
        return "Symplectic Group of rank %s over %s"%(self.degree()/2, self.base_ring())

class SymplecticGroup_finite_field(SymplecticGroup_generic, LinearGroup_finite_field):
    pass

