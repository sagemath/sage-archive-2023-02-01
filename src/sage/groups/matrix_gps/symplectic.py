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

from sage.rings.all import IntegerRing, is_FiniteField, Integer, FiniteField
from matrix_group import MatrixGroup_gap, MatrixGroup_gap_finite_field

def Sp(n, R, var='a'):
    """
    Return the symplectic group of degree n over R.

    EXAMPLES:
        sage: ?
    """
    if n%2!=0:
        raise ValueError, "n must be even"
    if isinstance(R, (int, long, Integer)):
        R = FiniteField(R, var)
    if is_FiniteField(R):
        return SymplecticGroup_finite_field(n, R)
    else:
        return SymplecticGroup_generic(n, R)

class SymplecticGroup_generic(MatrixGroup_gap):
    def _gap_init_(self):
        raise TypeError, 'no analogue of this symplectic group in GAP'

    def _latex_(self):
        """
        Return LaTeX representation of this group.

        EXAMPLES:
            sage: ?
        """
        return "\\text{Sp}_{%s}(%s)"%(self.degree(), self.field_of_definition()._latex_())

    def _repr_(self):
        """
        Return print representation of this group.
        """
        return "Symplectic Group of rank %s over %s"%(self.degree()/2, self.base_ring())

class SymplecticGroup_finite_field(SymplecticGroup_generic, MatrixGroup_gap_finite_field):
    def _gap_init_(self):
        """
        Return GAP string that evaluates to this group.

        EXAMPLES:
        sage: ?
        """
        return "Sp(%s, %s)"%(self.degree(), self.base_ring().order())



