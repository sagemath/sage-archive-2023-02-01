"""
Symplectic Linear Groups

AUTHORS:

- David Joyner (2006-03): initial version, modified from
  special_linear (by W. Stein)

EXAMPLES::

    sage: G = Sp(4,GF(7))
    sage: G._gap_init_()
    'Sp(4, 7)'
    sage: G
    Symplectic Group of rank 2 over Finite Field of size 7
    sage: G.random_element()
    [1 6 5 5]
    [2 1 4 5]
    [1 2 4 5]
    [4 0 2 2]
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

    EXAMPLES::

        sage: Sp(4,5)
        Symplectic Group of rank 2 over Finite Field of size 5
        sage: Sp(3,GF(7))
        Traceback (most recent call last):
        ...
        ValueError: the degree n (=3) must be even
    """
    if n%2!=0:
        raise ValueError, "the degree n (=%s) must be even"%n
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
        r"""
        Return LaTeX representation of this group.

        EXAMPLES::

            sage: latex(Sp(4,5))
            \text{Sp}_{4}(\Bold{F}_{5})
        """
        return "\\text{Sp}_{%s}(%s)"%(self.degree(), self.field_of_definition()._latex_())

    def _repr_(self):
        """
        Return print representation of this group.

        EXAMPLES::

            sage: Sp(2,4)
            Symplectic Group of rank 1 over Finite Field in a of size 2^2
        """
        return "Symplectic Group of rank %s over %s"%(self.degree()/2, self.base_ring())

class SymplecticGroup_finite_field(SymplecticGroup_generic, MatrixGroup_gap_finite_field):
    def _gap_init_(self):
        """
        Return GAP string that evaluates to this group.

        EXAMPLES::

            sage: Sp(2,4)._gap_init_()
            'Sp(2, 4)'
        """
        return "Sp(%s, %s)"%(self.degree(), self.base_ring().order())



