r"""
General Linear Groups

EXAMPLES:
    sage: GL(4,Q)
    General Linear Group of degree 4 over Rational Field
    sage: GL(1,Z)
    General Linear Group of degree 1 over Integer Ring
    sage: GL(100,RR)
    General Linear Group of degree 100 over Real Field with 53 bits of precision
    sage: GL(3,GF(49))
    General Linear Group of degree 3 over Finite Field in a of size 7^2

AUTHORS:
    -- David Joyner (2006-01)
    -- William Stein (2006-01)

"""

##TODO: Rework this and \code{special_linear} into MatrixGroup class for any
##field, wrapping all of GAP's matrix group commands in chapter 41
##Matrix Groups of the GAP reference manual.


#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

from sage.groups.group import Group
from sage.rings.all import Integer, is_Ring
from sage.rings.finite_field import is_FiniteField
from sage.matrix.matrix_space import MatrixSpace
from linear import LinearGroup_generic, LinearGroup_finite_field

def GL(n, R):
    """
    Return the general linear group of degree $n$ over the ring $R$.
    """
    if not is_Ring(R):
        raise TypeError, "R (=%) must be a ring"%R
    if is_FiniteField(R):
        return GeneralLinearGroup_finite_field(n, R)
    else:
        return GeneralLinearGroup_generic(n, R)

class GeneralLinearGroup_generic(LinearGroup_generic):
    def _gap_init_(self):
        return "GL(%s, %s)"%(self.degree(), self.base_ring().order())

    def _repr_(self):
        return "General Linear Group of degree %s over %s"%(self.degree(), self.base_ring())


class GeneralLinearGroup_finite_field(LinearGroup_finite_field, GeneralLinearGroup_generic):
    pass
