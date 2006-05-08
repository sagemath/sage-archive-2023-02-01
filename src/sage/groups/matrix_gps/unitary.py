r"""
 Unitary Groups $GU(n,q)$ and $SU(n,q)$

These are $n \times n$ unitary matrices with entries
in $GF(q^2)$.

AUTHOR:
    -- David Joyner: initial version (2006-3), modified from
                     special_linear (by W. Stein)

EXAMPLES:
    sage: G = SU(3,GF(5))
    sage: G
    Special Unitary Group of degree 3 over Finite Field of size 5
    sage: G._gap_init_()
    'SU(3, 5)'
    sage: G.random()
    [      3 2*a + 1   a + 1]
    [2*a + 2       4       a]
    [4*a + 2   a + 4       1]
    sage: G.base_ring()
    Finite Field of size 5
    sage: G.field_of_definition()
    Finite Field in a of size 5^2

"""

#*********************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*********************************************************************************

from sage.rings.all import IntegerRing, is_FiniteField
from linear import LinearGroup_generic, LinearGroup_finite_field

def GU(n, R):
    if is_FiniteField(R):
        return GeneralUnitaryGroup_finite_field(n, R)
    else:
        return GeneralUnitaryGroup_generic(n, R)

class GeneralUnitaryGroup_generic(LinearGroup_generic):
    def _gap_init_(self):
        return "GU(%s, %s)"%(self.degree(), self.base_ring().order())

    def _repr_(self):
        return "General Unitary Group of degree %s over %s"%(self.degree(), self.base_ring())

class GeneralUnitaryGroup_finite_field(GeneralUnitaryGroup_generic, LinearGroup_finite_field):
    pass

def SU(n, R):
    if is_FiniteField(R):
        return SpecialUnitaryGroup_finite_field(n, R)
    else:
        return SpecialUnitaryGroup_generic(n, R)

class SpecialUnitaryGroup_generic(LinearGroup_generic):
    def _gap_init_(self):
        return "SU(%s, %s)"%(self.degree(), self.base_ring().order())

    def _repr_(self):
        return "Special Unitary Group of degree %s over %s"%(self.degree(), self.base_ring())

class SpecialUnitaryGroup_finite_field(SpecialUnitaryGroup_generic, LinearGroup_finite_field):
    pass

