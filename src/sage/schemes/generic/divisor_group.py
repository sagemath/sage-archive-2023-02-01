"""
AUTHORS:
"""

#*******************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from sage.groups.group import AbelianGroup
from divisor import Divisor_curve
from sage.structure.parent_gens import ParentWithGens
from sage.rings.integer_ring import ZZ

class DivisorGroup(AbelianGroup):
    def __init__(self, scheme, value_ring=None):
        ParentWithGens.__init__(self, ZZ)
        self.__scheme = scheme
        if value_ring is None:
            value_ring = scheme.base_ring()
        self.__value_ring = value_ring

    def _repr_(self):
        return "Group of Divisors on %s"%self.__scheme

    def _coerce_impl(self, x):
        if x.parent() is self:
            return x
        raise TypeError

    def __cmp__(self, right):
        if not isinstance(right, DivisorGroup):
            return -1
        return cmp(self.__scheme, right.__scheme)

    def scheme(self):
        return self.__scheme

    def _an_element_impl(self):
        return self.__scheme.divisor([], check=False, reduce=False)

class DivisorGroup_curve(DivisorGroup):
    def __call__(self, v):
        return Divisor_curve(v)

