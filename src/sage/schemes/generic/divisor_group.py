"""
AUTHORS:
"""

#*******************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from sage.groups.group import AbelianGroup

class DivisorGroup(AbelianGroup):
    def __init__(self, scheme, value_ring=None):
        self.__scheme = scheme
        if value_ring is None:
            value_ring = scheme.base_ring()
        self.__value_ring = value_ring

    def _repr_(self):
        return "Group of Divisors on %s"%self.__scheme

    def __cmp__(self, right):
        if not isinstance(right, DivisorGroup):
            return -1
        return cmp(self.__scheme, right.__scheme)

    def scheme(self):
        return self.__scheme

class DivisorGroup_curve(DivisorGroup):
    def __call__(self, v):
        return Divisor_curve(v)

