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
from sage.structure.formal_sum import FormalSums_generic
from sage.rings.integer_ring import ZZ
import weakref

# we cache DivisorGroups by (scheme, value_ring) pairs
_DivisorGroup_cache = {}

def DivisorGroup(scheme, value_ring=None):
    if value_ring is None:
        value_ring = scheme.base_ring()
    key = (scheme, value_ring)
    if _DivisorGroup_cache.has_key(key):
        DG = _DivisorGroup_cache[key]
        if not DG() is None:
            return DG()

    from sage.schemes.plane_curves.curve import Curve_generic
    if isinstance(scheme, Curve_generic):
        DG = DivisorGroup_curve(scheme, value_ring)
    else:
        DG = DivisorGroup_generic(scheme, value_ring)
    _DivisorGroup_cache[key] = weakref.ref(DG)
    return DG

class DivisorGroup_generic(FormalSums_generic):
    def __init__(self, scheme, value_ring):
        FormalSums_generic.__init__(self)
        self.__scheme = scheme
        self.__value_ring = value_ring

    def _repr_(self):
        return "Group of Divisors on %s"%self.__scheme

    def _coerce_impl(self, x):
        if x.parent() is self:
            return x
        raise TypeError

    def __cmp__(self, right):
        if not isinstance(right, DivisorGroup_generic):
            return -1
        return cmp(self.__scheme, right.__scheme)

    def scheme(self):
        return self.__scheme

    def _an_element_impl(self):
        return self.__scheme.divisor([], check=False, reduce=False)

class DivisorGroup_curve(DivisorGroup_generic):
    def __call__(self, v):
        return Divisor_curve(v)
