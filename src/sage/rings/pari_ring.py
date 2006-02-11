"""
Ring of pari objects
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator

import sage.libs.pari.all as pari
import sage.rings.ring as ring
import ring_element
from coerce import bin_op

_obj = {}
class _uniq(object):
    def __new__(cls):
        if _obj.has_key(0):
            return _obj[0]
        O = object.__new__(cls)
        _obj[0] = O
        return O

class PariRing(ring.Ring, _uniq):
    """
    EXAMPLES:
        sage: R = PariRing(); R
        Pseudoring of all PARI objects.
        sage: loads(R.dumps()) == R
        True
    """
    def __repr__(self):
        return 'Pseudoring of all PARI objects.'

    def __call__(self, x):
        if isinstance(x, Pari):
            return x
        return Pari(x)

    def is_atomic_repr(self):
        return False

    def is_field(self):
        return False

    def characteristic(self):
        raise RuntimeError, "Not defined."
        #return 0

    def random(self, bound=0):
        return Pari(0)

    def __cmp__(self, other):
        if isinstance(other, PariRing):
            return 0
        return -1

    def zeta(self):
        return Pari(-1)


_inst = PariRing()

class Pari(ring_element.RingElement):
    """
    Element of Pari pseudo-ring.
    """
    def __init__(self, x):
        """
        EXAMPLES:
            sage: R = PariRing()
            sage: f = R('x^3 + 1/2')
            sage: f
            x^3 + 1/2
            sage: type(f)
            <class 'sage.rings.pari_ring.Pari'>
            sage: loads(f.dumps()) == f
            True
        """
        self.__x = pari.pari(x)

    def parent(self):
        return _inst

    def __repr__(self):
        return str(self.__x)

    def __add__(self, other):
        if not isinstance(other, Pari):
            return bin_op(self, other, operator.add)
        return Pari(self.__x + other.__x)

    def __radd__(self, other):
        if not isinstance(other, Pari):
            return bin_op(self, other, operator.add)
        return Pari(other.__x + self.__x )

    def __sub__(self, other):
        if not isinstance(other, Pari):
            return bin_op(self, other, operator.sub)
        return Pari(self.__x - other.__x)

    def __neg__(self):
        return Pari(-self.__x)

    def __mul__(self, other):
        if not isinstance(other, Pari):
            return bin_op(self, other, operator.mul)
        return Pari(self.__x * other.__x)

    def __pow__(self, other):
        if not isinstance(other, Pari):
            return bin_op(self, other, operator.pow)
        return Pari(self.__x ** other.__x)

    def __invert__(self):
        return Pari(~self.__x)

    def __div__(self, other):
        if not isinstance(other, Pari):
            return bin_op(self, other, operator.div)
        return self.__x * (~other.__x)

    def __cmp__(self, other):
        try:
            if not isinstance(other, Pari):
                other = Pari(other)
        except TypeError:
            return -1
        if self.__x < other.__x:
            return -1
        elif self.__x > other.__x:
            return 1
        return 0

    def __int__(self):
        return int(self.__x)



