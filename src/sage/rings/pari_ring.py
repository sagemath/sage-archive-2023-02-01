"""
Ring of pari objects
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
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
    def __init__(self):
        ring.Ring.__init__(self, self)
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

    def random_element(self, bound=0):
        return Pari(0)

    def random(self, bound=0):
        """
        Deprecated.  Use self.random_element() instead.
        """
        raise NotImplementedError, "Deprecated: use random_element() instead"

    def __cmp__(self, other):
        return cmp(type(self),type(other))

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

    def _add_(self, other):
        return Pari(self.__x + other.__x)

    def _sub_(self, other):
        return Pari(self.__x - other.__x)

    def _mul_(self, other):
        return Pari(self.__x * other.__x)

    def _div_(self, other):
        return self.__x * (~other.__x)

    def __neg__(self):
        return Pari(-self.__x)

    def __pow__(self, other):
        if not isinstance(other, Pari):
            return bin_op(self, other, operator.pow)
        return Pari(self.__x ** other.__x)

    def __invert__(self):
        return Pari(~self.__x)

    def __cmp__(self, other):
        return cmp(self.__x, other.__x)

    def __int__(self):
        return int(self.__x)



