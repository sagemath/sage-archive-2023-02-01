"""
Interval arithmetic
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@ucsd.edu>
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

import math
import operator

from coerce import bin_op
import sage.rings.ring as ring
import integer
import rational
import real_field
import ring_element

def _adapt(x):
    if isinstance(x, Interval):
        return x
    if isinstance(x, (int, long, float, integer.Integer,
                      rational.Rational, real_field.RealNumber)):
        return Interval(x)
    raise TypeError, "cannot construct interval from %s"%x

_obj = {}
class _uniq(object):
    def __new__(cls):
        if _obj.has_key(0):
            return _obj[0]
        O = object.__new__(cls)
        _obj[0] = O
        return O

class IntervalRing(ring.Ring, _uniq):
    """
    EXAMPLES:
        sage: R = IntervalRing(); R
        Float interval arithmetic pseudoring.
        sage: loads(R.dumps()) == R
        True
    """
    def _repr_(self):
        return 'Float interval arithmetic pseudoring.'

    def __call__(self, x, y=None):
        if y is None:
            return _adapt(x)
        else:
            return Interval(x,y)

    def _coerce_(self, x):
        if isinstance(x, (int, long, float, Interval, integer.Integer,
                          rational.Rational, real_field_element.RealNumber)):
            return self(x)
        raise TypeError

    def is_atomic_repr(self):
        return True

    def is_field(self):
        return True

    def characteristic(self):
        return 0

    def random(self, bound=0):
        return Interval(0)

    def __cmp__(self, other):
        if isinstance(other, IntervalRing):
            return 0
        return -1

    def zeta(self):
        return Interval(-1)

_inst = IntervalRing()  # make sure there is at least one instance.


class Interval(ring_element.RingElement):
    """
    EXAMPLES:
        sage: R = IntervalRing()
        sage: a = R(0,1); a
        [0.0, 1.0]
        sage: loads(a.dumps()) == a
        True
    """
    def __init__(self, min, max=None):
        if max == None: max = min
        self.__min, self.__max = float(min), float(max)
        assert self.__min <= self.__max
        ring_element.RingElement.__init__(self, _inst)

    def _repr_(self):
        return "[%s, %s]"%(self.__min, self.__max)

    def _latex_(self):
        return "[%s, %s]"%(self.__min, self.__max)

    def min(self):
        return self.__min
    lower = min

    def max(self):
        return self.__max
    upper = max

    def __add__(self, other):
        if not isinstance(other, Interval):
            return bin_op(self, other, operator.add)
        return Interval(self.__min+other.__min, self.__max + other.__max)

    def __sub__(self, other):
        if not isinstance(other, Interval):
            return bin_op(self, other, operator.sub)
        return Interval(self.__min+- other.__max, self.__max - other.__min)

    def __neg__(self):
        return Interval(-self.__max, -self.__min)

    def __mul__(self, other):
        if not isinstance(other, Interval):
            return bin_op(self, other, operator.mul)
        xl,xu,yl,yu = self.__min, self.__max, other.__min, other.__max
        return Interval(min(xl*yl,xl*yu,xu*yl,xu*yu), \
                        max(xl*yl,xl*yu,xu*yl,xu*yu))

    def __invert__(self):
        if 0 in self:
            raise ZeroDivisionError, "cannot invert interval %s"%self
        return Interval(1/self.__max, 1/self.__min)

    def __div__(self, other):
        if not isinstance(other, Interval):
            return bin_op(self, other, operator.div)
        return self * ~other

    def __contains__(self, x):
        x = float(x)
        return self.__min <= x and x <= self.__max

    def __cmp__(self, other):
        try:
            other = _adapt(other)
        except TypeError:
            return -1
        if self.__max < other.__min:
            return -1
        elif self.__min > other.__max:
            return 1
        elif self.__min == other.__min and self.__max == other.__max:
            return 0
        else:
            return -1
        #   raise ArithmeticError, "Unable to compare intervals %s and %s"%(self, other)

    def sqrt(self):
        return Interval(math.sqrt(self.__min), math.sqrt(self.__max))

    def length(self):
        return self.__max - self.__min

    def __int__(self):
        """
        If there is a unique integer in the interval, return it.
        Otherwise raise a ValueError exception.
        """
        x = math.ceil(self.__min)
        y = math.floor(self.__max)
        if x != y:
            raise ValueError, "Cannot coerce to int because there is not a unique integer in the interval %s"%self
        return int(x)

    def _integer_(self):
        """
        Used for coercion to the integers.

        EXAMPLES:
            sage: R = RealField(53)
            sage: I = IntervalRing()
            sage: a = I(R('1.6'), R('2.7'))
            sage: ZZ(a)
            2
            sage: a = I(R('2.1'), R('2.7'))
            sage: ZZ(a)
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce to int because there is not a unique integer in the interval [2.1, 2.7]
            sage: a = I(R('2.1'), R('5.7'))
            sage: ZZ(a)
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce to int because there is not a unique integer in the interval [2.1, 5.7]
        """
        return integer.Integer(int(self))

    def is_int(self):
        try:
            x = int(self)
        except ValueError:
            return False, None
        return True, x

