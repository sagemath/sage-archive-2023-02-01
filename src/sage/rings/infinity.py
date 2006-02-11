r"""
Infinity
"""

from infinity_element import InfinityElement
from ring_element import RingElement
from sage.misc.latex import latex

class Infinity(InfinityElement):
    def __init__(self):
        InfinityElement.__init__(self, None)

    #def __reduce__(self):
    #    """
    #    sage: loads(oo.dumps())
    #    Infinity
    #    """
    #    return Infinity, tuple([])

    def __repr__(self):
        return "Infinity"

    def _latex_(self):
        return "\\infty"

    def __add__(self, other):
        if isinstance(other, (RingElement, Infinity, int, long)):
            return self
        raise TypeError, "%s + %s not defined."%(self, other)

    def __radd__(self, left):
        return self+left

    def __sub__(self, other):
        if isinstance(other, (RingElement, int, long)) and not isinstance(other, InfinityElement):
            return self
        raise TypeError, "%s - %s not defined."%(self, other)

    def __rsub__(self, other):
        if isinstance(other, (RingElement, int, long)) and not isinstance(other, InfinityElement):
            return self
        raise TypeError, "%s - %s not defined."%(other, self)

    def __mul__(self, other):
        if isinstance(other, Infinity):
            return self
        if isinstance(other, (RingElement, int, long)):
            if other == 0:
                return other
            return self
        raise TypeError, "%s * %s not defined."%(self*other)

    def __rmul__(self, left):
        return self*left

    def __div__(self, other):
        if isinstance(other, (RingElement, int, long)) and not isinstance(other, InfinityElement):
            if other != 0:
                return self
        raise TypeError, "%s * %s not defined."%(self*other)

    def __cmp__(self, other):
        if isinstance(other, InfinityElement):
            return 0
        return 1  # infinity is bigger than everything but itself


infinity = Infinity()

def is_Infinity(x):
    return isinstance(x, Infinity)

