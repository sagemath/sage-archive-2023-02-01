r"""
Infinity Ring

The infinity ``ring'' is the set of two elements
\begin{verbatim}
        * infinity
        * A number less than infinity
\end{verbatim}

The rules for arithmetic are that the infinity ring does not canonically
coerce to any other ring, and all other rings canonically coerce to
the infinity ring, sending all elements to the single element
``a number less than infinity'' of the infinity ring.
Arithmetic and comparisons then takes place in the infinity ring,
where all arithmetic operations that are well defined are defined.

EXAMPLES:
We fetch the infinity ring and create some elements:
    sage: P = InfinityRing; P
    The Infinity Ring
    sage: P(5)
    A number less than infinity
    sage: P.ngens()
    1
    sage: oo = P.0; oo
    Infinity

We compare finite numbers with infinity.
    sage: 5 < oo
    True
    sage: 5 > oo
    False
    sage: oo < 5
    False
    sage: oo > 5
    True

We do arithmetic.
    sage: oo + 5
    Infinity

Note that many operations are not defined, since the result is
not well defined.
    sage: oo/0
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '/': 'The Infinity Ring' and 'Integer Ring'

What happened above is that 0 is canonically coerced to
"a number less than infinity" in the infinity ring, and the quotient
is then not well defined.

    sage: 0/oo
    A number less than infinity
    sage: oo * 0
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '*': 'The Infinity Ring' and 'The Infinity Ring'
    sage: oo/oo
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '/': 'The Infinity Ring' and 'The Infinity Ring'
"""

from ring_element import RingElement
from sage.rings.ring import Ring
from sage.structure.element import RingElement, InfinityElement
from sage.structure.parent_gens import ParentWithGens

class InfinityRing_class(Ring):
    def __init__(self):
        ParentWithGens.__init__(self, self, names=('oo',), normalize=False)

    def ngens(self):
        return 1

    def gen(self, n=0):
        try:
            return self._gen
        except AttributeError:
            self._gen = Infinity(self)
        return self._gen

    def less_than_infinity(self):
        try:
            return self._less_than_infinity
        except AttributeError:
            self._less_than_infinity = LessThanInfinity(self)
            return self._less_than_infinity

    def gens(self):
        return [self.gen()]

    def _repr_(self):
        return "The Infinity Ring"

    def __cmp__(self, right):
        if isinstance(right, InfinityRing_class):
            return 0
        return cmp(type(self), type(right))

    def __call__(self, x):
        if isinstance(x, InfinityElement):
            if x.parent() is self:
                return x
            else:
                return self.gen()
        elif isinstance(x, RingElement) or isinstance(x, (int,long,float,complex)):
            return self.less_than_infinity()
        else:
            raise TypeError

    def _coerce_impl(self, x):
        if isinstance(x, InfinityElement):
            x = Infinity()
            x._set_parent(self)
            return x
        elif isinstance(x, RingElement):
            return less_than_infinity
        else:
            raise TypeError

InfinityRing = InfinityRing_class()

class LessThanInfinity(RingElement):
    def __init__(self, parent=InfinityRing):
        RingElement.__init__(self, parent)

    def _repr_(self):
        return "A number less than infinity"

    def _latex_(self):
        return "(<\\infty)"

    def _add_(self, other):
        if isinstance(other, Infinity):
            return other
        return self

    def _sub_(self, other):
        if isinstance(other, Infinity):
            return other
        return self

    def _mul_(self, other):
        if isinstance(other, Infinity):
            raise TypeError, "oo times number < oo not defined"
        return self

    def _div_(self, other):
        if isinstance(other, Infinity):
            return self
        raise TypeError, "quotient of oo by number < oo not defined"

    def __cmp__(self, other):
        if isinstance(other, Infinity):
            return -1
        return 0


class Infinity(InfinityElement):
    def __init__(self, parent=InfinityRing):
        InfinityElement.__init__(self, parent)

    def _repr_(self):
        return "Infinity"

    def lcm(self, x):
        """
        Return the least common multiple of oo and x, which
        is by definition oo unless x is 0.

        EXAMPLES:
            sage: oo.lcm(0)
            0
            sage: oo.lcm(oo)
            Infinity
            sage: oo.lcm(10)
            Infinity
        """
        if x == 0:
            return x
        else:
            return self

    def _latex_(self):
        return "\\infty"

    def _add_(self, other):
        return self

    def _sub_(self, other):
        if not isinstance(other, Infinity):
            return self
        raise TypeError, "oo - oo not defined"

    def _mul_(self, other):
        if isinstance(other, Infinity):
            return self
        raise TypeError, "oo times smaller number not defined"

    def __cmp__(self, other):
        if isinstance(other, Infinity):
            return 0
        return 1


infinity = InfinityRing.gen(0)
less_than_infinity = InfinityRing.less_than_infinity()

def is_Infinity(x):
    return isinstance(x, Infinity)

