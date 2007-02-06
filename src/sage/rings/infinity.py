r"""
Infinity Ring

The infinity ``ring'' is the set of two elements
\begin{verbatim}
        * infinity
        * A number less than infinity
\end{verbatim}

This isn't really a ring, but a formal construction that is incredibly
useful in much of the implementation of SAGE.

The rules for arithmetic are that the infinity ring does not canonically
coerce to any other ring, and all other rings canonically coerce to
the infinity ring, sending all elements to the single element
``a number less than infinity'' of the infinity ring.
Arithmetic and comparisons then takes place in the infinity ring,
where all arithmetic operations that are well defined are defined.

EXAMPLES:
We fetch the infinity ring and create some elements:
    sage: P = InfinityRing(); P
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
    ArithmeticError: quotient of oo by number < oo not defined

What happened above is that 0 is canonically coerced to
"a number less than infinity" in the infinity ring, and the quotient
is then not well defined.

    sage: 0/oo
    A number less than infinity
    sage: oo * 0
    Traceback (most recent call last):
    ...
    ArithmeticError: oo times smaller number not defined

    sage: oo/oo
    Traceback (most recent call last):
    ...
    ArithmeticError: oo / oo not defined


Saving and loading:
    sage: R = loads(dumps(InfinityRing())); R
    The Infinity Ring
    sage: R == InfinityRing()
    True
    sage: R is InfinityRing()
    False
"""

from ring_element import RingElement
from sage.rings.ring import Ring
from sage.structure.element import RingElement, InfinityElement
from sage.structure.parent_gens import ParentWithGens

class InfinityRing_class(Ring):
    """
    The infinity "ring", which contains oo and one formal quantity
    that is "less than infinity".
    """
    def __init__(self):
        ParentWithGens.__init__(self, self, names=('oo',), normalize=False)

    def ngens(self):
        """
        Return the number of generators (1) of the infinity ring.

        EXAMPLES:
            sage: InfinityRing().ngens()
            1
        """
        return 1

    def gen(self, n=0):
        """
        Return the "generator" of the infinity ring.  By convention this
        is infinity.

        EXAMPLES:
            sage: InfinityRing().gen()
            Infinity
        """
        try:
            return self._gen
        except AttributeError:
            self._gen = Infinity(self)
        return self._gen

    def less_than_infinity(self):
        """
        EXAMPLES:
            sage: InfinityRing().less_than_infinity()
            A number less than infinity
        """
        try:
            return self._less_than_infinity
        except AttributeError:
            self._less_than_infinity = LessThanInfinity(self)
            return self._less_than_infinity

    def gens(self):
        """
        EXAMPLES:
            sage: InfinityRing().gens()
            [Infinity]
        """
        return [self.gen()]

    def _repr_(self):
        """
        EXAMPLES:
            sage: InfinityRing()
            The Infinity Ring
        """
        return "The Infinity Ring"

    def __cmp__(self, right):
        """
        The only ring that is equal to the infinity ring is another copy of the
        infinity ring.

        Other rings compare based on their underlying types.

        EXAMPLES:
            sage: InfinityRing() == InfinityRing()
            True
            sage: InfinityRing() == loads(dumps(InfinityRing()))
            True
            sage: InfinityRing() == ZZ
            False
        """
        if isinstance(right, InfinityRing_class):
            return 0
        return cmp(type(self), type(right))

    def __call__(self, x):
        """
        Coerce x into the infinity ring.

        If x is infinity, then x coerce to the infinity element.  All other
        elements of rings coerce to the "less than infinity" element.

        EXAMPLES:
            sage: R = InfinityRing()
            sage: R(5)
            A number less than infinity
            sage: R(oo)
            Infinity
            sage: R(RDF(1.459))
            A number less than infinity
            sage: R(RDF(1)/RDF(0))       # real double field
            Infinity
            sage: R(RR(1)/RR(0))   # mpr real field
            Infinity
        """
        if isinstance(x, InfinityElement):
            if x.parent() is self:
                return x
            else:
                return self.gen()
        elif isinstance(x, (int,long,float,complex)):
            return self.less_than_infinity()
        elif isinstance(x, RingElement):
            if hasattr(x, 'is_infinity') and x.is_infinity():
                    return self.gen()
            return self.less_than_infinity()
        else:
            raise TypeError, 'no coercion of non-ring element to infinity ring'

    def _coerce_impl(self, x):
        """
        Canonical coercion of elements into self.

        EXAMPLES:
            sage: R = InfinityRing()
            sage: R._coerce_(oo)
            Infinity
            sage: R._coerce_(5)
            A number less than infinity
            sage: R._coerce_('hello')
            Traceback (most recent call last):
            ...
            TypeError
        """
        if isinstance(x, InfinityElement):
            x = Infinity()
            x._set_parent(self)
            return x
        elif isinstance(x, RingElement):
            return self.less_than_infinity()
        else:
            raise TypeError

theInfinityRing = InfinityRing_class()
def InfinityRing():
    r"""
    Return the infinity ring.

    The infinity ``ring'' is the set of two elements "infinity" and
    "a number less than infinity".

    EXAMPLES:
        sage: InfinityRing()
        The Infinity Ring
    """
    return theInfinityRing

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
            raise ArithmeticError, "oo times number < oo not defined"
        return self

    def _div_(self, other):
        if isinstance(other, Infinity):
            return self
        raise ArithmeticError, "quotient of numbers < oo not defined"

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
        raise ArithmeticError, "oo - oo not defined"

    def _mul_(self, other):
        if isinstance(other, Infinity):
            return self
        raise ArithmeticError, "oo times smaller number not defined"

    def _div_(self, other):
        if isinstance(other, Infinity):
            raise ArithmeticError, "oo / oo not defined"
        raise ArithmeticError, "quotient of oo by number < oo not defined"

    def __cmp__(self, other):
        if isinstance(other, Infinity):
            return 0
        return 1


infinity = theInfinityRing.gen(0)
less_than_infinity = theInfinityRing.less_than_infinity()

def is_Infinity(x):
    return isinstance(x, Infinity)

