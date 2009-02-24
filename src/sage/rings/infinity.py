r"""
Infinity Rings

The unsigned infinity "ring" is the set of two elements

::

            * infinity
            * A number less than infinity

The rules for arithmetic are that the unsigned infinity ring does
not canonically coerce to any other ring, and all other rings
canonically coerce to the unsigned infinity ring, sending all
elements to the single element "a number less than infinity" of the
unsigned infinity ring. Arithmetic and comparisons then take place
in the unsigned infinity ring, where all arithmetic operations that
are well-defined are defined.

The infinity "ring" is the set of five elements

::

            * plus infinity
            * a positive finite element
            * zero
            * a negative finite element
            * negative infinity

The infinity ring coerces to the unsigned infinity ring, sending
the infinite elements to infinity and the non-infinite elements to
"a number less than infinity." Any ordered ring coerces to the
infinity ring in the obvious way.

Note: the shorthand oo is predefined in Sage to be the same as
+Infinity in the infinity ring. It is considered equal to, but not
the same as Infinity in the UnsignedInfinityRing::

    sage: oo
    +Infinity
    sage: oo is InfinityRing.0
    True
    sage: oo is UnsignedInfinityRing.0
    False
    sage: oo == UnsignedInfinityRing.0
    True

EXAMPLES: We fetch the unsigned infinity ring and create some
elements::

    sage: P = UnsignedInfinityRing; P
    The Unsigned Infinity Ring
    sage: P(5)
    A number less than infinity
    sage: P.ngens()
    1
    sage: unsigned_oo = P.0; unsigned_oo
    Infinity

We compare finite numbers with infinity::

    sage: 5 < unsigned_oo
    True
    sage: 5 > unsigned_oo
    False
    sage: unsigned_oo < 5
    False
    sage: unsigned_oo > 5
    True

We do arithmetic::

    sage: unsigned_oo + 5
    Infinity

We make 1 / unsigned_oo return the integer 0 so that arithmetic of
the following type works::

    sage: (1/unsigned_oo) + 2
    2
    sage: 32/5 - (2.439/unsigned_oo)
    32/5

Note that many operations are not defined, since the result is not
well-defined::

    sage: unsigned_oo/0
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '/': 'The Unsigned Infinity Ring' and 'Integer Ring'

What happened above is that 0 is canonically coerced to "a number
less than infinity" in the unsigned infinity ring, and the quotient
is then not well-defined.

::

    sage: 0/unsigned_oo
    0
    sage: unsigned_oo * 0
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '*': 'The Unsigned Infinity Ring' and 'Integer Ring'
    sage: unsigned_oo/unsigned_oo
    Traceback (most recent call last):
    ...
    TypeError: infinity 'ring' has no fraction field

In the infinity ring, we can negate infinity, multiply positive
numbers by infinity, etc.

::

    sage: P = InfinityRing; P
    The Infinity Ring
    sage: P(5)
    A positive finite number

The symbol oo is predefined as a shorthand for +Infinity::

    sage: oo
    +Infinity

We compare finite and infinite elements::

    sage: 5 < oo
    True
    sage: P(-5) < P(5)
    True
    sage: P(2) < P(3)
    False
    sage: -oo < oo
    True

We can do more arithmetic than in the unsigned infinity ring::

    sage: 2 * oo
    +Infinity
    sage: -2 * oo
    -Infinity
    sage: 1 - oo
    -Infinity
    sage: 1 / oo
    0
    sage: -1 / oo
    0

We make 1 / oo and 1 / -oo return the integer 0 instead of the
infinity ring Zero so that arithmetic of the following type works::

    sage: (1/oo) + 2
    2
    sage: 32/5 - (2.439/-oo)
    32/5

If we try to subtract infinities or multiply infinity by zero we
still get an error::

    sage: oo - oo
    Traceback (most recent call last):
    ...
    SignError: cannot add infinity to minus infinity
    sage: 0 * oo
    Traceback (most recent call last):
    ...
    SignError: cannot multiply infinity by zero
    sage: P(2) + P(-3)
    Traceback (most recent call last):
    ...
    SignError: cannot add positive finite value to negative finite value

TESTS::

    sage: P = InfinityRing
    sage: P == loads(dumps(P))
    True

::

    sage: P(2) == loads(dumps(P(2)))
    True

The following is assumed in a lot of code (i.e., "is" is used for
testing whether something is infinity), so make sure it is
satisfied::

    sage: loads(dumps(infinity)) is infinity
    True
"""

from sage.rings.ring_element import RingElement
from sage.rings.ring import Ring
from sage.structure.element import RingElement, InfinityElement, PlusInfinityElement, MinusInfinityElement
from sage.structure.parent_gens import ParentWithGens
#import sage.rings.real_double
#import sage.rings.real_mpfr
import sage.rings.integer
import sage.rings.rational

from sage.rings.integer_ring import ZZ

_obj = {}
class _uniq0(object):
    def __new__(cls):
        if _obj.has_key(0):
            return _obj[0]
        O = Ring.__new__(cls)
        _obj[0] = O
        return O
class _uniq1(object):
    def __new__(cls):
        if _obj.has_key(1):
            return _obj[1]
        O = InfinityElement.__new__(cls)
        _obj[1] = O
        return O
class _uniq2(object):
    def __new__(cls):
        if _obj.has_key(2):
            return _obj[2]
        O = Ring.__new__(cls)
        _obj[2] = O
        return O
class _uniq3(object):
    def __new__(cls):
        if _obj.has_key(3):
            return _obj[3]
        O = MinusInfinityElement.__new__(cls)
        _obj[3] = O
        return O
class _uniq4(object):
    def __new__(cls):
        if _obj.has_key(4):
            return _obj[4]
        O = PlusInfinityElement.__new__(cls)
        _obj[4] = O
        return O

class UnsignedInfinityRing_class(_uniq0, Ring):
    def __init__(self):
        ParentWithGens.__init__(self, self, names=('oo',), normalize=False)

    def ngens(self):
        return 1

    def fraction_field(self):
        raise TypeError, "infinity 'ring' has no fraction field"

    def gen(self, n=0):
        try:
            return self._gen
        except AttributeError:
            self._gen = UnsignedInfinity()
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
        return "The Unsigned Infinity Ring"

    def __cmp__(self, right):
        if isinstance(right, UnsignedInfinityRing_class):
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
            x = UnsignedInfinity()
            x._set_parent(self)
            return x
        elif isinstance(x, RingElement):
            return less_than_infinity
        else:
            raise TypeError

UnsignedInfinityRing = UnsignedInfinityRing_class()

class LessThanInfinity(RingElement):
    def __init__(self, parent=UnsignedInfinityRing):
        RingElement.__init__(self, parent)

    def _repr_(self):
        return "A number less than infinity"

    def _latex_(self):
        return "(<\\infty)"

    def _add_(self, other):
        if isinstance(other, UnsignedInfinity):
            return other
        return self

    def _sub_(self, other):
        if isinstance(other, UnsignedInfinity):
            return other
        return self

    def _mul_(self, other):
        if isinstance(other, UnsignedInfinity):
            raise TypeError, "oo times number < oo not defined"
        return self

    def _div_(self, other):
        if isinstance(other, UnsignedInfinity):
            return ZZ(0)
        raise TypeError, "quotient of oo by number < oo not defined"

    def __cmp__(self, other):
        if isinstance(other, UnsignedInfinity):
            return -1
        return 0


class UnsignedInfinity(_uniq1, InfinityElement):
    def __init__(self):
        InfinityElement.__init__(self, UnsignedInfinityRing)

    def _repr_(self):
        return "Infinity"

    def _maxima_init_(self):
        return "inf"

    def lcm(self, x):
        """
        Return the least common multiple of oo and x, which is by
        definition oo unless x is 0.

        EXAMPLES::

            sage: oo = UnsignedInfinityRing.gen(0)
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
        if not isinstance(other, UnsignedInfinity):
            return self
        raise TypeError, "oo - oo not defined"

    def _mul_(self, other):
        if isinstance(other, UnsignedInfinity):
            return self
        raise TypeError, "oo times smaller number not defined"

    def __cmp__(self, other):
        if isinstance(other, UnsignedInfinity):
            return 0
        return 1


unsigned_infinity = UnsignedInfinityRing.gen(0)
less_than_infinity = UnsignedInfinityRing.less_than_infinity()

def is_Infinite(x):
    return isinstance(x, InfinityElement)

class SignError(Exception):
    pass

class InfinityRing_class(_uniq2, Ring):
    def __init__(self):
        ParentWithGens.__init__(self, self, names=('oo',), normalize=False)

    def fraction_field(self):
        raise TypeError, "infinity 'ring' has no fraction field"

    def ngens(self):
        return 2

    def gen(self, n=0):
        try:
            if n == 0:
                return self._gen0
            elif n == 1:
                return self._gen1
            else:
                raise IndexError, "n must be 0 or 1"
        except AttributeError:
            if n == 0:
                self._gen0 = PlusInfinity()
                return self._gen0
            elif n == 1:
                self._gen1 = MinusInfinity()
                return self._gen1

    def gens(self):
        return [self.gen(0), self.gen(1)]

    def _repr_(self):
        return "The Infinity Ring"

    def __cmp__(self, right):
        if isinstance(right, InfinityRing_class):
            return 0
        return cmp(type(self), type(right))

    def __call__(self, x):
        if isinstance(x, PlusInfinityElement):
            return self.gen(0)
        elif isinstance(x, MinusInfinityElement):
            return self.gen(1)
        elif isinstance(x, InfinityElement):
            return self.gen(0)
        elif isinstance(x, (sage.rings.integer.Integer, sage.rings.rational.Rational, sage.rings.real_double.RealDoubleElement, sage.rings.real_mpfr.RealNumber)) or isinstance(x, (int,long,float)):
            if x < 0:
                return FiniteNumber(self, -1)
            elif x > 0:
                return FiniteNumber(self, 1)
            else:
                return FiniteNumber(self, 0)
        else:
            raise TypeError

    def _coerce_impl(self, x):
        if isinstance(x, PlusInfinityElement):
            if x.parent() is self:
                return x
            else:
                return self.gen(0)
        elif isinstance(x, MinusInfinityElement):
            if x.parent() is self:
                return x
            else:
                return self.gen(1)
        elif isinstance(x, (sage.rings.integer.Integer, sage.rings.rational.Rational, sage.rings.real_double.RealDoubleElement, sage.rings.real_mpfr.RealNumber)) or isinstance(x, (int,long,float)):
            if x < 0:
                return FiniteNumber(self, -1)
            elif x > 0:
                return FiniteNumber(self, 1)
            else:
                return FiniteNumber(self, 0)
        else:
            raise TypeError

class FiniteNumber(RingElement):
    def __init__(self, parent, x):
        RingElement.__init__(self, parent)
        self.value = x

    def __cmp__(self, other):
        if isinstance(other, PlusInfinity):
            return -1
        if isinstance(other, MinusInfinity):
            return 1
        return self.value.__cmp__(other.value)

    def _add_(self, other):
        if isinstance(other, InfinityElement):
            return other
        if self.value * other.value < 0:
            raise SignError, "cannot add positive finite value to negative finite value"
        return FiniteNumber(self.parent(), self.value)

    def _mul_(self, other):
        if other.is_zero():
            if isinstance(self, InfinityElement):
                raise SignError, "cannot multiply infinity by zero"
            return ZZ(0)
        if self.value < 0:
            if isinstance(other, InfinityElement):
                return -other
            return FiniteNumber(self.parent(), self.value * other.value)
        if self.value > 0:
            if isinstance(other, InfinityElement):
                return other
            return FiniteNumber(self.parent(), self.value * other.value)
        if self.value == 0:
            if isinstance(other, InfinityElement):
                raise SignError, "cannot multiply infinity by zero"
            return ZZ(0)

    def _div_(self, other):
        return self._mul_(other.__invert__())

    def _sub_(self, other):
        return self._add_(other._neg_())

    def __invert__(self):
        if self.value == 0:
            raise ZeroDivisionError, "Cannot divide by zero"
        return self

    def _neg_(self):
        return FiniteNumber(self.parent(), -self.value)

    def _repr_(self):
        if self.value < 0:
            return "A negative finite number"
        if self.value > 0:
            return "A positive finite number"
        return "Zero"

    def _latex_(self):
        return self._repr_()

    def __abs__(self):
        if self.value == 0:
            return FiniteNumber(self.parent(), 0)
        return FiniteNumber(self.parent(), 1)

    def sqrt(self):
        if self._value < 0:
            raise SignError, "cannot take square root of a negative number"
        return self

    def square_root(self):
        if self._value < 0:
            raise SignError, "cannot take square root of a negative number"
        return self

class MinusInfinity(_uniq3, MinusInfinityElement):
    def __init__(self):
        InfinityElement.__init__(self, InfinityRing)

    def __cmp__(self, other):
        if isinstance(other, MinusInfinity):
            return 0
        return -1

    def _repr_(self):
        return "-Infinity"

    def _maxima_init_(self):
        """
        EXAMPLES::

            sage: maxima(-oo)
            minf
        """
        return "minf"

    def _latex_(self):
        return "-\\infty"

    def __float__(self):
        r"""
        Generate a floating-point -infinity.  The printing of
        floating-point -infinity varies across platforms.

        EXAMPLES:
            sage: RDF(-infinity)
            -inf
            sage: float(-infinity) # random
            -inf
            sage: CDF(-infinity)
            -inf
            sage: (-infinity).__float__() # random
            -inf
        """
        # Evidently there is no standard way to generate an infinity
        # in Python (before Python 2.6).
        from sage.rings.all import RR
        return float(RR(self))

    def _add_(self, other):
        if isinstance(other, PlusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _mul_(self, other):
        if other < 0:
            return -self
        if other > 0:
            return self
        raise SignError, "cannot multiply infinity by zero"

    def _sub_(self, other):
        if isinstance(other, MinusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _div_(self, other):
        return self * other.__invert__()

    def _neg_(self):
        return self.parent().gen(0)

    def __abs__(self):
        return self.parent().gen(0)

    def lcm(self, x):
        """
        Return the least common multiple of -oo and x, which is by
        definition oo unless x is 0.

        EXAMPLES::

            sage: moo = InfinityRing.gen(1)
            sage: moo.lcm(0)
            0
            sage: moo.lcm(oo)
            +Infinity
            sage: moo.lcm(10)
            +Infinity
        """
        if x == 0:
            return x
        else:
            return -self

    def sqrt(self):
        raise SignError, "cannot take square root of negative infinity"

    def square_root(self):
        raise SignError, "cannot take square root of negative infinity"

class PlusInfinity(_uniq4, PlusInfinityElement):
    def __init__(self):
        InfinityElement.__init__(self, InfinityRing)

    def __cmp__(self, other):
        if isinstance(other, PlusInfinity):
            return 0
        return 1

    def __repr__(self):
        return "+Infinity"

    def _maxima_init_(self):
        """
        EXAMPLES::

            sage: maxima(oo)
            inf
        """
        return "inf"

    def _latex_(self):
        return "+\\infty"

    def __float__(self):
        r"""
        Generate a floating-point infinity.  The printing of
        floating-point infinity varies across platforms.

        EXAMPLES:
            sage: RDF(infinity)
            inf
            sage: float(infinity) # random
            inf
            sage: CDF(infinity)
            inf
            sage: infinity.__float__() # random
            inf
        """
        # Evidently there is no standard way to generate an infinity
        # in Python (before Python 2.6).
        from sage.rings.all import RR
        return float(RR(self))

    def _add_(self, other):
        if isinstance(other, MinusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _mul_(self, other):
        if other < 0:
            return -self
        if other > 0:
            return self
        raise SignError, "cannot multiply infinity by zero"

    def _sub_(self, other):
        if isinstance(other, PlusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _div_(self, other):
        return self * other.__invert__()

    def _neg_(self):
        return self.parent().gen(1)

    def __abs__(self):
        return self

    def lcm(self, x):
        """
        Return the least common multiple of oo and x, which is by
        definition oo unless x is 0.

        EXAMPLES::

            sage: oo = InfinityRing.gen(0)
            sage: oo.lcm(0)
            0
            sage: oo.lcm(oo)
            +Infinity
            sage: oo.lcm(10)
            +Infinity
        """
        if x == 0:
            return x
        else:
            return self

    def sqrt(self):
        return self

    def square_root(self):
        return self

    def _sympy_(self):
        """
        Converts oo to sympy oo.

        Then you don't have to worry which oo you use, like in these
        examples:

        EXAMPLE::

            sage: import sympy
            sage: bool(oo == sympy.oo) # indirect doctest
            True
            sage: bool(SR(oo) == sympy.oo)
            True
        """
        import sympy
        return sympy.oo

InfinityRing = InfinityRing_class()
infinity = InfinityRing.gen(0)
Infinity = infinity
minus_infinity = InfinityRing.gen(1)
















