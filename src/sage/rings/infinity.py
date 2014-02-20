r"""
Infinity Rings

The unsigned infinity "ring" is the set of two elements

1. infinity

2. A number less than infinity

The rules for arithmetic are that the unsigned infinity ring does
not canonically coerce to any other ring, and all other rings
canonically coerce to the unsigned infinity ring, sending all
elements to the single element "a number less than infinity" of the
unsigned infinity ring. Arithmetic and comparisons then take place
in the unsigned infinity ring, where all arithmetic operations that
are well-defined are defined.

The infinity "ring" is the set of five elements

1. plus infinity

2. a positive finite element

3. zero

4. a negative finite element

5. negative infinity

The infinity ring coerces to the unsigned infinity ring, sending
the infinite elements to infinity and the non-infinite elements to
"a number less than infinity." Any ordered ring coerces to the
infinity ring in the obvious way.

.. NOTE::

    The shorthand ``oo`` is predefined in Sage to be the same as
    ``+Infinity`` in the infinity ring. It is considered equal to, but not
    the same as ``Infinity`` in the
    :class:`UnsignedInfinityRing<UnsignedInfinityRing_class>`.

EXAMPLES:

We fetch the unsigned infinity ring and create some elements::

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

Demonstrating the shorthand ``oo`` versus ``Infinity``::

    sage: oo
    +Infinity
    sage: oo is InfinityRing.0
    True
    sage: oo is UnsignedInfinityRing.0
    False
    sage: oo == UnsignedInfinityRing.0
    True

We do arithmetic::

    sage: unsigned_oo + 5
    Infinity

We make ``1 / unsigned_oo`` return the integer 0 so that arithmetic of
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
    ValueError: unsigned oo times smaller number not defined

What happened above is that 0 is canonically coerced to "a number
less than infinity" in the unsigned infinity ring, and the quotient
is then not well-defined.

::

    sage: 0/unsigned_oo
    0
    sage: unsigned_oo * 0
    Traceback (most recent call last):
    ...
    ValueError: unsigned oo times smaller number not defined
    sage: unsigned_oo/unsigned_oo
    Traceback (most recent call last):
    ...
    ValueError: unsigned oo times smaller number not defined

In the infinity ring, we can negate infinity, multiply positive
numbers by infinity, etc.

::

    sage: P = InfinityRing; P
    The Infinity Ring
    sage: P(5)
    A positive finite number

The symbol ``oo`` is predefined as a shorthand for ``+Infinity``::

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

We make ``1 / oo`` and ``1 / -oo`` return the integer 0 instead of the
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
testing whether something is infinity), so make sure it is satisfied::

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
class _uniq(object):
    def __new__(cls, *args):
        """
        This ensures uniqueness of these objects.

        EXAMPLE::

            sage: sage.rings.infinity.UnsignedInfinityRing_class() is sage.rings.infinity.UnsignedInfinityRing_class()
            True
        """
        if cls in _obj:
            return _obj[cls]
        _obj[cls] = O = cls.__bases__[-1].__new__(cls, *args)
        return O

class AnInfinity:

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: [x._repr_() for x in [unsigned_infinity, oo, -oo]]
            ['Infinity', '+Infinity', '-Infinity']
        """
        return self._sign_char + "Infinity"

    def _giac_init_(self):
        """
        TESTS::

            sage: [x._giac_init_() for x in [unsigned_infinity, oo, -oo]]
            ['infinity', '+infinity', '-infinity']
        """
        return self._sign_char + "infinity"

    def _maxima_init_(self):
        """
        TESTS::

            sage: maxima(-oo)
            minf
            sage: [x._maxima_init_() for x in [unsigned_infinity, oo, -oo]]
            ['inf', 'inf', 'minf']
        """
        if self._sign < 0:
            return 'minf'
        else:
            return 'inf'

    def _pari_(self):
        """
        Convert ``self`` to a Pari object.

        This always raises an exception since Pari does not have
        infinities.

        TESTS::

            sage: pari(-oo) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: cannot convert infinity to Pari
            sage: pari(oo)
            Traceback (most recent call last):
            ...
            TypeError: cannot convert infinity to Pari
        """
        raise TypeError, 'cannot convert infinity to Pari'

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(oo) # indirect doctest
            +\infty
            sage: [x._latex_() for x in [unsigned_infinity, oo, -oo]]
            ['\\infty', '+\\infty', '-\\infty']
        """
        return self._sign_char + "\\infty"

    def __cmp__(self, other):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: oo == oo
            True
            sage: oo < oo
            False
            sage: -oo < oo
            True
            sage: -oo < 3 < oo
            True

            sage: unsigned_infinity == 3
            False
            sage: unsigned_infinity == unsigned_infinity
            True
            sage: unsigned_infinity == oo
            True
        """
        if isinstance(other, AnInfinity):
            return cmp(self._sign, other._sign)
        elif self._sign:
            return self._sign
        else:
            return 1

    def __abs__(self):
        """
        EXAMPLES::

            sage: [abs(x) for x in [UnsignedInfinityRing.gen(), oo, -oo]]
            [Infinity, +Infinity, +Infinity]
        """
        return -self if self._sign < 0 else self

    def _add_(self, other):
        """
        Add ``self`` to ``other``.

        EXAMPLES::

            sage: -oo + -oo # indirect doctest
            -Infinity
            sage: -oo + 3
            -Infinity
            sage: oo + -100
            +Infinity
            sage: oo + -oo
            Traceback (most recent call last):
            ...
            SignError: cannot add infinity to minus infinity

            sage: unsigned_infinity = UnsignedInfinityRing.gen()
            sage: unsigned_infinity + unsigned_infinity
            Infinity
            sage: unsigned_infinity + 88/3
            Infinity
        """
        if isinstance(other, AnInfinity):
            if self._sign != other._sign:
                raise SignError, "cannot add infinity to minus infinity"
        return self

    def _sub_(self, other):
        """
        EXAMPLES::

            sage: -oo - oo # indirect doctest
            -Infinity
            sage: oo - -oo
            +Infinity
            sage: oo - 4
            +Infinity
            sage: -oo - 1
            -Infinity
            sage: oo - oo
            Traceback (most recent call last):
            ...
            SignError: cannot add infinity to minus infinity
            sage: unsigned_infinity - 4
            Infinity
            sage: unsigned_infinity - unsigned_infinity
            Traceback (most recent call last):
            ...
            ValueError: oo - oo not defined
        """
        if isinstance(other, AnInfinity):
            if self._sign == 0:
                raise ValueError, "oo - oo not defined"
            elif self._sign == other._sign:
                raise SignError, "cannot add infinity to minus infinity"
        return self

    def _mul_(self, other):
        """
        EXAMPLES::

            sage: oo * 19 # indirect doctest
            +Infinity
            sage: oo * oo
            +Infinity
            sage: -oo * oo
            -Infinity
            sage: -oo * 4
            -Infinity
            sage: -oo * -2/3
            +Infinity
            sage: -oo * 0
            Traceback (most recent call last):
            ...
            SignError: cannot multiply infinity by zero
        """
        if other < 0:
            return -self
        if other > 0:
            return self
        raise SignError, "cannot multiply infinity by zero"

    def _div_(self, other):
        """
        EXAMPLES::

            sage: 1.5 / oo # indirect doctest
            0
            sage: oo / -4
            -Infinity
            sage: oo / oo
            Traceback (most recent call last):
            ...
            SignError: cannot multiply infinity by zero
        """
        return self * ~other

    def __float__(self):
        r"""
        Generate a floating-point infinity.  The printing of
        floating-point infinity varies across platforms.

        EXAMPLES::

            sage: RDF(infinity)
            +infinity
            sage: float(infinity) # random
            +infinity
            sage: CDF(infinity)
            +infinity
            sage: infinity.__float__() # random
            +infinity

            sage: RDF(-infinity)
            -infinity
            sage: float(-infinity) # random
            -inf
            sage: CDF(-infinity)
            -infinity
            sage: (-infinity).__float__() # random
            -inf

        """
        # Evidently there is no standard way to generate an infinity
        # in Python (before Python 2.6).
        from sage.rings.all import RR
        return float(RR(self))

    def lcm(self, x):
        """
        Return the least common multiple of ``oo`` and ``x``, which
        is by definition oo unless ``x`` is 0.

        EXAMPLES::

            sage: oo.lcm(0)
            0
            sage: oo.lcm(oo)
            +Infinity
            sage: oo.lcm(-oo)
            +Infinity
            sage: oo.lcm(10)
            +Infinity
            sage: (-oo).lcm(10)
            +Infinity
        """
        if x == 0:
            return x
        else:
            return abs(self)


class UnsignedInfinityRing_class(_uniq, Ring):

    def __init__(self):
        """
        Initialize ``self``.

        TESTS::

            sage: sage.rings.infinity.UnsignedInfinityRing_class() is sage.rings.infinity.UnsignedInfinityRing_class() is UnsignedInfinityRing
            True
        """
        ParentWithGens.__init__(self, self, names=('oo',), normalize=False)

    def ngens(self):
        """
        The unsigned infinity ring has one "generator."

        EXAMPLES::

            sage: UnsignedInfinityRing.ngens()
            1
            sage: len(UnsignedInfinityRing.gens())
            1
        """
        return 1

    def fraction_field(self):
        """
        The unsigned infinity ring isn't an integral domain.

        EXAMPLES::

            sage: UnsignedInfinityRing.fraction_field()
            Traceback (most recent call last):
            ...
            TypeError: infinity 'ring' has no fraction field
        """
        raise TypeError, "infinity 'ring' has no fraction field"

    def gen(self, n=0):
        """
        The "generator" of ``self`` is the infinity object.

        EXAMPLES::

            sage: UnsignedInfinityRing.gen()
            Infinity
            sage: UnsignedInfinityRing.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: UnsignedInfinityRing only has one generator
        """
        if n == 0:
            try:
                return self._gen
            except AttributeError:
                self._gen = UnsignedInfinity()
                return self._gen
        else:
            raise IndexError, "UnsignedInfinityRing only has one generator"

    def gens(self):
        """
        The "generator" of ``self`` is the infinity object.

        EXAMPLES::

            sage: UnsignedInfinityRing.gens()
            [Infinity]
        """
        return [self.gen()]

    def less_than_infinity(self):
        """
        This is the element that represents a finite value.

        EXAMPLES::

            sage: UnsignedInfinityRing.less_than_infinity()
            A number less than infinity
            sage: UnsignedInfinityRing(5) is UnsignedInfinityRing.less_than_infinity()
            True
        """
        try:
            return self._less_than_infinity
        except AttributeError:
            self._less_than_infinity = LessThanInfinity(self)
            return self._less_than_infinity

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: UnsignedInfinityRing._repr_()
            'The Unsigned Infinity Ring'
        """
        return "The Unsigned Infinity Ring"

    def __cmp__(self, right):
        """
        Compare ``self`` to ``right``.

        TESTS::

            sage: infinity == UnsignedInfinityRing.gen()
            True
            sage: UnsignedInfinityRing(3) == UnsignedInfinityRing(-19.5)
            True
        """
        if isinstance(right, UnsignedInfinityRing_class):
            return 0
        return cmp(type(self), type(right))

    def _element_constructor_(self, x):
        """
        TESTS::

            sage: UnsignedInfinityRing(2) # indirect doctest
            A number less than infinity
            sage: UnsignedInfinityRing(I)
            A number less than infinity
            sage: UnsignedInfinityRing(infinity)
            Infinity
            sage: UnsignedInfinityRing(float('+inf'))
            Infinity
        """
        if isinstance(x, InfinityElement):
            if x.parent() is self:
                return x
            else:
                return self.gen()
        elif isinstance(x, RingElement) or isinstance(x, (int,long)):
            return self.less_than_infinity()
        elif isinstance(x, (float, complex)):
            if x == float('+inf') or x == float('-inf'):
                return self.gen()
            return self.less_than_infinity()
        else:
            raise TypeError

    def _coerce_map_from_(self, R):
        """
        EXAMPLES::

            sage: UnsignedInfinityRing.has_coerce_map_from(int) # indirect doctest
            True
            sage: UnsignedInfinityRing.has_coerce_map_from(CC)
            True
            sage: UnsignedInfinityRing.has_coerce_map_from(QuadraticField(-163, 'a'))
            True
            sage: UnsignedInfinityRing.has_coerce_map_from(QQ^3)
            False
            sage: UnsignedInfinityRing.has_coerce_map_from(SymmetricGroup(13))
            False
        """
        return isinstance(R, Ring) or R in (int, long, float, complex)

UnsignedInfinityRing = UnsignedInfinityRing_class()

class LessThanInfinity(_uniq, RingElement):
    def __init__(self, parent=UnsignedInfinityRing):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: sage.rings.infinity.LessThanInfinity() is UnsignedInfinityRing(5)
            True
        """
        RingElement.__init__(self, parent)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: UnsignedInfinityRing(5)._repr_()
            'A number less than infinity'
        """
        return "A number less than infinity"

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: UnsignedInfinityRing(5)._latex_()
            '(<\\infty)'
        """
        return "(<\\infty)"

    def _add_(self, other):
        """
        EXAMPLES::

            sage: UnsignedInfinityRing(5) + UnsignedInfinityRing(-3) # indirect doctest
            A number less than infinity
            sage: UnsignedInfinityRing(5) + unsigned_infinity
            Infinity
        """
        if isinstance(other, UnsignedInfinity):
            return other
        return self

    def _sub_(self, other):
        """
        EXAMPLES::

            sage: UnsignedInfinityRing(5) - UnsignedInfinityRing(-3) # indirect doctest
            A number less than infinity
            sage: UnsignedInfinityRing(5) - unsigned_infinity
            Infinity
        """
        if isinstance(other, UnsignedInfinity):
            return other
        return self

    def _mul_(self, other):
        """
        EXAMPLES::

            sage: UnsignedInfinityRing(4) * UnsignedInfinityRing(-3) # indirect doctest
            A number less than infinity
            sage: 5 * unsigned_infinity
            Traceback (most recent call last):
            ...
            ValueError: oo times number < oo not defined
            sage: unsigned_infinity * unsigned_infinity
            Infinity
        """
        if isinstance(other, UnsignedInfinity):
            raise ValueError, "oo times number < oo not defined"
        return self

    def _div_(self, other):
        """
        Can't eliminate possibility of zero division....

        EXAMPLES::

            sage: UnsignedInfinityRing(2) / UnsignedInfinityRing(5) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: quotient of number < oo by number < oo not defined
            sage: 1 / unsigned_infinity
            0
        """
        if isinstance(other, UnsignedInfinity):
            return ZZ(0)
        raise ValueError, "quotient of number < oo by number < oo not defined"

    def __cmp__(self, other):
        """
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: 1 == unsigned_infinity
            False
        """
        if isinstance(other, UnsignedInfinity):
            return -1
        return 0


class UnsignedInfinity(_uniq, AnInfinity, InfinityElement):

    _sign = 0
    _sign_char = ''

    def __init__(self):
        """
        Initialize ``self``.

        TESTS::

            sage: sage.rings.infinity.UnsignedInfinity() is sage.rings.infinity.UnsignedInfinity() is unsigned_infinity
            True
        """
        InfinityElement.__init__(self, UnsignedInfinityRing)

    def _mul_(self, other):
        """
        Can't rule out an attempt at multiplication by 0.

        EXAMPLES::

            sage: unsigned_infinity * unsigned_infinity # indirect doctest
            Infinity
            sage: unsigned_infinity * 0
            Traceback (most recent call last):
            ...
            ValueError: unsigned oo times smaller number not defined
            sage: unsigned_infinity * 3
            Traceback (most recent call last):
            ...
            ValueError: unsigned oo times smaller number not defined
        """
        if isinstance(other, UnsignedInfinity):
            return self
        raise ValueError, "unsigned oo times smaller number not defined"


unsigned_infinity = UnsignedInfinityRing.gen(0)
less_than_infinity = UnsignedInfinityRing.less_than_infinity()

def is_Infinite(x):
    """
    This is a type check for infinity elements.

    EXAMPLES::

        sage: sage.rings.infinity.is_Infinite(oo)
        True
        sage: sage.rings.infinity.is_Infinite(-oo)
        True
        sage: sage.rings.infinity.is_Infinite(unsigned_infinity)
        True
        sage: sage.rings.infinity.is_Infinite(3)
        False
        sage: sage.rings.infinity.is_Infinite(RR(infinity))
        False
        sage: sage.rings.infinity.is_Infinite(ZZ)
        False
    """
    return isinstance(x, InfinityElement)

class SignError(ArithmeticError):
    """
    Sign error exception.
    """
    pass

class InfinityRing_class(_uniq, Ring):
    def __init__(self):
        """
        Initialize ``self``.

        TEST::

            sage: sage.rings.infinity.InfinityRing_class() is sage.rings.infinity.InfinityRing_class() is InfinityRing
            True
        """
        ParentWithGens.__init__(self, self, names=('oo',), normalize=False)

    def fraction_field(self):
        """
        This isn't really a ring, let alone an integral domain.

        TEST::

            sage: InfinityRing.fraction_field()
            Traceback (most recent call last):
            ...
            TypeError: infinity 'ring' has no fraction field
        """
        raise TypeError, "infinity 'ring' has no fraction field"

    def ngens(self):
        """
        The two generators are plus and minus infinity.

        EXAMPLES::

            sage: InfinityRing.ngens()
            2
            sage: len(InfinityRing.gens())
            2
        """
        return 2

    def gen(self, n=0):
        """
        The two generators are plus and minus infinity.

        EXAMPLES::

            sage: InfinityRing.gen(0)
            +Infinity
            sage: InfinityRing.gen(1)
            -Infinity
            sage: InfinityRing.gen(2)
            Traceback (most recent call last):
            ...
            IndexError: n must be 0 or 1
        """
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
        """
        The two generators are plus and minus infinity.

        EXAMPLES::

            sage: InfinityRing.gens()
            [+Infinity, -Infinity]
        """
        return [self.gen(0), self.gen(1)]

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TEST::

            sage: InfinityRing._repr_()
            'The Infinity Ring'
        """
        return "The Infinity Ring"

    def __cmp__(self, right):
        """
        Compare ``self`` to ``right``.

        TESTS::

            sage: InfinityRing == InfinityRing
            True
            sage: InfinityRing == UnsignedInfinityRing
            False
        """
        if isinstance(right, InfinityRing_class):
            return 0
        return cmp(type(self), type(right))

    def _element_constructor_(self, x):
        """
        TESTS::

            sage: InfinityRing(-oo) # indirect doctest
            -Infinity
            sage: InfinityRing(3)
            A positive finite number
            sage: InfinityRing(-1.5)
            A negative finite number
            sage: [InfinityRing(a) for a in [-2..2]]
            [A negative finite number, A negative finite number, Zero, A positive finite number, A positive finite number]
            sage: K.<a> = QuadraticField(3)
            sage: InfinityRing(a)
            A positive finite number
            sage: InfinityRing(a - 2)
            A negative finite number

        Check that :trac:`14045` is fixed::

            sage: InfinityRing(float('+inf'))
            +Infinity
            sage: InfinityRing(float('-inf'))
            -Infinity
            sage: oo > float('+inf')
            False
            sage: oo == float('+inf')
            True
        """
        if isinstance(x, PlusInfinityElement):
            return self.gen(0)
        elif isinstance(x, MinusInfinityElement):
            return self.gen(1)
        elif isinstance(x, InfinityElement):
            return self.gen(0)
        elif (isinstance(x, (sage.rings.integer.Integer,
                             sage.rings.rational.Rational,
                             sage.rings.real_double.RealDoubleElement,
                             sage.rings.real_mpfr.RealNumber))
                or isinstance(x, (int,long,float))
                or sage.rings.real_mpfr.RealField(sage.rings.real_mpfr.mpfr_prec_min())(x)):
            if x == float('+inf'):
                return self.gen(0)
            elif x == float('-inf'):
                return self.gen(1)
            elif x < 0:
                return FiniteNumber(self, -1)
            elif x > 0:
                return FiniteNumber(self, 1)
            elif x == 0:
                return FiniteNumber(self, 0)
            else:
                raise TypeError
        else:
            raise TypeError

    def _coerce_map_from_(self, R):
        """
        There is a coercion from anything that has a coercion into the reals.

        EXAMPLES::

            sage: InfinityRing.has_coerce_map_from(int) # indirect doctest
            True
            sage: InfinityRing.has_coerce_map_from(AA)
            True
            sage: InfinityRing.has_coerce_map_from(RDF)
            True
            sage: InfinityRing.has_coerce_map_from(CC)
            False
        """
        return sage.rings.real_mpfr.RealField(sage.rings.real_mpfr.mpfr_prec_min()).has_coerce_map_from(R)

class FiniteNumber(RingElement):
    def __init__(self, parent, x):
        """
        Initialize ``self``.

        TESTS::

            sage: sage.rings.infinity.FiniteNumber(InfinityRing, 1)
            A positive finite number
            sage: sage.rings.infinity.FiniteNumber(InfinityRing, -1)
            A negative finite number
            sage: sage.rings.infinity.FiniteNumber(InfinityRing, 0)
            Zero
        """
        RingElement.__init__(self, parent)
        self.value = x

    def __cmp__(self, other):
        """
        Compare ``self`` and ``other``.

        EXAMPLES::

            sage: P = InfinityRing
            sage: -oo < P(-5) < P(0) < P(1.5) < oo
            True
            sage: P(1) < P(100)
            False
            sage: P(-1) == P(-100)
            True
        """
        if isinstance(other, PlusInfinity):
            return -1
        if isinstance(other, MinusInfinity):
            return 1
        return cmp(self.value, other.value)

    def _add_(self, other):
        """
        EXAMPLES::

            sage: P = InfinityRing
            sage: 4 + oo # indirect doctest
            +Infinity
            sage: P(4) + P(2)
            A positive finite number
            sage: P(-1) + P(1)
            Traceback (most recent call last):
            ...
            SignError: cannot add positive finite value to negative finite value
        """
        if isinstance(other, InfinityElement):
            return other
        if self.value * other.value < 0:
            raise SignError, "cannot add positive finite value to negative finite value"
        return FiniteNumber(self.parent(), self.value)

    def _mul_(self, other):
        """
        EXAMPLES::

            sage: P = InfinityRing
            sage: 0 * oo # indirect doctest
            Traceback (most recent call last):
            ...
            SignError: cannot multiply infinity by zero
            sage: -1 * oo
            -Infinity
            sage: -2 * oo
            -Infinity
            sage: 3 * oo
            +Infinity
            sage: -oo * oo
            -Infinity
            sage: P(0) * 3
            0
            sage: P(-3) * P(2/3)
            A negative finite number
        """
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
        """
        EXAMPLES::

            sage: P = InfinityRing
            sage: 1 / oo # indirect doctest
            0
            sage: oo / 4
            +Infinity
            sage: oo / -4
            -Infinity
            sage: P(1) / P(-4)
            A negative finite number
        """
        return self * ~other

    def _sub_(self, other):
        """
        EXAMPLES::

            sage: P = InfinityRing
            sage: 4 - oo # indirect doctest
            -Infinity
            sage: 5 - -oo
            +Infinity
            sage: P(44) - P(4)
            Traceback (most recent call last):
            ...
            SignError: cannot add positive finite value to negative finite value
            sage: P(44) - P(-1)
            A positive finite number
        """
        return self._add_(-other)

    def __invert__(self):
        """
        EXAMPLES::

            sage: P = InfinityRing
            sage: ~P(2)
            A positive finite number
            sage: ~P(-7)
            A negative finite number
            sage: ~P(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot divide by zero
        """
        if self.value == 0:
            raise ZeroDivisionError, "Cannot divide by zero"
        return self

    def _neg_(self):
        """
        EXAMPLES::

            sage: a = InfinityRing(5); a
            A positive finite number
            sage: -a # indirect doctest
            A negative finite number
            sage: -(-a) == a
            True
            sage: -InfinityRing(0)
            Zero
        """
        return FiniteNumber(self.parent(), -self.value)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: InfinityRing(-2)._repr_()
            'A negative finite number'
            sage: InfinityRing(7)._repr_()
            'A positive finite number'
            sage: InfinityRing(0)._repr_()
            'Zero'
        """
        if self.value < 0:
            return "A negative finite number"
        if self.value > 0:
            return "A positive finite number"
        return "Zero"

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        TESTS::

            sage: a = InfinityRing(pi); a
            A positive finite number
            sage: a._latex_()
            'A positive finite number'
            sage: [latex(InfinityRing(a)) for a in [-2..2]]
            [A negative finite number, A negative finite number, Zero, A positive finite number, A positive finite number]
        """
        return self._repr_()

    def __abs__(self):
        """
        EXAMPLES::

            sage: abs(InfinityRing(-3))
            A positive finite number
            sage: abs(InfinityRing(3))
            A positive finite number
            sage: abs(InfinityRing(0))
            Zero
        """
        if self.value == 0:
            return FiniteNumber(self.parent(), 0)
        return FiniteNumber(self.parent(), 1)

    def sqrt(self):
        """
        EXAMPLES::

            sage: InfinityRing(7).sqrt()
            A positive finite number
            sage: InfinityRing(0).sqrt()
            Zero
            sage: InfinityRing(-.001).sqrt()
            Traceback (most recent call last):
            ...
            SignError: cannot take square root of a negative number
        """
        if self.value < 0:
            raise SignError, "cannot take square root of a negative number"
        return self

class MinusInfinity(_uniq, AnInfinity, MinusInfinityElement):

    _sign = -1
    _sign_char = '-'

    def __init__(self):
        """
        Initialize ``self``.

        TESTS::

            sage: sage.rings.infinity.MinusInfinity() is sage.rings.infinity.MinusInfinity() is -oo
            True
        """
        InfinityElement.__init__(self, InfinityRing)

    def _neg_(self):
        """
        EXAMPLES::

            sage: -(-oo) # indirect doctest
            +Infinity
        """
        return self.parent().gen(0)

    def sqrt(self):
        """
        EXAMPLES::

            sage: (-oo).sqrt()
            Traceback (most recent call last):
            ...
            SignError: cannot take square root of negative infinity
        """
        raise SignError, "cannot take square root of negative infinity"

    def _sympy_(self):
        """
        Converts ``-oo`` to sympy ``-oo``.

        Then you don't have to worry which ``oo`` you use, like in these
        examples:

        EXAMPLES::

            sage: import sympy
            sage: bool(-oo == -sympy.oo)
            True
            sage: bool(SR(-oo) == -sympy.oo)
            True
            sage: bool((-oo)._sympy_() == -sympy.oo)
            True

        """
        import sympy
        return -sympy.oo


class PlusInfinity(_uniq, AnInfinity, PlusInfinityElement):

    _sign = 1
    _sign_char = '+'

    def __init__(self):
        """
        Initialize ``self``.

        TESTS::

            sage: sage.rings.infinity.PlusInfinity() is sage.rings.infinity.PlusInfinity() is oo
            True
        """
        InfinityElement.__init__(self, InfinityRing)

    def _neg_(self):
        """
        TESTS::

            sage: -oo # indirect doctest
            -Infinity
        """
        return self.parent().gen(1)

    def sqrt(self):
        """
        The square root of ``self``.

        The square root of infinity is infinity.

        EXAMPLES::

            sage: oo.sqrt()
            +Infinity
        """
        return self

    def _sympy_(self):
        """
        Converts ``oo`` to sympy ``oo``.

        Then you don't have to worry which ``oo`` you use, like in these
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

