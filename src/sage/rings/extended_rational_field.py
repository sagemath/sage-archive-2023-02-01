#*****************************************************************************
#
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>
#                          William Stien <wstein@gmail.com>
#                          David Roe <roed314@gmail.com>
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

import sage.rings.rational_field
import sage.rings.rational
import sage.rings.integer
import sage.rings.infinity
import sage.structure.element
import sage.rings.extended_integer_ring
import field
from sage.structure.parent_gens import ParentWithGens
from sage.categories.morphism import Morphism


Rational = sage.rings.rational.Rational
RationalField = sage.rings.rational_field.RationalField
QQ = sage.rings.rational_field.QQ
Integer = sage.rings.integer.Integer
InfinityElement = sage.structure.element.InfinityElement
PlusInfinityElement = sage.structure.element.PlusInfinityElement
MinusInfinityElement = sage.structure.element.MinusInfinityElement
InfinityRing = sage.rings.infinity.InfinityRing
ExtendedIntegerRing = sage.rings.extended_integer_ring.ExtendedIntegerRing
IntegerPlusInfinity = sage.rings.extended_integer_ring.IntegerPlusInfinity
IntegerMinusInfinity = sage.rings.extended_integer_ring.IntegerMinusInfinity
SignError = sage.rings.infinity.SignError

import sage.rings.number_field.number_field_base as number_field_base

class Q_to_ExtendedQ(Morphism):
    def __init__(self, ExtQQ=None):
        """
        EXAMPLES:
            sage: from sage.rings.extended_rational_field import Q_to_ExtendedQ
            sage: f = Q_to_ExtendedQ()
            sage: loads(dumps(f))
            Natural morphism:
              From: Rational Field
              To:   Extended Rational Field
        """
        import sage.categories.homset
        if ExtQQ is None:
            ExtQQ = ExtendedRationalField
        Morphism.__init__(self, sage.categories.homset.Hom(QQ, ExtQQ))
        self._repr_type_str = "Natural"

    def _call_(self, x):
        """
        Returns the image of x under self.

        EXAMPLES:
            sage: from sage.rings.extended_rational_field import Q_to_ExtendedQ
            sage: f = Q_to_ExtendedQ()
            sage: f(QQ(2)) #indirect doctest
            2
            sage: type(_)
            <class 'sage.rings.extended_rational_field.ExtendedRational'>
            sage: f(2)  #indirect doctest
            2
        """
        return ExtendedRational(x)



_obj = {}
class _uniq0(object):
    def __new__(cls):
        if _obj.has_key(0):
            return _obj[0]
        O = number_field_base.NumberField.__new__(cls)
        _obj[0] = O
        return O

class _uniq1(object):
    def __new__(cls):
        if _obj.has_key(1):
            return _obj[1]
        O = PlusInfinityElement.__new__(cls)
        _obj[1] = O
        return O

class _uniq2(object):
    def __new__(cls):
        if _obj.has_key(2):
            return _obj[2]
        O = MinusInfinityElement.__new__(cls)
        _obj[2] = O
        return O

class ExtendedRationalField_class(_uniq0, RationalField):
    def __init__(self):
        """
        TESTS:
            sage: E = ExtendedRationalField
            sage: E == loads(dumps(E))
            True
        """
        ParentWithGens.__init__(self, self)
        self._assign_names(('x'),normalize=False) # ???
        self._populate_coercion_lists_(coerce_list=[Q_to_ExtendedQ(self)])

    def _repr_(self):
        """
        EXAMPLES:
            sage: ExtendedRationalField._repr_()
            'Extended Rational Field'
        """
        return "Extended Rational Field"

    def _latex_(self):
        """
        EXAMPLES:
            sage: latex(ExtendedRationalField) #indirect doctest
            \mathbf{Q}\cup\{\pm\infty\}
        """
        return "\\mathbf{Q}\\cup\\{\\pm\\infty\\}"

    def _element_constructor_(self, x, base = 0):
        """
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(oo)
            +Infinity
            sage: type(_)
            <class 'sage.rings.extended_rational_field.RationalPlusInfinity'>
            sage: E(-oo)
            -Infinity
            sage: E(0)
            0
            sage: E(2)
            2

        """
        if isinstance(x, sage.rings.infinity.MinusInfinity):
            return self.gen(2)
        if isinstance(x, sage.structure.element.InfinityElement):
            return self.gen(1)
        if isinstance(x, sage.rings.infinity.FiniteNumber):
            if x == 0:
                return ExtendedRational(0)
            raise TypeError, "cannot coerce unknown finite number into the extended rationals"
        return ExtendedRational(x, base)

    def _coerce_map_from_(self, S):
        """
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E.coerce_map_from(ZZ) #indirect doctest
            Composite map:
              From: Integer Ring
              To:   Extended Rational Field
              Defn:   Natural morphism:
                      From: Integer Ring
                      To:   Rational Field
                    then
                      Natural morphism:
                      From: Rational Field
                      To:   Extended Rational Field
            sage: E.coerce_map_from(QQ) #indirect doctest
            Natural morphism:
              From: Rational Field
              To:   Extended Rational Field
            sage: E.coerce_map_from(ExtendedIntegerRing)
            Conversion map:
              From: Extended Integer Ring
              To:   Extended Rational Field

            sage: E.coerce(int(2))
            2
            sage: E.coerce(1/2)
            1/2
            sage: E.coerce(oo)
            +Infinity
            sage: R.<x> = QQ[]
            sage: E.coerce(R(1))
            Traceback (most recent call last):
            ...
            TypeError: no cannonical coercion from Univariate Polynomial Ring in x over Rational Field to Extended Rational Field

        TESTS:
            sage: ExtendedRationalField(2)*ExtendedIntegerRing(2)
            4
            sage: type(_)
            <class 'sage.rings.extended_rational_field.ExtendedRational'>

        """
        if S is QQ:
            return Q_to_ExtendedQ()
        elif S == ExtendedIntegerRing or S == InfinityRing:
            return self._generic_convert_map(S)

    def _is_valid_homomorphism(self, codomain, im_gens):
        """
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E._is_valid_homomorphism(E, (0,0,0))
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __iter__(self):
        """
        Returns a generator for the elements of the extended rational
        field.

        EXAMPLES:
            sage: it = iter(ExtendedRationalField)
            sage: [it.next() for _ in range(10)]
            [0, 1, -1, +Infinity, -Infinity, 2, 1/2, -2, -1/2, 3]
        """
        yield self(0)
        yield self(1)
        yield self(-1)
        yield self.gen(1)
        yield self.gen(2)
        from integer_ring import IntegerRing
        for n in IntegerRing():
            m = abs(n)
            for d in abs(n).coprime_integers(m):
                yield n/d
                yield d/n

    def complex_embedding(self, prec=53):
        """
        EXAMPLES:
            sage: ExtendedRationalField.complex_embedding()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def gens(self):
        """
        Returns the generators of the extended rational field.

        EXAMPLES:
            sage: ExtendedRationalField.gens()
            (1, +Infinity, -Infinity)
        """
        return (self(1), self.gen(1), self.gen(2), )

    def gen(self, n=0):
        r"""
        Returns the $n^{th} generator of the extended rational field.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E.gen()
            1
            sage: E.gen(0)
            1
            sage: E.gen(1)
            +Infinity
            sage: E.gen(2)
            -Infinity
            sage: E.gen(3)
            Traceback (most recent call last):
            ...
            IndexError: n must be 0, 1, or 2
        """
        if n == 0:
            return self(1)
        elif n == 1:
            try:
                return self.gen1
            except AttributeError:
                self.gen1 = RationalPlusInfinity()
                return self.gen1
        elif n == 2:
            try:
                return self.gen2
            except AttributeError:
                self.gen2 = RationalMinusInfinity()
                return self.gen2
        else:
            raise IndexError, "n must be 0, 1, or 2"

    def is_prime_field(self):
        """
        Returns whether or not the extended rational field is a
        prime field (which it is not).

        EXAMPLES:
            sage: ExtendedRationalField.is_prime_field()
            False
        """
        return False

    def ngens(self):
        """
        Returns the number of generators of the extended rational field.

        EXAMPLES:
            sage: ExtendedRationalField.ngens()
            3
        """
        return 3


ExtendedRationalField = ExtendedRationalField_class()

class ExtendedRational(Rational):
    def __init__(self, x = None, base = 0):
        """
        Constructor for elements of the extended rational field.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: a = E(2)
            sage: a == loads(dumps(a))
            True
        """
        Rational.__init__(self, x, base)
        self._set_parent(ExtendedRationalField)

    def _rational_(self):
        """
        Returns a rational version of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(2)._rational_()
            2
            sage: type(_)
            <type 'sage.rings.rational.Rational'>
            sage: type(QQ(E(1/2)))
            <type 'sage.rings.rational.Rational'>
        """
        return Rational.copy(self)

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: cmp(E(2),E(3))     #indirect doctest
            -1
            sage: cmp(E(2),E(-2))    #indirect doctest
            1
            sage: cmp(E(2),E(-oo))   #indirect doctest
            1
            sage: cmp(E(2), E(oo))   #indirect doctest
            -1
            sage: cmp(E(2), E(2))    #indirect doctest
            0
            sage: cmp(E(oo), E(oo))  #indirect doctest
            0
            sage: cmp(E(oo), E(-oo)) #indirect doctest
            1
        """
        if isinstance(other, InfinityElement):
            return -other.__cmp__(self)
        return cmp(Rational(self), Rational(other))

    def copy(self):
        """
        Returns a copy of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: a = E(2)
            sage: b = a.copy()
            sage: a == b
            True
            sage: a is b
            False
        """
        return self.parent()(Rational.copy(self))

    def lcm(self, other):
        """
        Returns the least common multiple of self and other.  If other is plus or
        minus infinity, then the lcm is +Infinity.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(2).lcm(3)
            6
            sage: E(2).lcm(oo)
            +Infinity
            sage: E(2).lcm(-oo)
            +Infinity
        """
        if isinstance(self, InfinityElement) or isinstance(other, InfinityElement):
            return self.parent().gen(1)
        else:
            return self.parent()(Rational.lcm(self, QQ(other)))

    def sqrt(self):
        """
        Returns the square root of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(2).sqrt()
            sqrt(2)
            sage: E(-2).sqrt()
            sqrt(2)*I
            sage: E(4).sqrt()
            2
            sage: type(_)
            <class 'sage.rings.extended_rational_field.ExtendedRational'>
            sage: sqrt(E(2))
            sqrt(2)

        """
        value = Rational.sqrt(self)
        try:
            return self.parent()(value)
        except TypeError:
            return value

    def nth_root(self, n):
        """
        Computes the nth root of self, or raises a \exception{ValueError}
        if self is not a perfect nth power.

        INPUT:
            n -- integer (must fit in C int type)

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(8).nth_root(3)
            2
            sage: E(7).nth_root(3)
            Traceback (most recent call last):
            ...
            ValueError: not a perfect nth power
        """
        return self.parent()(Rational.nth_root(self, n))

    def _add_(self, right):
        """
        Returns the sum of self and right.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(2) + E(4) #indirect doctest
            6
            sage: E(2) + 4
            6
            sage: E(2) + E(oo)
            +Infinity
            sage: E(2) + E(-oo)
            -Infinity
        """
        if isinstance(right, InfinityElement):
            return right
        return self.parent()(Rational(self) + Rational(right))

    def _sub_(self, right):
        """
        Returns the difference between self and right.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(2) - E(4) #indirect doctest
            -2
            sage: E(2) - E(oo)
            -Infinity
            sage: E(2) - E(-oo)
            +Infinity
        """
        if isinstance(right, InfinityElement):
            return -right
        return self.parent()(Rational(self) - Rational(right))

    def _neg_(self):
        """
        Returns the negation of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: -E(2) #indirect doctest
            -2
            sage: -E(oo)
            -Infinity
            sage: -E(-oo)
            +Infinity
        """
        return self.parent()(-Rational(self))

    def _mul_(self, right):
        """
        Returns the product of self and right.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(2)*E(4) #indirect doctest
            8
            sage: E(2)*E(oo)
            +Infinity
            sage: E(2)*E(-oo)
            -Infinity
            sage: E(-2)*E(oo)
            -Infinity
            sage: E(-2)*E(-oo)
            +Infinity
            sage: E(4)*2
            8
            sage: E(0)*E(oo)
            Traceback (most recent call last):
            ...
            SignError: cannot multiply infinity by zero
        """
        if isinstance(right, InfinityElement):
            return right._mul_(self)
        return self.parent()(Rational(self) * Rational(right))

    def _div_(self, right):
        """
        Returns self / right.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(2)/4 #indirect doctest
            1/2
            sage: E(2)/E(oo)
            0
            sage: E(2)/E(-oo)
            0
            sage: E(-2)/E(oo)
            0
        """
        if isinstance(right, InfinityElement):
            return self.parent()(0)
        return self.parent()(Rational(self) / Rational(right))

    def __invert__(self):
        """
        Returns the inverse of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: ~E(2) #indirect doctest
            1/2
            sage: ~E(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
        """
        return self.parent()(~Rational(self))

    def __pow__(self, n):
        """
        Returns self^n.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(2)^2 #indirect doctest
            4
            sage: E(2)^1
            2
            sage: E(2)^(1/3)
            2^(1/3)
            sage: E(1)^E(oo)
            1
            sage: E(2)^E(oo)
            +Infinity
            sage: E(0)^E(oo)
            Traceback (most recent call last):
            ...
            SignError: 0^infinity is not defined
            sage: E(-1/2)^E(oo)
            0
            sage: E(-2)^E(oo)
            Traceback (most recent call last):
            ...
            SignError: negative^infinity is not defined

            sage: E(2)^E(-oo)
            0
            sage: E(1)^E(-oo)
            1
            sage: E(1/2)^E(-oo)
            +Infinity
            sage: E(-1)^E(-oo)
            Traceback (most recent call last):
            ...
            SignError: x^(-infinity) not defined for -1 <= x <= 0
            sage: E(-2)^E(-oo)
            0
        """
        if isinstance(n, InfinityElement):
            if isinstance(n, PlusInfinityElement):
                if self > 1:
                    return self.parent().gen(1)
                elif self == 1:
                    return self
                elif self > 0:
                    return self.parent()(0)
                elif self == 0:
                    raise SignError, "0^infinity is not defined"
                elif self > -1:
                    return self.parent()(0)
                else:
                    raise SignError, "negative^infinity is not defined"
            elif isinstance(n, MinusInfinityElement):
                if self > 1:
                    return self.parent()(0)
                elif self == 1:
                    return self
                elif self > 0:
                    return self.parent().gen(1)
                elif self >= -1:
                    raise SignError, "x^(-infinity) not defined for -1 <= x <= 0"
                else:
                    return self.parent()(0)
            else:
                raise TypeError, "cannot raise n to an unsigned infinite power."
        value = Rational(self)**n
        try:
            return self.parent()(value)
        except TypeError:
            return value

    def __abs__(self):
        """
        Returns the absolute value of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: abs(E(-2)) #indirect doctest
            2
            sage: abs(E(2))  #indirect doctest
            2
        """
        return self.parent()(Rational(self).__abs__())

    def numerator(self):
        """
        Returns the numerator of self as an extended integer.  If you want an actual integer,
        use numer instead.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(2).numerator()
            2
            sage: E(3/4).numerator()
            3
        """
        return ExtendedIntegerRing(Rational(self).numerator())

    def denominator(self):
        """
        Returns the denominator of self as an extended integer.  If you want an actual integer,
        use denom instead.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(2).denominator()
            1
            sage: E(3/4).denominator()
            4
        """
        return ExtendedIntegerRing(Rational(self).denominator())

    def floor(self):
        """
        Returns the floor of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(4/3).floor()
            1
            sage: floor(E(4/3))
            1
            sage: E(-4/3).floor()
            -2
        """
        return ExtendedIntegerRing(Rational(self).floor())

    def ceil(self):
        """
        Returns the ceiling of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(4/3).ceil()
            2
            sage: ceil(E(4/3))
            2
            sage: E(-4/3).ceil()
            -1
        """
        return ExtendedIntegerRing(Rational(self).ceil())

    def __lshift__(self, n):
        r"""
        Returns $self*2^n$.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(1/2) << 3 #indirect doctest
            4
            sage: E(-2) << 2
            -8
        """
        return self.parent()(Rational(self).__lshift__(n))

    def __rshift__(self, n):
        r"""
        Returns $self*2^{-n}$.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(2) >> 1 #indirect doctest
            1
            sage: E(-8) >> 2
            -2
        """
        return self.parent()(Rational(self).__rshift__(n))


class RationalPlusInfinity(_uniq1, PlusInfinityElement):
    def __init__(self):
        """
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: a = E(oo)
            sage: a == loads(dumps(a))
            True
        """
        PlusInfinityElement.__init__(self, ExtendedRationalField)

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: cmp(E(oo), 3) #indirect doctest
            1
            sage: cmp(E(oo), E(oo)) #indirect doctest
            0
        """
        if isinstance(other, RationalPlusInfinity):
            return 0
        return 1

    def __repr__(self):
        """
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: repr(E(oo)) #indirect doctest
            '+Infinity'
        """
        return "+Infinity"

    def _latex_(self):
        """
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: latex(E(oo)) #indirect doctest
            +\infty
        """
        return "+\\infty"

    def _add_(self, other):
        """
        Returns self + other.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(oo) + 2 #indirect doctest
            +Infinity
            sage: E(oo) + E(-oo)
            Traceback (most recent call last):
            ...
            SignError: cannot add infinity to minus infinity
        """

        if isinstance(other, RationalMinusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _mul_(self, other):
        """
        Returns self*other.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(oo)*2   #indirect doctest
            +Infinity
            sage: E(oo)*E(2)
            +Infinity
            sage: E(oo)*E(-2)
            -Infinity
            sage: E(oo)*E(oo)
            +Infinity
            sage: E(oo)*E(-oo)
            -Infinity
            sage: E(oo)*E(0)
            Traceback (most recent call last):
            ...
            SignError: cannot multiply infinity by zero
        """

        if other < 0:
            return -self
        if other > 0:
            return self
        raise SignError, "cannot multiply infinity by zero"

    def _sub_(self, other):
        """
        Returns self-other.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(oo) - E(2) #indirect doctest
            +Infinity
            sage: E(oo) - E(-oo)
            +Infinity
            sage: E(oo) - E(oo)
            Traceback (most recent call last):
            ...
            SignError: cannot subtract infinity from infinity
        """

        if isinstance(other, RationalPlusInfinity):
            raise SignError, "cannot subtract infinity from infinity"
        return self

    def _div_(self, other):
        """
        Returns self/other.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(oo)/E(2) #indirect doctest
            +Infinity
        """

        return self * other.__invert__()

    def _neg_(self):
        """
        Returns the negation of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: -E(oo) #indirect doctest
            -Infinity
        """
        return self.parent().gen(2)

    def __invert__(self):
        """
        Returns 1 / self.
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: ~E(oo) #indirect doctest
            0
        """

        return ExtendedRational(0)

    def __abs__(self):
        """
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: abs(E(oo)) #indirect doctest
            +Infinity
        """
        return self

    def __pow__(self, right):
        """
        Returns self^right.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(oo)^2 #indirect doctest
            +Infinity
            sage: E(oo)^(-2)
            0
            sage: E(oo)^0
            Traceback (most recent call last):
            ...
            SignError: Cannot raise infinity to the zeroth power
            sage: E(oo)^E(oo)
            +Infinity
            sage: R.<x> = QQ[]
            sage: E(oo)^x
            Traceback (most recent call last):
            ...
            TypeError: cannot exponentiate
        """

        if not isinstance(right, (int, long, Integer, Rational, PlusInfinityElement, MinusInfinityElement)):
            raise TypeError, "cannot exponentiate"
        if right < 0:
            return ExtendedRational(0)
        elif right > 0:
            return self
        else:
            raise SignError, "Cannot raise infinity to the zeroth power"

    def lcm(self, x):
        """
        Return the least common multiple of oo and x, which
        is by definition oo unless x is 0.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(oo).lcm(0)
            0
            sage: E(oo).lcm(oo)
            +Infinity
            sage: E(oo).lcm(10)
            +Infinity
        """
        if x == 0:
            return x
        else:
            return self

    def sqrt(self):
        """
        Returns the square root of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(oo).sqrt()
            +Infinity
        """
        return self

    def nth_root(self, n):
        r"""
        Returns the $n^{th}$ root of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(oo).nth_root(4)
            +Infinity
            sage: E(oo).nth_root(2)
            +Infinity
        """
        return self

    def floor(self):
        """
        Returns the floor of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(oo).floor()
            +Infinity
        """
        return IntegerPlusInfinity()

    def ceil(self):
        """
        Returns the ceiling of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(oo).ceil()
            +Infinity
        """
        return IntegerPlusInfinity()

    def numerator(self):
        """
        Returns the numerator of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(oo).numerator()
            +Infinity
        """
        return IntegerPlusInfinity()

    def denominator(self):
        """
        Returns the denominator of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(oo).denominator()
            1
        """
        return ExtendedIntegerRing(1)

class RationalMinusInfinity(_uniq2, MinusInfinityElement):
    def __init__(self):
        """
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: a = E(-oo)
            sage: a == loads(dumps(a))
            True
        """
        InfinityElement.__init__(self, ExtendedRationalField)

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: cmp(E(-oo), 2) #indirect doctest
            -1
            sage: cmp(E(-oo), E(-oo)) #indirect doctest
            0
        """

        if isinstance(other, RationalMinusInfinity):
            return 0
        return -1

    def __repr__(self):
        """
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(-oo).__repr__()
            '-Infinity'
        """
        return "-Infinity"

    def _latex_(self):
        """
        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: latex(E(-oo)) #indirect doctest
            -\infty
        """

        return "-\\infty"

    def _add_(self, other):
        """
        Returns self + other.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(-oo) + 2 #indirect doctest
            -Infinity
            sage: E(-oo) + E(-oo)
            -Infinity
            sage: E(-oo) + E(oo)
            Traceback (most recent call last):
            ...
            SignError: cannot add infinity to minus infinity
        """

        if isinstance(other, RationalPlusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _mul_(self, other):
        """
        Returns self*other.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(-oo)*2 #indirect doctest
            -Infinity
            sage: E(-oo)*-2
            +Infinity
            sage: E(-oo)*E(-oo)
            +Infinity
            sage: E(-oo)*E(oo)
            -Infinity
            sage: E(-oo)*0
            Traceback (most recent call last):
            ...
            SignError: cannot multiply minus infinity by zero
        """

        if other < 0:
            return -self
        if other > 0:
            return self
        raise SignError, "cannot multiply minus infinity by zero"

    def _sub_(self, other):
        """
        Returns self - other.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(-oo) - 2 #indirect doctest
            -Infinity
            sage: E(-oo) - E(oo)
            -Infinity
            sage: E(-oo) - E(-oo)
            Traceback (most recent call last):
            ...
            SignError: cannot subtract minus infinity from minus infinity
        """
        if isinstance(other, RationalMinusInfinity):
            raise SignError, "cannot subtract minus infinity from minus infinity"
        return self

    def _div_(self, other):
        """
        Returns self / other.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(-oo)/2 #indirect doctest
            -Infinity
        """
        return self * other.__invert__()

    def _neg_(self):
        """
        Returns the negation of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: -E(-oo) #indirect doctest
            +Infinity
        """
        return self.parent().gen(1)

    def __invert__(self):
        """
        Returns 1 / self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: ~E(-oo) #indirect doctest
            0
        """
        return ExtendedRational(0)

    def __abs__(self):
        """
        Returns the absolute value of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: abs(E(-oo)) #indirect doctest
            +Infinity
        """
        return -self

    def __pow__(self, right):
        """
        Returns self^right.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(-oo)^2 #indirect doctest
            +Infinity
            sage: E(-oo)^(1/3)
            -Infinity
            sage: E(-oo)^(2/3)
            +Infinity
            sage: E(-oo)^(1/4)
            Traceback (most recent call last):
            ...
            SignError: Cannot take an even root of minus infinity
            sage: E(-oo)^(-2)
            0
            sage: E(-oo)^4
            +Infinity
            sage: E(-oo)^0
            Traceback (most recent call last):
            ...
            SignError: Cannot raise minus infinity to the zeroth power
            sage: R.<x> = QQ[]
            sage: E(-oo)^x
            Traceback (most recent call last):
            ...
            TypeError: cannot exponentiate
        """
        if isinstance(right, Rational):
            if right.denominator() % 2 == 0:
                raise SignError, "Cannot take an even root of minus infinity"
        elif not isinstance(right, (int, long, Integer, PlusInfinityElement, MinusInfinityElement)):
            raise TypeError, "cannot exponentiate"
        if right < 0:
            return ExtendedRational(0)
        elif right > 0:
            right = QQ(right)
            if right.numerator() % 2 == 0:
                return -self
            else:
                return self
        else:
            raise SignError, "Cannot raise minus infinity to the zeroth power"

    def lcm(self, x):
        """
        Return the least common multiple of -oo and x, which
        is by definition oo unless x is 0.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: moo = E(-oo)
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
        """
        Returns the square root of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(-oo).sqrt()
            Traceback (most recent call last):
            ...
            SignError: cannot take the square root of minus infinity
        """
        raise SignError, "cannot take the square root of minus infinity"

    def nth_root(self, n):
        """
        Returns the nth root of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(-oo).nth_root(3)
            -Infinity
            sage: E(-oo).nth_root(4)
            Traceback (most recent call last):
            ...
            SignError: cannot take an even root of minus infinity
        """
        if n % 2 == 0:
            raise SignError, "cannot take an even root of minus infinity"
        return self

    def floor(self):
        """
        Returns the floor of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(-oo).floor()
            -Infinity
        """
        return IntegerMinusInfinity()

    def ceil(self):
        """
        Returns the ceiling of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(-oo).ceil()
            -Infinity
        """

        return IntegerMinusInfinity()

    def numerator(self):
        """
        Returns the numerator of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(-oo).numerator()
            -Infinity
        """
        return IntegerMinusInfinity()

    def denominator(self):
        """
        Returns the denominator of self.

        EXAMPLES:
            sage: E = ExtendedRationalField
            sage: E(-oo).denominator()
            1
        """
        return ExtendedIntegerRing(1)
