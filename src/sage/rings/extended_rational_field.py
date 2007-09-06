import sage.rings.rational_field
import sage.rings.rational
import sage.rings.integer
import sage.rings.infinity
import sage.structure.element
import sage.rings.extended_integer_ring
import field
from sage.structure.parent_gens import ParentWithGens

Rational = sage.rings.rational.Rational
RationalField = sage.rings.rational_field.RationalField
Integer = sage.rings.integer.Integer
InfinityElement = sage.structure.element.InfinityElement
PlusInfinityElement = sage.structure.element.PlusInfinityElement
MinusInfinityElement = sage.structure.element.MinusInfinityElement
ExtendedIntegerRing = sage.rings.extended_integer_ring.ExtendedIntegerRing
IntegerPlusInfinity = sage.rings.extended_integer_ring.IntegerPlusInfinity
IntegerMinusInfinity = sage.rings.extended_integer_ring.IntegerMinusInfinity
SignError = sage.rings.infinity.SignError

import sage.rings.number_field.number_field_base as number_field_base

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
        ParentWithGens.__init__(self, self)
        self._assign_names(('x'),normalize=False)

    def _repr_(self):
        return "Extended Rational Field"

    def _latex_(self):
        return "\\mathbf{Q}\\cup\\{\\pm\\infty\\}"

    def __call__(self, x, base = 0):
        if isinstance(x, sage.rings.infinity.MinusInfinity):
            return self.gen(2)
        if isinstance(x, sage.structure.element.InfinityElement):
            return self.gen(1)
        if isinstance(x, sage.rings.infinity.FiniteNumber):
            if x == 0:
                return ExtendedRational(0)
            raise TypeError, "cannot coerce unknown finite number into the extended rationals"
        return ExtendedRational(x, base)

    def _coerce_impl(self, x):
        if isinstance(x, (int, long, sage.rings.integer.Integer, Rational)):
            return self(x)
        if isinstance(x, (IntegerPlusInfinity, IntegerMinusInfinity)):
            return self(x)
        raise TypeError, "no implicit coercion of element to the rational numbers"

    def _is_valid_homomorphism(self, codomain, im_gens):
        raise NotImplementedError

    def __iter__(self):
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
        raise NotImplementedError

    def gens(self):
        return (self(1), self.gen(1), self.gen(2), )

    def gen(self, n=0):
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
            raise IndexError, "n must be 0, 1 or 2"

    def is_prime_field(self):
        return False

    def ngens(self):
        return 3

    def numberfield(self, poly_var, nf_var):
        raise NotImplementedError

ExtendedRationalField = ExtendedRationalField_class()

class ExtendedRational(Rational):
    def __init__(self, x = None, base = 0):
        Rational.__init__(self, x, base)
        self._set_parent(ExtendedRationalField)

    def __cmp__(self, other):
        if isinstance(other, InfinityElement):
            return -other.__cmp__(self)
        return cmp(Rational(self), Rational(other))

    def copy(self):
        return self.parent()(Rational.copy(self))

    def lcm(self, other):
        if isinstance(other, InfinityElement):
            return self.parent().gen(1)
        else:
            return self.parent()(Rational.lcm(self, other))

    def square_root(self):
        return self.parent()(Rational.square_root(self))

    def nth_root(self):
        return self.parent()(Rational.nth_root(self))

    def _add_(self, right):
        if isinstance(right, InfinityElement):
            return right
        return self.parent()(Rational(self) + Rational(right))

    def _sub_(self, right):
        if isinstance(right, InfinityElement):
            return -right
        return self.parent()(Rational(self) - Rational(right))

    def _neg_(self):
        return self.parent()(-Rational(self))

    def _mul_(self, right):
        if isinstance(right, InfinityElement):
            return right._mul_(self)
        return self.parent()(Rational(self) * Rational(right))

    def _div_(self, right):
        if isinstance(right, InfinityElement):
            return self.parent()(0)
        return self.parent()(Rational(self) / Rational(right))

    def __invert__(self):
        return self.parent()(~Rational(self))

    def __pow__(self, n):
        if isinstance(n, InfinityElement):
            if isinstance(n, PlusInfinityElement):
                if self > 1:
                    return self.parent().gen(1)
                elif self == 1:
                    return self
                elif self > 0:
                    return self.parent()(0)
                elif self == 0:
                    raise SignError, "0^infinity not defined"
                elif self > -1:
                    return self.parent()(0)
                else:
                    raise SignError, "negative^infinity not defined"
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
        return self.parent()(Rational(self)**n)

    def __abs__(self):
        return self.parent()(Rational(self).__abs__())

    def numerator(self):
        """
        Returns the numerator of self as an extended integer.  If you want an actual integer, use numer instead.
        """
        return ExtendedIntegerRing(Rational(self).numerator())

    def denominator(self):
        """
        Returns the denominator of self as an extended integer.  If you want an actual integer, use denom instead.
        """
        return ExtendedIntegerRing(Rational(self).denominator())

    def floor(self):
        return ExtendedIntegerRing(Rational(self).floor())

    def ceil(self):
        return ExtendedIntegerRing(Rational(self).ceil())

    def __lshift__(self, n):
        return self.parent()(Rational(self).__lshift__(n))

    def __rshift__(self, n):
        return self.parent()(Rational(self).__rshift__(n))


class RationalPlusInfinity(_uniq1, PlusInfinityElement):
    def __init__(self):
        PlusInfinityElement.__init__(self, ExtendedRationalField)

    def __cmp__(self, other):
        if isinstance(other, RationalPlusInfinity):
            return 0
        return 1

    def __repr__(self):
        return "+Infinity"

    def _latex_(self):
        return "+\\infty"

    def _add_(self, other):
        if isinstance(other, RationalMinusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _mul_(self, other):
        if other < 0:
            return -self
        if other > 0:
            return self
        raise TypeError, "cannot multiply infinity by zero"

    def _sub_(self, other):
        if isinstance(other, RationalPlusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _div_(self, other):
        return self * other.__invert__()

    def _neg_(self):
        return self.parent().gen(2)

    def __invert__(self):
        return ExtendedRational(0)

    def __abs__(self):
        return self

    def __pow__(self, right):
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

    def nth_root(self, n):
        return self

    def floor(self):
        return IntegerPlusInfinity()

    def ceil(self):
        return IntegerPlusInfinity()

    def numerator(self):
        return IntegerPlusInfinity()

    def denominator(self):
        return ExtendedIntegerRing(1)

class RationalMinusInfinity(_uniq2, MinusInfinityElement):
    def __init__(self):
        InfinityElement.__init__(self, ExtendedRationalField)

    def __cmp__(self, other):
        if isinstance(other, RationalMinusInfinity):
            return 0
        return -1

    def _repr_(self):
        return "-Infinity"

    def _latex_(self):
        return "-\\infty"

    def _add_(self, other):
        if isinstance(other, RationalPlusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _mul_(self, other):
        if other < 0:
            return -self
        if other > 0:
            return self
        raise SignError, "cannot multiply infinity by zero"

    def _sub_(self, other):
        if isinstance(other, RationalMinusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _div_(self, other):
        return self * other.__invert__()

    def _neg_(self):
        return self.parent().gen(1)

    def __invert__(self):
        return ExtendedRational(0)

    def __abs__(self):
        return self.parent().gen(1)

    def __pow__(self, right):
        if isinstance(right, Rational):
            if right.denominator() % 2 == 0:
                raise SignError, "Cannot take an even root of negative infinity"
        elif not isinstance(right, (int, long, Integer, PlusInfinityElement, MinusInfinityElement)):
            raise TypeError, "cannot exponentiate"
        if right < 0:
            return ExtendedRational(0)
        elif right > 0:
            return self
        else:
            raise SignError, "Cannot raise negative infinity to the zeroth power"

    def lcm(self, x):
        """
        Return the least common multiple of -oo and x, which
        is by definition oo unless x is 0.

        EXAMPLES:
            sage: moo = ExtendedRationalField.gen(2)
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

    def nth_root(self, n):
        if n % 2 == 0:
            raise SignError, "cannot take an even root of negative infinity"
        return self

    def floor(self):
        return IntegerMinusInfinity()

    def ceil(self):
        return IntegerMinusInfinity()

    def numerator(self):
        return IntegerMinusInfinity()

    def denominator(self):
        return ExtendedIntegerRing(1)

