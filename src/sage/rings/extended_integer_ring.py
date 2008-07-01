import sage.rings.integer_ring
import sage.rings.rational
import sage.rings.integer
import sage.rings.infinity
import sage.structure.element

ParentWithGens = sage.structure.parent_gens.ParentWithGens
Rational = sage.rings.rational.Rational
Integer = sage.rings.integer.Integer
IntegerWrapper = sage.rings.integer.IntegerWrapper
IntegerRing_class = sage.rings.integer_ring.IntegerRing_class
InfinityElement = sage.structure.element.InfinityElement
PlusInfinityElement = sage.structure.element.PlusInfinityElement
MinusInfinityElement = sage.structure.element.MinusInfinityElement
InfinityElement = sage.structure.element.InfinityElement
SignError = sage.rings.infinity.SignError
ZZ = sage.rings.integer_ring.ZZ

_obj = {}
class _uniq0(object):
    def __new__(cls):
        if _obj.has_key(0):
            return _obj[0]
        O = IntegerRing_class.__new__(cls)
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

class ExtendedIntegerRing_class(_uniq0, IntegerRing_class):
    def __init__(self):
        ParentWithGens.__init__(self, self)
        self._assign_names(('x'),normalize=False)

    def _repr_(self):
        return "Extended Integer Ring"

    def _latex_(self):
        return "\\mathbf{Z}\\cup\\{\\pm\\infty\\}"

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: cmp(ExtendedIntegerRing, ExtendedRationalField) #random due to architecture dependence
            1
            sage: cmp(ExtendedIntegerRing, ExtendedIntegerRing)
            0
        """
        return cmp(other.__class__, ExtendedIntegerRing_class)


    def __call__(self, x, base = 0):
        if isinstance(x, sage.rings.infinity.MinusInfinity):
            return self.gen(2)
        if isinstance(x, sage.structure.element.InfinityElement):
            return self.gen(1)
        if isinstance(x, sage.rings.infinity.FiniteNumber):
            if x == 0:
                return ExtendedInteger(0)
            raise TypeError, "cannot coerce unknown finite number into the extended integers"
        return ExtendedInteger(x, base)

    def _coerce_impl(self, x):
        if isinstance(x, (int, long, sage.rings.integer.Integer)):
            return self(x)
        raise TypeError, "no implicit coercion of element to the extended integer ring"

    def _is_valid_homomorphism(self, codomain, im_gens):
        raise NotImplementedError

    def __iter__(self):
        yield self(0)
        yield self(1)
        yield self(-1)
        yield self.gen(1)
        yield self.gen(2)
        n = self(2)
        while True:
            yield n
            yield -n
            n = n + 1

    def random_element(self, x=None, y=None):
        if y is None:
            if x is None:
                a = ZZ.random_element(-3, 4)
                if a == -3:
                    return self.gen(2)
                elif a == 3:
                    return self.gen(1)
                else:
                    return self(a)
            else:
                a = ZZ.random_element(x + 1)
                if a == x:
                    return self.gen(1)
                else:
                    return self(a)
        else:
            a = ZZ.random_element(x - 1, y + 1)
            if a == x - 1:
                return self.gen(2)
            elif a == y:
                return self.gen(1)
            else:
                return self(a)

    def fraction_field(self):
        from sage.rings.extended_rational_field import ExtendedRationalField
        return ExtendedRationalField

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
                self.gen1 = IntegerPlusInfinity()
                return self.gen1
        elif n == 2:
            try:
                return self.gen2
            except AttributeError:
                self.gen2 = IntegerMinusInfinity()
                return self.gen1
        else:
            raise IndexError, "n must be 0, 1 or 2"

    def ngens(self):
        return 3

    def numberfield(self, poly_var, nf_var):
        raise NotImplementedError

    def zeta(self, n=2):
        if n == 1:
            return self(1)
        elif n == 2:
            return self(-1)
        else:
            raise ValueError, "no nth root of unity in the extended integer ring"

ExtendedIntegerRing = ExtendedIntegerRing_class()

class ExtendedInteger(IntegerWrapper):
    def __init__(self, x = None, base = 0):
        IntegerWrapper.__init__(self, x, base)
        self._set_parent(ExtendedIntegerRing)

    def __cmp__(self, other):
        if isinstance(other, InfinityElement):
            return -other.__cmp__(self)
        return cmp(Integer(self), Integer(other))

    def __xor__(self, other):
        if isinstance(other, InfinityElement):
            return other.__xor__(self)
        return self.parent()(Integer.__xor__(self, other))

    def copy(self):
        return self.parent()(Integer.copy(self))

    def lcm(self, other):
        if isinstance(other, InfinityElement):
            return self.parent().gen(1)
        else:
            return self.parent()(Integer.lcm(self, other))

    def gcd(self, other):
        if isinstance(other, InfinityElement):
            return self
        else:
            return self.parent()(Integer.gcd(self, other))

    def square_root(self):
        return self.parent()(Integer.square_root(self))

    def nth_root(self, n, report_exact=0):
        x, exact = Integer.nth_root(self, n, report_exact)
        if report_exact:
            return self.parent()(x), exact
        else:
            return self.parent()(x)

    def _add_(self, right):
        if isinstance(right, InfinityElement):
            return right
        return self.parent()(Integer(self) + Integer(right))

    def _sub_(self, right):
        if isinstance(right, InfinityElement):
            return -right
        return self.parent()(Integer(self) - Integer(right))

    def _neg_(self):
        return self.parent()(-Integer(self))

    def _mul_(self, right):
        if isinstance(right, InfinityElement):
            return right._mul_(self)
        return self.parent()(Integer(self) * Integer(right))

    def _div_(self, right):
        from sage.rings.extended_rational_field import ExtendedRationalField
        if isinstance(right, InfinityElement):
            return ExtendedRationalField(0)
        return ExtendedRationalField(Integer(self) / Integer(right))

    def __floordiv__(self, right):
        """
        EXAMPLES:
            sage: R = ExtendedIntegerRing
            sage: R(3) // R(17)
            0
            sage: R(17) // R(3)
            5
            sage: R(25) // Infinity
            0
        """
        if isinstance(right, InfinityElement):
            return self.parent()(0)
        else:
            return self.parent()(Integer(self).__floordiv__(Integer(right)))

    def __invert__(self):
        from sage.rings.extended_rational_field import ExtendedRationalField
        return ExtendedRationalField(~Integer(self))

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
        return self.parent()(Integer(self)**n)

    def __abs__(self):
        return self.parent()(Integer(self).__abs__())

    def __mod__(self, right):
        if isinstance(right, InfinityElement):
            return self
        return self.parent()(Integer(self) % Integer(right))

    def quo_rem(self, other):
        if isinstance(other, InfinityElement):
            return self.parent()(0), self
        else:
            a, b = Integer(self).quo_rem(other)
            return self.parent()(a), self.parent()(b)

    def powermod(self, exp, mod):
        if isinstance(mod, InfinityElement):
            return self.__pow__(exp)
        if isinstance(exp, InfinityElement):
            raise TypeError, "result not well defined"
        return self.parent()(Integer(self).powermod(exp, mod))

    def coprime_integers(self, m):
        if isinstance(m, InfinityElement):
            raise NotImplementedError
        return [self.parent()(x) for x in Integer.coprime_integers(self, m)]

    def divides(self, n):
        if isinstance(n, InfinityElement):
            return True
        return Integer(self).divides(n)

    def numerator(self):
        return self

    def denominator(self):
        return self.parent()(1)

    def factorial(self):
        return self.parent()(Integer(self).factorial())

    def floor(self):
        return self

    def ceil(self):
        return self

    def squarefree_part(self):
        return self.parent()(Integer(self).squarefree_part())

    def next_prime(self):
        return self.parent()(Integer(self).next_prime())

    def isqrt(self):
        return self.parent()(Integer.isqrt(self))

    def _xgcd(self, n):
        a, b, c = Integer._xgcd(self, n)
        return self.parent()(a), self.parent()(b), self.parent()(c)

    def __lshift__(self, n):
        return self.parent()(Integer(self).__lshift__(n))

    def __rshift__(self, n):
        return self.parent()(Integer(self).__rshift__(n))

    def __and__(self, n):
        if isinstance(n, InfinityElement):
            return self
        return self.parent()(Integer.__and__(self, n))

    def __or__(self, n):
        if isinstance(n, InfinityElement):
            return n
        return self.parent()(Integer.__or__(self, n))

    def inverse_mod(self, n):
        if isinstance(n, InfinityElement):
            return self.__invert__()
        return self.parent()(Integer.inverse_mod(self, n))

class IntegerPlusInfinity(_uniq1, PlusInfinityElement):
    def __init__(self):
        PlusInfinityElement.__init__(self, ExtendedIntegerRing)

    def __cmp__(self, other):
        if isinstance(other, IntegerPlusInfinity):
            return 0
        return 1

    def __repr__(self):
        return "+Infinity"

    def _latex_(self):
        return "+\\infty"

    def _add_(self, other):
        if isinstance(other, IntegerMinusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _mul_(self, other):
        if other < 0:
            return -self
        if other > 0:
            return self
        raise TypeError, "cannot multiply infinity by zero"

    def _sub_(self, other):
        if isinstance(other, IntegerPlusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _div_(self, other):
        return self * other.__invert__()

    def _neg_(self):
        return self.parent().gen(2)

    def __invert__(self):
        return ExtendedInteger(0)

    def __abs__(self):
        return self

    def __pow__(self, right):
        if not isinstance(right, (int, long, Integer, Rational, PlusInfinityElement, MinusInfinityElement)):
            raise TypeError, "cannot exponentiate"
        if right < 0:
            return ExtendedInteger(0)
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
        return self

    def ceil(self):
        return self

    def numerator(self):
        return self

    def denominator(self):
        return self.parent()(1)

    def __xor__(self, other):
        if isinstance(other, InfinityElement):
            return self.parent()(0)
        if other < 0:
            return -self
        return self

    def gcd(self, other):
        if isintance(other, InfinityElement) or other == 0:
            return self
        return other.__abs__()

    def __floordiv__(self, other):
        """
        EXAMPLES:
            sage: inf = ExtendedIntegerRing(Infinity)
            sage: inf // 3
            +Infinity
            sage: inf // -3
            -Infinity
            sage: inf // inf
            Traceback (most recent call last):
            ...
            ValueError: cannot divide Infinity by Infinity
            sage: inf // ExtendedIntegerRing(-Infinity)
            Traceback (most recent call last):
            ...
            ValueError: cannot divide Infinity by Infinity
        """
        if isinstance(other, InfinityElement):
            raise ValueError, "cannot divide Infinity by Infinity"
        elif other > ZZ(0):
            return IntegerPlusInfinity()
        elif other < ZZ(0):
            return IntegerMinusInfinity()
        else:
            raise TypeError, "cannot divide Infinity by object of type %s"%type(other)

    def __mod__(self, right):
        raise ValueError, "remainder not well defined"

    def quo_rem(self, other):
        raise ValueError, "remainder not well defined"

    def powermod(self, exp, mod):
        raise ValueError, "remainder not well defined"

    def coprime_integers(self, m):
        return [self.parent()(1)]

    def divides(self, m):
        return isinstance(m, InfinityElement)

    def factorial(self):
        return self

    def squarefree_part(self):
        raise ValueError, "square free part not defined"

    def next_prime(self):
        raise ValueError, "no primes after infinity"

    def isqrt(self):
        return self

    def _xgcd(self, n):
        raise ValueError, "xgcd not well defined"

    def __lshift__(self, n):
        if isinstance(n, PlusInfinityElement):
            raise SignError, "infinite shift not well defined"
        return self

    def __rshift__(self, n):
        if isinstance(n, MinusInfinityElement):
            raise SignError, "infinite shift not well defined"
        return self

    def __and__(self, n):
        return n.__abs__()

    def __or__(self, n):
        if n < 0:
            return -self
        return self

    def inverse_mod(self, n):
        raise ValueError, "modular inverse not well defined"

class IntegerMinusInfinity(_uniq2, MinusInfinityElement):
    def __init__(self):
        InfinityElement.__init__(self, ExtendedIntegerRing)

    def __cmp__(self, other):
        if isinstance(other, IntegerMinusInfinity):
            return 0
        return -1

    def _repr_(self):
        return "-Infinity"

    def _latex_(self):
        return "-\\infty"

    def _add_(self, other):
        if isinstance(other, IntegerPlusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _mul_(self, other):
        if other < 0:
            return -self
        if other > 0:
            return self
        raise SignError, "cannot multiply infinity by zero"

    def _sub_(self, other):
        if isinstance(other, IntegerMinusInfinity):
            raise SignError, "cannot add infinity to minus infinity"
        return self

    def _div_(self, other):
        return self * other.__invert__()

    def _neg_(self):
        return self.parent().gen(1)

    def __invert__(self):
        return ExtendedInteger(0)

    def __abs__(self):
        return self.parent().gen(1)

    def __pow__(self, right):
        if isinstance(right, Rational):
            if right.denominator() % 2 == 0:
                raise SignError, "Cannot take an even root of negative infinity"
        elif not isinstance(right, (int, long, Integer, PlusInfinityElement, MinusInfinityElement)):
            raise TypeError, "cannot exponentiate"
        if right < 0:
            return ExtendedInteger(0)
        elif right > 0:
            return self
        else:
            raise SignError, "Cannot raise negative infinity to the zeroth power"

    def lcm(self, x):
        """
        Return the least common multiple of -oo and x, which
        is by definition oo unless x is 0.

        EXAMPLES:
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

    def nth_root(self, n):
        if n % 2 == 0:
            raise SignError, "cannot take an even root of negative infinity"
        return self

    def floor(self):
        return self

    def ceil(self):
        return self

    def numerator(self):
        return self

    def denominator(self):
        return self.parent()(1)

    def __xor__(self, other):
        if isinstance(other, InfinityElement):
            return self.parent()(0)
        if other < 0:
            return -self
        return self

    def gcd(self, other):
        if isintance(other, InfinityElement) or other == 0:
            return self.parent().gen(1)
        return other.__abs__()

    def __floordiv__(self, other):
        """
        EXAMPLES:
            sage: inf = ExtendedIntegerRing(-Infinity)
            sage: inf // 3
            -Infinity
            sage: inf // -3
            +Infinity
            sage: inf // inf
            Traceback (most recent call last):
            ...
            ValueError: cannot divide Infinity by Infinity
            sage: inf // ExtendedIntegerRing(Infinity)
            Traceback (most recent call last):
            ...
            ValueError: cannot divide Infinity by Infinity
        """
        if isinstance(other, InfinityElement):
            raise ValueError, "cannot divide Infinity by Infinity"
        elif other > ZZ(0):
            return IntegerMinusInfinity()
        elif other < ZZ(0):
            return IntegerPlusInfinity()
        else:
            raise TypeError, "cannot divide Infinity by object of type %s"%type(other)

    def __mod__(self, right):
        raise ValueError, "remainder not well defined"

    def quo_rem(self, other):
        raise ValueError, "remainder not well defined"

    def powermod(self, exp, mod):
        raise ValueError, "remainder not well defined"

    def coprime_integers(self, m):
        return [self.parent()(1)]

    def divides(self, m):
        return isinstance(m, InfinityElement)

    def factorial(self):
        raise ValueError, "-Infinity! not defined"

    def squarefree_part(self):
        raise ValueError, "square free part not defined"

    def next_prime(self):
        raise ValueError, "no primes after infinity"

    def isqrt(self):
        raise SignError, "square root not defined"

    def _xgcd(self, n):
        raise ValueError, "xgcd not well defined"

    def __lshift__(self, n):
        if isinstance(n, MinusInfinityElement):
            raise SignError, "infinite shift not well defined"
        return self

    def __rshift__(self, n):
        if isinstance(n, PlusInfinityElement):
            raise SignError, "infinite shift not well defined"
        return self

    def __and__(self, n):
        return n

    def __or__(self, n):
        return self

    def inverse_mod(self, n):
        raise ValueError, "modular inverse not well defined"
