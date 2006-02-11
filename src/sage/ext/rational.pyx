"""
Rational Numbers
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

include "interrupt.pxi"  # ctrl-c interrupt block support

import operator

import sage.misc.misc as misc
import sage.rings.rational_field
import sage.rings.coerce
import sage.libs.pari.all

cimport integer
import integer

import mpfr

cimport arith
import arith
cdef arith.arith_int ai
ai = arith.arith_int()

include "gmp.pxi"

cdef class Rational(element.FieldElement)

cdef public void set_from_mpq(Rational self, mpq_t value):
    mpq_set(self.value, value)

cdef public void set_from_Rational(Rational self, Rational other):
    mpq_set(self.value, other.value)

cdef public void set_from_Integer(Rational self, integer.Integer other):
    mpq_set_z(self.value, other.value)

cdef object Rational_mul_(Rational a, Rational b):
    cdef Rational x
    x = Rational()

    _sig_on
    mpq_mul(x.value, a.value, b.value)
    _sig_off

    return x

cdef object Rational_div_(Rational a, Rational b):
    cdef Rational x
    x = Rational()

    _sig_on
    mpq_div(x.value, a.value, b.value)
    _sig_off

    return x

cdef Rational_add_(Rational self, Rational other):
    cdef Rational x
    x = Rational()
    _sig_on
    mpq_add(x.value, self.value, other.value)
    _sig_off
    return x

cdef Rational_sub_(Rational self, Rational other):
    cdef Rational x
    x = Rational()

    _sig_on
    mpq_sub(x.value, self.value, other.value)
    _sig_off

    return x

cdef class Rational(element.FieldElement):

    def __new__(self, x=None):
        mpq_init(self.value)

    def __init__(self, x=None):
        if not (x is None):
            self.__set_value(x)

    def parent(self):
        return sage.rings.rational_field.Q

    def __reduce__(self):
        return sage.rings.rational.make_rational, (self.str(32),)

    def _reduce_set(self, s):
        mpq_set_str(self.value, s, 32)

    def __set_value(self, x):
        cdef int n
        if isinstance(x, Rational):
            set_from_Rational(self, x)
        elif isinstance(x, int):
            i = x
            mpq_set_si(self.value, i, 1)

        elif isinstance(x, long):
            s = "%x"%x
            mpq_set_str(self.value, s, 16)

        elif isinstance(x, integer.Integer):
            set_from_Integer(self, x)

        elif isinstance(x, str):
            if not mpq_set_str(self.value, x, 0) == 0:
                raise TypeError, "unable to convert %s to a rational"%x
            mpq_canonicalize(self.value)

        elif hasattr(x, "_rational_"):
            set_from_Rational(self, x._rational_())

        elif isinstance(x, tuple) and len(x) == 2:
            s = "%s/%s"%x
            if not mpq_set_str(self.value, s, 0) == 0:
                raise TypeError, "unable to convert %s to a rational"%s
            mpq_canonicalize(self.value)

        elif isinstance(x, list) and len(x) == 1:

            self.__set_value(x[0])

        elif isinstance(x, sage.libs.pari.all.gen):
            # TODO: figure out how to convert to pari integer in base 16
            s = str(x)
            if mpq_set_str(self.value, s, 0):
                raise TypeError, "Unable to coerce %s (%s) to Rational"%(x,type(x))

        else:
            raise TypeError, "Unable to coerce %s (%s) to Rational"%(x,type(x))

    cdef void set_from_mpq(Rational self, mpq_t value):
        mpq_set(self.value, value)

    cdef cmp(Rational self, Rational x):
        cdef int i
        i = mpq_cmp(self.value, x.value)
        if i < 0:
            return -1
        elif i == 0:
            return 0
        else:
            return 1

    def __cmp__(self, x):
        return self.cmp(x)

    def __richcmp__(Rational self, right, int op):
        cdef int n
        if not isinstance(right, Rational):
            try:
                n = sage.rings.coerce.cmp(self, right)
            except TypeError:
                n = -1
        else:
            n = self.cmp(right)
        return self._rich_to_bool(op, n)

    def _pari_(self):
        return self.numerator()._pari_()/self.denominator()._pari_()


    def copy(self):
        cdef Rational z
        z = Rational()
        mpq_set(z.value, self.value)
        return z

    def  __dealloc__(self):
        mpq_clear(self.value)

    def __repr__(self):
        return self.str()

    def _latex_(self):
        if self.denom() == 1:
            return str(self.numer())
        else:
            if self < 0:
                return "-\\frac{%s}{%s}"%(-self.numer(), self.denom())
            else:
               return "\\frac{%s}{%s}"%(self.numer(), self.denom())

    def _mpfr_(self, R):
        return R(self.numerator()) / R(self.denominator())

    def _im_gens_(self, codomain, im_gens):
        return codomain._coerce_(self)

    def valuation(self, p):
        return self.numerator().valuation(p) - self.denominator().valuation(p)

    def sqrt(self, bits=53):
        """
        Returns the positive square root of self as a real number to
        the given number of bits of precision if self is nonnegative,
        and raises a \\exception{ValueError} exception otherwise.
        """
        R = mpfr.RealField(bits)
        return self._mpfr_(R).sqrt()

    def square_root(self):
        """
        Return the positive rational square root of self, or raises a
        ValueError if self is not a perfect square.
        """
        # TODO -- this could be quicker maybe, by using GMP directly.
        return self.numerator().square_root() / self.denominator().square_root()

    def str(self, int base=10):
        if base < 2 or base > 36:
            raise ValueError, "base (=%s) must be between 2 and 36"%base
        cdef size_t n
        cdef char *s

        n = mpz_sizeinbase (mpq_numref(self.value), base) \
            + mpz_sizeinbase (mpq_denref(self.value), base) + 3
        if n > 4200000 and base != 2 and base != 4 and base != 8 and base != 16 and base != 32:
            raise RuntimeError, "String representation of rationals with more than 4200000 digits is not available in bases other than a power of 2. This is because of a bug in GMP.  To obtain the base-b expansion of x use x.str(b)."
        s = <char *>PyMem_Malloc(n)
        if s == NULL:
            raise MemoryError, "Unable to allocate enough memory for the string representation of an integer."

        _sig_on
        mpq_get_str(s, base, self.value)
        _sig_off
        k = PyString_FromString(s)
        PyMem_Free(s)
        return k

    def __float__(self):
        return mpq_get_d(self.value)

    def __hash__(self):
        return hash(mpq_get_d(self.value))

    def __getitem__(self, int n):
        if n == 0:
            return self
        raise IndexError, "index n (=%s) out of range; it must be 0"%n

    def set_si(self, signed long int n):
        mpq_set_si(self.value, n, 1)

    def set_str(self, s, base=10):
        valid = mpq_set_str(self.value, s, base)
        if valid != 0:
            raise ValueError, "invalid literal:" + s

    def __add_(Rational self, Rational other):
        cdef Rational x
        x = Rational()
        _sig_on
        mpq_add(x.value, self.value, other.value)
        _sig_off
        return x

    def __add__(x, y):
        if isinstance(x, Rational) and isinstance(y, Rational):
            return x.__add_(y)
        return sage.rings.coerce.bin_op(x, y, operator.add)

    def __sub_(Rational self, Rational other):
        cdef Rational x
        x = Rational()
        _sig_on
        mpq_sub(x.value, self.value, other.value)
        _sig_off
        return x

    def __sub__(x, y):
        if isinstance(x, Rational) and isinstance(y, Rational):
            return x.__sub_(y)
        return sage.rings.coerce.bin_op(x, y, operator.sub)

    def __mul_(Rational self, Rational other):
        cdef Rational x
        x = Rational()
        _sig_on
        mpq_mul(x.value, self.value, other.value)
        _sig_off
        return x

    def __mul__(x, y):
        if isinstance(x, Rational) and isinstance(y, Rational):
            return x.__mul_(y)
        return sage.rings.coerce.bin_op(x, y, operator.mul)

    def __div_(Rational self, Rational other):
        if not other:
            raise ZeroDivisionError, "Rational division by zero"
        cdef Rational x
        x = Rational()
        _sig_on
        mpq_div(x.value, self.value, other.value)
        _sig_off
        return x

    def __div__(x, y):
        if isinstance(x, Rational) and isinstance(y, Rational):
            return x.__div_(y)
        return sage.rings.coerce.bin_op(x, y, operator.div)

    def __invert__(self):
        if self.is_zero():
            raise ZeroDivisionError, "rational division by zero"
        cdef Rational x
        x = Rational()
        _sig_on
        mpq_inv(x.value, self.value)
        _sig_off
        return x

    def __pow__(self, n, dummy):
        cdef Rational _self, x
        if not isinstance(self, Rational):
            return self.__pow__(float(n))
        _self = self
        if n < 0:
            x = _self**(-n)
            return x.__invert__()
        cdef unsigned int _n
        _n = n
        x = Rational()
        cdef mpz_t num, den

        _sig_on
        mpz_init(num)
        mpz_init(den)
        mpz_pow_ui(num, mpq_numref(_self.value), _n)
        mpz_pow_ui(den, mpq_denref(_self.value), _n)
        mpq_set_num(x.value, num)
        mpq_set_den(x.value, den)
        mpz_clear(num)
        mpz_clear(den)
        _sig_off

        return x

    def __pos__(self):
        return self

    def __neg__(self):
        cdef Rational x
        x = Rational()
        mpq_neg(x.value, self.value)
        return x

    def __nonzero__(self):
        return not self.numerator().is_zero()

    def __abs__(self):
        cdef Rational x
        x = Rational()
        mpq_abs(x.value, self.value)
        return x

    def mod_ui(Rational self, unsigned long int n):
        cdef unsigned int num, den, a
        cdef Rational x

        # Documentation from GMP manual:
        # "For the ui variants the return value is the remainder, and
        # in fact returning the remainder is all the div_ui functions do."
        _sig_on
        num = mpz_fdiv_ui(mpq_numref(self.value), n)
        den = mpz_fdiv_ui(mpq_denref(self.value), n)
        _sig_off
        return int((num * ai.inverse_mod_int(den, n)) % n)

    def __mod__(Rational self, other):
        other = integer.Integer(other)
        if not other:
            raise ZeroDivisionError, "Rational modulo by zero"
        n = self.numer() % other
        d = self.denom() % other
        _sig_on
        d = d.inverse_mod(other)
        _sig_off
        return (n*d)%other

    def numer(self):
        cdef integer.Integer n
        n = integer.Integer()
        n.set_from_mpz(mpq_numref(self.value))
        return n

    def numerator(self):
        """
        Return the numerator of this rational number.

        EXAMPLE:
            sage: x = 5/11
            sage: x.numerator()
            5

            sage: x = 9/3
            sage: x.numerator()
            3
        """
        return self.numer()

    def __int__(self):
        if mpz_cmp_si(mpq_denref(self.value),1) != 0:
            raise ValueError, "cannot convert %s to an int"%self
        return int(self.numerator())

    def __long__(self):
        return int(self)

    def denom(self):
        cdef integer.Integer n
        n = integer.Integer()
        n.set_from_mpz(mpq_denref(self.value))
        return n

    def denominator(self):
        """
        Return the denominator of this rational number.

        EXAMPLES:
            sage: x = 5/11
            sage: x.denominator()
            11
            sage: x = 9/3
            sage: x.denominator()
            1
        """
        return self.denom()

    def factor(self):
        return sage.rings.rational_field.factor(self)

    def floor(self):
        return self.numer().__floordiv__(self.denom())

    def height(self):
        """
        The max absolute value of the numerator and denominator of self,
        as an Integer.
        """
        x = abs(self.numer())
        if x > self.denom():
            return x
        return self.denom()

    def _lcm(self, Rational other):
        """
        Returns the least common multiple, in the rational numbers,
        of self and other.  This function returns either 0 or 1 (as
        a rational number).
        """
        if mpz_cmp_si(mpq_numref(self.value), 0) == 0 and \
               mpz_cmp_si(mpq_numref(other.value), 0) == 0:
            return Rational(0)
        return Rational(1)

    def _gcd(self, Rational other):
        """
        Returns the least common multiple, in the rational numbers,
        of self and other.  This function returns either 0 or 1 (as
        a rational number).
        """
        if mpz_cmp_si(mpq_numref(self.value), 0) == 0 and \
               mpz_cmp_si(mpq_numref(other.value), 0) == 0:
            return Rational(0)
        return Rational(1)


    def additive_order(self):
        """
        Return the additive order of self.

        EXAMPLES:
            sage: QQ(0).additive_order()
            1
            sage: QQ(1).additive_order()
            Infinity
        """
        import sage.rings.infinity
        if self.is_zero():
            return integer.Integer(1)
        else:
            return sage.rings.infinity.infinity


    def multiplicative_order(self):
        """
        Return the multiplicative order of self, if self is a unit, or raise
        \code{ArithmeticError} otherwise.

        EXAMPLES:
            sage: QQ(1).multiplicative_order()
            1
            sage: QQ('1/-1').multiplicative_order()
            2
            sage: QQ(0).multiplicative_order()
            Traceback (most recent call last):
            ...
            ArithmeticError: no power of 0 is a unit
            sage: QQ('2/3').multiplicative_order()
            Traceback (most recent call last):
            ...
            ArithmeticError: no power of 2/3 is a unit
        """
        if mpz_cmp_si(mpq_numref(self.value), 1) == 0:
            return integer.Integer(1)
        elif mpz_cmp_si(mpq_numref(self.value), -1) == 0:
            return integer.Integer(2)
        else:
            raise ArithmeticError, "no power of %s is a unit"%self

    def is_one(self):
        return mpz_cmp_si(mpq_numref(self.value), 1) == 0

    def is_zero(self):
        return mpz_cmp_si(mpq_numref(self.value), 0) == 0


def rational_reconstruction(integer.Integer a, integer.Integer m):
    """
    Find the rational reconstruction of a mod m, if it exists.
    INPUT:
        a -- Integer
        m -- Integer
    OUTPUT:
        x -- rings.rational.Rational
    """
    cdef Rational x
    #if not isinstance(a, integer.Integer) or not isinstance(m, integer.Integer):
    #    a = integer.Integer(a)
    #    m = integer.Integer(m)
    x = Rational()
    mpq_rational_reconstruction(x.value, a.get_value()[0], m.get_value()[0])
    return x

# def test(n):
#     cdef mpz_t a, m
#     mpz_init(a); mpz_init(m)
#     mpz_set_str(a, "133173434946179363436783721312585228097", 0)
#     mpz_set_str(m, "10000000000000000000000000000000000000000", 0)
#     cdef int k
#     cdef mpq_t x
#     mpq_init(x)
#     import time
#     t = time.clock()
#     for k from 0 <= k < n:
#         mpq_rational_reconstruction(x, a, m)
#     print time.clock() - t
