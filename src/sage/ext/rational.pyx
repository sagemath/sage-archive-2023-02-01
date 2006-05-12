"""
Rational Numbers

AUTHORS:
    -- William Stein (2005): first version
    -- William Stein (2006-02-22): floor and ceil (pure fast GMP versions).
    -- Gonzalo Tornaria and William Stein (2006-03-02): greatly improved
                    python/GMP conversion; hashing
    -- William Stein and Naqi Jaffery (2006-03-06): height, sqrt examples,
          and improve behavior of sqrt.
"""


###########################################################################
#       Copyright (C) 2004, 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

include "interrupt.pxi"  # ctrl-c interrupt block support

import operator

from sage.misc.mathml import mathml

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

cdef extern from "mpz_pylong.h":
    cdef mpz_get_pylong(mpz_t src)
    cdef int mpz_set_pylong(mpz_t dst, src) except -1
    cdef long mpz_pythonhash(mpz_t src)

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
        """
        EXAMPLES:
            sage: a = long(901824309821093821093812093810928309183091832091)
            sage: b = QQ(a); b
            901824309821093821093812093810928309183091832091
            sage: QQ(b)
            901824309821093821093812093810928309183091832091
            sage: QQ(int(93820984323))
            93820984323
            sage: QQ(ZZ(901824309821093821093812093810928309183091832091))
            901824309821093821093812093810928309183091832091
            sage: QQ('-930482/9320842317')
            -930482/9320842317
            sage: QQ((-930482, 9320842317))
            -930482/9320842317
            sage: QQ([9320842317])
            9320842317
            sage: QQ(pari(39029384023840928309482842098430284398243982394))
            39029384023840928309482842098430284398243982394
            sage: QQ('sage')
            Traceback (most recent call last):
            ...
            TypeError: unable to convert sage to a rational
        """
        cdef int n
        if isinstance(x, Rational):
            set_from_Rational(self, x)

        elif isinstance(x, int):
            i = x
            mpq_set_si(self.value, i, 1)

        elif isinstance(x, long):
            mpz_set_pylong(mpq_numref(self.value), x)
            #s = "%x"%x
            #mpq_set_str(self.value, s, 16)

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

        elif isinstance(x, sage.libs.pari.all.pari_gen):
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

    def _mathml_(self):
        if self.denom() == 1:
            return '<mn>%s</mn>'%(self.numer())
        else:
            t = ''
            if self < 0:
                t = t + '<mo>-</mo>'
            t = + '<mfrac><mrow>%s<mrow>%s</mrow></mfrac>'%(
                mathml(abs(self.numer())), mathml(self.denom()))
            return t

    def _mpfr_(self, R):
        return R(self.numerator()) / R(self.denominator())

    def _im_gens_(self, codomain, im_gens):
        return codomain._coerce_(self)

    def lcm(self, Rational other):
        """
        Return the least common multiple of self and other.

        Our hopefully interesting notion of LCM for rational numbers
        is illustrated in the examples below.

        EXAMPLES:
            sage: lcm(2/3,1/5)
            2
            sage: lcm(2/3,7/5)
            14
            sage: lcm(1/3,1/5)
            1
            sage: lcm(1/3,1/6)
            1/3
        """
        d = self.denom()*other.denom()
        self_d = self.numer()*other.denom()
        other_d = other.numer()*self.denom()
        return self_d.lcm(other_d) / d

    def gcd(self, Rational other):
        """
        Return the least common multiple of self and other.

        Our hopefully interesting notion of GCD for rational numbers
        is illustrated in the examples below.

        EXAMPLES:
            sage: gcd(2/3,1/5)
            1/15
            sage: gcd(2/3,7/5)
            1/15
            sage: gcd(1/3,1/6)
            1/6
            sage: gcd(6/7,9/7)
            3/7
        """
        d = self.denom()*other.denom()
        self_d = self.numer()*other.denom()
        other_d = other.numer()*self.denom()
        return self_d.gcd(other_d) / d

    def valuation(self, p):
        return self.numerator().valuation(p) - self.denominator().valuation(p)

    def sqrt(self, bits=None):
        r"""
        Returns the positive square root of self as a real number to
        the given number of bits of precision if self is nonnegative,
        and raises a \exception{ValueError} exception otherwise.

        INPUT:
            bits -- number of bits of precision.
                    If bits is not specified, the number of
                    bits of precision is at least twice the
                    number of bits of self (the precision
                    is always at least 53 bits if not specified).
        OUTPUT:
            integer, real number, or complex number.

        EXAMPLES:
            sage: x = 23/2
            sage: x.sqrt()
            3.3911649915626341
            sage: x = 32/5
            sage: x.sqrt()
            2.5298221281347035
            sage: x = 16/9
            sage: x.sqrt()
            4/3
            sage: x.sqrt(53)
            1.3333333333333333
            sage: x = 9837/2
            sage: x.sqrt()
            70.132018365365752
            sage: x = 645373/45
            sage: x.sqrt()
            119.75651223303984
            sage: x = -12/5
            sage: x.sqrt()
            1.5491933384829668*I

        AUTHOR:
            -- Naqi Jaffery (2006-03-05): examples
        """
        if bits is None:
            try:
                return self.square_root()
            except ValueError:
                pass
            bits = max(53, 2*(mpz_sizeinbase(self.value, 2)+2))

        if self < 0:
            x = sage.rings.complex_field.ComplexField(bits)(self)
            return x.sqrt()
        else:
            R = mpfr.RealField(bits)
            return self._mpfr_(R).sqrt()

    def square_root(self):
        r"""
        Return the positive rational square root of self, or raises a
        \exception{ValueError} if self is not a perfect square.

        EXAMPLES:
            sage: x = 125/5
            sage: x.square_root()
            5
            sage: x = 64/4
            sage: x.square_root()
            4
            sage: x = 1000/10
            sage: x.square_root()
            10
            sage: x = 81/3
            sage: x.square_root()
            Traceback (most recent call last):
            ...
            ValueError: self (=27) is not a perfect square

        AUTHOR:
            -- Naqi Jaffery (2006-03-05): examples
        """
        # TODO -- this could be quicker, by using GMP directly.
        return self.numerator().square_root() / self.denominator().square_root()

    def str(self, int base=10):
        if base < 2 or base > 36:
            raise ValueError, "base (=%s) must be between 2 and 36"%base
        cdef size_t n
        cdef char *s

        n = mpz_sizeinbase (mpq_numref(self.value), base) \
            + mpz_sizeinbase (mpq_denref(self.value), base) + 3
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
        cdef long n, d
        n = mpz_pythonhash(mpq_numref(self.value))
        d = mpz_pythonhash(mpq_denref(self.value))
        if d == 1:
            return n
        n = n ^ d
        if n == -1:
            return -2
        return n

##         cdef char *s
##         cdef int h

##         s = mpz_get_str(NULL, 16, mpq_numref(self.value))
##         h = hash(long(s,16))
##         free(s)
##         if mpz_cmp_si(mpq_denref(self.value), 1) == 0:
##             return h
##         else:
##             s = mpz_get_str(NULL, 16, mpq_denref(self.value))
##             h = h ^ hash(long(s,16))     # xor
##             if h == -1:    # -1 is not a valid return value
##                 h = -2
##             free(s)
##         return h

        #cdef int n
        #n = mpz_get_si(mpq_numref(self.value)) * \
        #    mpz_get_si(mpq_denref(self.value))
        #if n == -1:
        #    return -2     # since -1 is not an allowed Python hash for C ext -- it's an error indicator.


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
        """
        Return the numerator of this rational number.

        EXAMPLE:
            sage: x = -5/11
            sage: x.numer()
            -5
        """
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
        """
        Return coercion of self to Python int.

        This takes the floor of self if self has a denominator (which
        is consistent with Python's long(floats)).

        EXAMPLES:
            sage: int(7/3)
            2
            sage: int(-7/3)
            -3
        """
        return int(self.__long__())

    def __long__(self):
        """
        Return coercion of self to Python long.

        This takes the floor of self if self has a denominator (which
        is consistent with Python's long(floats)).

        EXAMPLES:
            sage: long(7/3)
            2L
            sage: long(-7/3)
            -3L
        """
        cdef mpz_t x
        if mpz_cmp_si(mpq_denref(self.value),1) != 0:
            mpz_init(x)
            mpz_fdiv_q(x, mpq_numref(self.value), mpq_denref(self.value))
            n = mpz_get_pylong(x)
            mpz_clear(x)
            return n
        else:
            return mpz_get_pylong(mpq_numref(self.value))

    def denom(self):
        """
        self.denom(): Return the denominator of this rational number.

        EXAMPLES:
            sage: x = 5/13
            sage: x.denom()
            13
            sage: x = -9/3
            sage: x.denom()
            1
        """
        cdef integer.Integer n
        n = integer.Integer()
        n.set_from_mpz(mpq_denref(self.value))
        return n

    def denominator(self):
        """
        self.denominator(): Return the denominator of this rational number.

        EXAMPLES:
            sage: x = -5/11
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
        """
        self.floor(): Return the floor of this rational number.

        EXAMPLES:
            sage: n = 5/3; n.floor()
            1
            sage: n = -17/19; n.floor()
            -1
            sage: n = -7/2; n.floor()
            -4
            sage: n = 7/2; n.floor()
            3
            sage: n = 10/2; n.floor()
            5
        """
        cdef integer.Integer n
        n = integer.Integer()
        mpz_fdiv_q(n.value, mpq_numref(self.value), mpq_denref(self.value))
        return n

    def ceil(self):
        """
        self.ceil(): Return the ceiling of this rational number.

        If this rational number is an integer, this returns this
        number, otherwise it returns the floor of this number +1.

        EXAMPLES:
            sage: n = 5/3; n.ceil()
            2
            sage: n = -17/19; n.ceil()
            0
            sage: n = -7/2; n.ceil()
            -3
            sage: n = 7/2; n.ceil()
            4
            sage: n = 10/2; n.ceil()
            5
        """
        cdef integer.Integer n
        n = integer.Integer()
        mpz_cdiv_q(n.value, mpq_numref(self.value), mpq_denref(self.value))
        return n

    def height(self):
        """
        The max absolute value of the numerator and denominator of self,
        as an Integer.

        EXAMPLES:
            sage: a = 2/3
            sage: a.height()
            3
            sage: a = 34/3
            sage: a.height()
            34
            sage: a = -97/4
            sage: a.height()
            97

        AUTHOR:
            -- Naqi Jaffery (2006-03-05): examples
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

    def _lshift(self, unsigned long int exp):
        r"""
        Return $self/2^exp$
        """
        cdef Rational x
        x = Rational()
        _sig_on
        mpq_mul_2exp(x.value,self.value,exp)
        _sig_off
        return x

    def __lshift__(x,y):
        if isinstance(x, Rational) and isinstance(y, Rational):
            return x._lshift(y)
        return sage.rings.coerce.bin_op(x, y, operator.lshift)

    def _rshift(self, unsigned long int exp):
        r"""
        Return $self/2^exp$
        """
        cdef Rational x
        x = Rational()
        _sig_on
        mpq_div_2exp(x.value,self.value,exp)
        _sig_off
        return x

    def __rshift__(x,y):
        if isinstance(x, Rational) and isinstance(y, Rational):
            return x._rshift(y)
        return sage.rings.coerce.bin_op(x, y, operator.rshift)

    ##################################################
    # Support for interfaces
    ##################################################

    def _pari_(self):
        return self.numerator()._pari_()/self.denominator()._pari_()

    def _interface_init_(self):
        """
        EXAMPLES:
            sage: kash(3/1).Type()              # optional
            elt-fld^rat
            sage: magma(3/1).Type()             # optional
            FldRatElt
        """
        return '%s/%s'%(self.numerator(), self.denominator())



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
