r"""
Rational Numbers

AUTHORS:
    -- William Stein (2005): first version
    -- William Stein (2006-02-22): floor and ceil (pure fast GMP versions).
    -- Gonzalo Tornaria and William Stein (2006-03-02): greatly improved
                    python/GMP conversion; hashing
    -- William Stein and Naqi Jaffery (2006-03-06): height, sqrt examples,
          and improve behavior of sqrt.
    -- David Harvey (2006-09-15): added nth_root
"""


###########################################################################
#       Copyright (C) 2004, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

include "../ext/interrupt.pxi"  # ctrl-c interrupt block support
include "../ext/gmp.pxi"
include "../ext/stdsage.pxi"


import operator

from sage.misc.mathml import mathml

import sage.misc.misc as misc
import sage.rings.rational_field
import sage.libs.pari.all

cimport integer
import integer

from sage.structure.element cimport RingElement, ModuleElement
from sage.structure.element import bin_op

import sage.rings.real_mpfr

cimport sage.ext.arith
import  sage.ext.arith

cdef sage.ext.arith.arith_int ai
ai = sage.ext.arith.arith_int()

cdef extern from "../ext/mpz_pylong.h":
    cdef mpz_get_pylong(mpz_t src)
    cdef int mpz_set_pylong(mpz_t dst, src) except -1
    cdef long mpz_pythonhash(mpz_t src)

cdef class Rational(sage.structure.element.FieldElement)

cdef public void set_from_mpq(Rational self, mpq_t value):
    mpq_set(self.value, value)

cdef public void set_from_Rational(Rational self, Rational other):
    mpq_set(self.value, other.value)

cdef public void set_from_Integer(Rational self, integer.Integer other):
    mpq_set_z(self.value, other.value)

cdef object Rational_mul_(Rational a, Rational b):
    cdef Rational x
    x = <Rational> PY_NEW(Rational)

    _sig_on
    mpq_mul(x.value, a.value, b.value)
    _sig_off

    return x

cdef object Rational_div_(Rational a, Rational b):
    cdef Rational x
    x = <Rational> PY_NEW(Rational)

    _sig_on
    mpq_div(x.value, a.value, b.value)
    _sig_off

    return x

cdef Rational_add_(Rational self, Rational other):
    cdef Rational x
    x = <Rational> PY_NEW(Rational)
    _sig_on
    mpq_add(x.value, self.value, other.value)
    _sig_off
    return x

cdef Rational_sub_(Rational self, Rational other):
    cdef Rational x
    x = <Rational> PY_NEW(Rational)

    _sig_on
    mpq_sub(x.value, self.value, other.value)
    _sig_off

    return x

cdef object the_rational_ring
the_rational_ring = sage.rings.rational_field.Q

cdef class Rational(sage.structure.element.FieldElement):
    """
    A Rational number.

    Rational numbers are implemented using the GMP C library.

    EXAMPLES:
        sage: a = -2/3
        sage: type(a)
        <type 'sage.rings.rational.Rational'>
        sage: parent(a)
        Rational Field
        sage: Rational('1/0')
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 1/0 to a rational
        sage: Rational(1.5)
        3/2
        sage: Rational('9/6')
        3/2
    """
    def __new__(self, x=None, int base=0):
        global the_rational_ring
        mpq_init(self.value)
        self._parent = the_rational_ring

    def __init__(self, x=None, int base=0):
        if not (x is None):
            self.__set_value(x, base)

    def __reduce__(self):
        return sage.rings.rational.make_rational, (self.str(32),)

    def __index__(self):
        """
        Needed so integers can be used as list indices.

        EXAMPLES:
            sage: v = [1,2,3,4,5]
            sage: v[3/1]
            4
            sage: v[3/2]
            Traceback (most recent call last):
            ...
            TypeError: rational is not an integer
        """
        if self.denominator() == 1:
            return int(self)
        raise TypeError, "rational is not an integer"

    def _reduce_set(self, s):
        mpq_set_str(self.value, s, 32)

    def __set_value(self, x, int base):
        cdef int n

        if isinstance(x, Rational):
            set_from_Rational(self, x)

        elif isinstance(x, int):
            i = x
            mpq_set_si(self.value, i, 1)

        elif isinstance(x, long):
            mpz_set_pylong(mpq_numref(self.value), x)

        elif isinstance(x, integer.Integer):
            set_from_Integer(self, x)

        elif isinstance(x, sage.rings.real_mpfr.RealNumber):

            if x == 0:
                mpq_set_si(self.value, 0, 1)
                return
            if not base:
                v = x._pari_().contfrac().contfracpnqn()
                s = '%s/%s'%(hex(v[0,0]), hex(v[1,0]))
                n = mpq_set_str(self.value, s, 16)
                # this mpq_denref business is to program around a bug in GMP, where
                # it doesn't correctly return an error code in mpq_set_str (or canonicalalize)
                # for strings with a 0 in the denominator.   In fact, it crashes horribly.
                if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                    raise TypeError, "unable to convert %s to a rational"%x
                mpq_canonicalize(self.value)
            else:
                xstr = x.str(base)
                if '.' in xstr:
                    exp = (len(xstr) - (xstr.index('.') +1))
                    p = base**exp
                    pstr = '1'+'0'*exp
                    s = xstr.replace('.','') +'/'+pstr
                    n = mpq_set_str( self.value, s, base)
                    if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                        raise TypeError, "unable to convert %s to a rational"%x
                    mpq_canonicalize(self.value)
                else:
                    n = mpq_set_str(self.value, xstr, base)
                    if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                        raise TypeError, "unable to convert %s to a rational"%x
                    mpq_canonicalize(self.value)

        elif isinstance(x, str):
            n = mpq_set_str(self.value, x, base)
            if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                raise TypeError, "unable to convert %s to a rational"%x
            mpq_canonicalize(self.value)

        elif hasattr(x, "_rational_"):
            set_from_Rational(self, x._rational_())

        elif isinstance(x, tuple) and len(x) == 2:
            s = "%s/%s"%x
            if x[1] == 0:
                raise ValueError, "denominator must not be 0"
            n = mpq_set_str(self.value, s, 0)
            if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                raise TypeError, "unable to convert %s to a rational"%s
            mpq_canonicalize(self.value)

        elif isinstance(x, list) and len(x) == 1:
            self.__set_value(x[0], base)

        elif isinstance(x, sage.libs.pari.all.pari_gen):
            # TODO: figure out how to convert to pari integer in base 16
            s = str(x)
            n = mpq_set_str(self.value, s, 0)
            if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                raise TypeError, "Unable to coerce %s (%s) to Rational"%(x,type(x))

        else:

            raise TypeError, "Unable to coerce %s (%s) to Rational"%(x,type(x))

    cdef void set_from_mpq(Rational self, mpq_t value):
        mpq_set(self.value, value)


    def list(self):
        """
        Return a list with the rational element in it, to be
        compatible with the method for number fields.

        EXAMPLES:
        sage: m = 5/3
        sage: m.list()
        [5/3]
        """
        return [ self ]


    def __richcmp__(left, right, int op):
        """
        EXAMPLES:
            sage: 1/3 < 2/3
            True
            sage: 2/3 < 1/3
            False
            sage: 4/5 < 2.0
            True
            sage: 4/5 < 0.8
            False
        """
        return (<sage.structure.element.Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, sage.structure.element.Element right) except -2:
        cdef int i
        i = mpq_cmp((<Rational>left).value, (<Rational>right).value)
        if i < 0: return -1
        elif i == 0: return 0
        else: return 1

    def copy(self):
        cdef Rational z
        z = <Rational> PY_NEW(Rational)
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
            t = t + '<mfrac><mrow>%s</mrow><mrow>%s</mrow></mfrac>'%(
                mathml(abs(self.numer())), mathml(self.denom()))
            return t

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
            3.39116499156263
            sage: x = 32/5
            sage: x.sqrt()
            2.52982212813470
            sage: x = 16/9
            sage: x.sqrt()
            4/3
            sage: x.sqrt(53)
            1.33333333333333
            sage: x = 9837/2
            sage: x.sqrt()
            70.1320183653658
            sage: x = 645373/45
            sage: x.sqrt()
            119.756512233040
            sage: x = -12/5
            sage: x.sqrt()
            1.54919333848297*I

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
            R = sage.rings.real_mpfr.RealField(bits)
            return R(self).sqrt()

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

    def nth_root(self, int n):
        r"""
        Computes the nth root of self, or raises a \exception{ValueError}
        if self is not a perfect nth power.

        INPUT:
            n -- integer (must fit in C int type)

        AUTHOR:
           -- David Harvey (2006-09-15)

        EXAMPLES:
          sage: (25/4).nth_root(2)
          5/2
          sage: (125/8).nth_root(3)
          5/2
          sage: (-125/8).nth_root(3)
          -5/2
          sage: (25/4).nth_root(-2)
          2/5

          sage: (9/2).nth_root(2)
          Traceback (most recent call last):
          ...
          ValueError: not a perfect nth power

          sage: (-25/4).nth_root(2)
          Traceback (most recent call last):
          ...
          ValueError: cannot take even root of negative number

        """
        # TODO -- this could be quicker, by using GMP directly.
        cdef integer.Integer num
        cdef integer.Integer den
        cdef int negative

        if n > 0:
            negative = 0
        elif n < 0:
            n = -n
            negative = 1
        else:
            raise ValueError, "n cannot be zero"

        num, exact = self.numerator().nth_root(n, 1)
        if not exact:
            raise ValueError, "not a perfect nth power"

        den, exact = self.denominator().nth_root(n, 1)
        if not exact:
            raise ValueError, "not a perfect nth power"

        if negative:
            return den / num
        else:
            return num / den


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
        k = <object> PyString_FromString(s)
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
        if valid != 0 or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
            mpq_set_si(self.value, 0, 1)  # so data is valid -- but don't waste time making backup.
            raise ValueError, "invalid literal (%s); object set to 0"%s
        mpq_canonicalize(self.value)

    ################################################################
    # Optimized arithmetic
    ################################################################
    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_add(x.value, self.value, (<Rational>right).value)
        return x

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        # self and right are guaranteed to be Integers
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_sub(x.value, self.value, (<Rational>right).value)
        return x

    cdef ModuleElement _neg_c_impl(self):
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_neg(x.value, self.value)
        return x

    cdef RingElement _mul_c_impl(self, RingElement right):
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        if mpz_sizeinbase (mpq_numref(self.value), 2)  > 100000 or \
             mpz_sizeinbase (mpq_denref(self.value), 2) > 100000:
            # We only use the signal handler (to enable ctrl-c out) in case
            # self is huge, so the product might actually take a while to compute.
            _sig_on
            mpq_mul(x.value, self.value, (<Rational>right).value)
            _sig_off
        else:
            mpq_mul(x.value, self.value, (<Rational>right).value)
        return x

    cdef RingElement _div_c_impl(self, RingElement right):
        """
        EXAMPLES:
            sage: 2/3
            2/3
            sage: 3/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Rational division by zero
        """
        if mpq_cmp_si((<Rational> right).value, 0, 1) == 0:
            raise ZeroDivisionError, "Rational division by zero"
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_div(x.value, self.value, (<Rational>right).value)
        return x

    ################################################################
    # Other arithmetic operations.
    ################################################################

    def __invert__(self):
        if self.is_zero():
            raise ZeroDivisionError, "rational division by zero"
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_inv(x.value, self.value)
        return x

    def __pow__(self, n, dummy):
        """
        Raise self to the integer power n.

        EXAMPLES:
            sage: (2/3)^5
            32/243
            sage: (-1/1)^(1/3)
            Traceback (most recent call last):
            ...
            TypeError: exponent (=1/3) must be an integer.
            Coerce your numbers to real or complex numbers first.
        """
        cdef Rational _self, x
        if not isinstance(self, Rational):
            return self.__pow__(float(n))
        _self = self
        if n < 0:
            x = _self**(-n)
            return x.__invert__()
        cdef unsigned int _n
        try:
            _n = integer.Integer(n)
        except TypeError:
            raise TypeError, "exponent (=%s) must be an integer.\nCoerce your numbers to real or complex numbers first."%n
        x = <Rational> PY_NEW(Rational)
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
        x = <Rational> PY_NEW(Rational)
        mpq_neg(x.value, self.value)
        return x

    def __nonzero__(self):
        return not self.numerator().is_zero()

    def __abs__(self):
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_abs(x.value, self.value)
        return x

    def mod_ui(Rational self, unsigned long int n):
        cdef unsigned int num, den, a

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

    def norm(self):
        """
        Returns the norm from Q to Q of x (which is just x). This was
        added for compatibility with NumberFields.

        EXAMPLES:
            sage: (1/3).norm()
             1/3

        AUTHOR:
          -- Craig Citro
        """
        return self

    def trace(self):
        """
        Returns the trace from Q to Q of x (which is just x). This was
        added for compatibility with NumberFields.

        EXAMPLES:
            sage: (1/3).trace()
             1/3

        AUTHOR:
          -- Craig Citro
        """
        return self

    def charpoly(self, var):
        """
        Return the characteristic polynomial of this rational number.
        This will always be just x - self; this is really here
        so that code written for number fields won't crash when
        applied to rational numbers.

        EXAMPLES:
            sage: (1/3).charpoly('x')
             x - 1/3

        AUTHOR:
          -- Craig Citro
        """
        QQ = self.parent()
        return QQ[var]([-self,1])

    def minpoly(self):
        """
        Return the minimal polynomial of this rational number.
        This will always be just x - self; this is really here
        so that code written for number fields won't crash when
        applied to rational numbers.

        EXAMPLES:
            sage: (1/3).minpoly()
             x - 1/3

        AUTHOR:
          -- Craig Citro
        """
        QQ = self.parent()
        return QQ['x']([-self,1])

    cdef integer.Integer _integer_c(self):
        if not mpz_cmp_si(mpq_denref(self.value), 1) == 0:
            raise TypeError, "no coercion of this rational to integer"
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
        n.set_from_mpz(mpq_numref(self.value))
        return n

    def _integer_(self):
        return self._integer_c()

    def numer(self):
        """
        Return the numerator of this rational number.

        EXAMPLE:
            sage: x = -5/11
            sage: x.numer()
            -5
        """
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
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
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
        n.set_from_mpz(mpq_numref(self.value))
        return n


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
        n = PY_NEW(integer.Integer)
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
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
        n.set_from_mpz(mpq_denref(self.value))
        return n

    def factor(self):
        return sage.rings.rational_field.factor(self)

    def floor(self):
        """
        self.floor(): Return the floor of this rational number as an integer.

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

    cdef _lshift(self, unsigned long int exp):
        r"""
        Return $self/2^exp$
        """
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        _sig_on
        mpq_mul_2exp(x.value,self.value,exp)
        _sig_off
        return x

    def __lshift__(x,y):
        if isinstance(x, Rational) and isinstance(y, Rational):
            return x._lshift(y)
        return bin_op(x, y, operator.lshift)

    cdef _rshift(self, unsigned long int exp):
        r"""
        Return $self/2^exp$
        """
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        _sig_on
        mpq_div_2exp(x.value,self.value,exp)
        _sig_off
        return x

    def __rshift__(x,y):
        if isinstance(x, Rational) and isinstance(y, Rational):
            return x._rshift(y)
        return bin_op(x, y, operator.rshift)

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



def pyrex_rational_reconstruction(integer.Integer a, integer.Integer m):
    """
    Find the rational reconstruction of a mod m, if it exists.
    INPUT:
        a -- Integer
        m -- Integer
    OUTPUT:
        x -- rings.rational.Rational
    """
    cdef Rational x
    x = <Rational> PY_NEW(Rational)
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


def make_rational(s):
    r = Rational()
    r._reduce_set(s)
    return r
