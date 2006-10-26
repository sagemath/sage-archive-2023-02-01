"""
Real Numbers

AUTHORS: Kyle Schalm <kschalm@math.utexas.edu> (2005-09)
         William Stein <wstein@gmail.com>: bug fixes, examples, maintenance
         Didier Deshommes <dfdeshom@gmail.com> (2006-03-19): examples
         David Harvey (2006-09-20): compatibility with Element._parent
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005-2006 William Stein <wstein@gmail.com>
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

#*****************************************************************************
# general TODOs:
#
# more type conversions and coercion. examples:
# sage: R(1/2)
# TypeError: Unable to convert x (='1/2') to mpfr.
#
# sage: 1 + R(42)
# _49 = 1
#
# ................
# Need to revisit c_parent() business; i.e. need to find a way to do type
# conversions in pyrex as fast as possible while still being compatible with
# Element._parent conventions
#*****************************************************************************

import math # for log
import sys

include '../ext/interrupt.pxi'

cimport sage.rings.ring
import  sage.rings.ring

cimport sage.structure.element
import  sage.structure.element

import sage.rings.coerce
import operator

import sage.rings.rational
import sage.rings.complex_field
import sage.rings.integer

import sage.rings.infinity



#*****************************************************************************
# Headers.  When you past things in here from mpfr, be sure
# to remove const's, since those aren't allowed in pyrex.  Also, it can be
# challenging figuring out how to modify things from mpfr.h to be valid pyrex
# code.    Note that what is here is only used for generating the C code.
# The C compiler doesn't see any of this -- it only sees mpfr.h and stdlib.h
#*****************************************************************************

cdef class RealNumber(sage.structure.element.RingElement)

#*****************************************************************************
#
#       Implementation
#
#*****************************************************************************

#*****************************************************************************
#
#       External Python access to constants
#
#*****************************************************************************

def mpfr_prec_min():
    """
    Return the mpfr variable MPFR_PREC_MIN.
    """
    return MPFR_PREC_MIN

def mpfr_prec_max():
    return MPFR_PREC_MAX

#*****************************************************************************
#
#       Real Field
#
#*****************************************************************************
# The real field is in Pyrex, so mpfr elements will have access to
# their parent via direct C calls, which will be faster.

_rounding_modes = ['RNDN', 'RNDZ', 'RNDU', 'RNDD']

cdef class RealField(sage.rings.ring.Field):
    """
    RealField(prec, sci_not, rnd):

    INPUT:
        prec -- (integer) precision; default = 53
                prec is the number of bits used to represent the
                mantissa of a floating-point number.  The
                precision can be any integer between mpfr_prec_min()
                and mpfr_prec_max(). In the current implementation,
                mpfr_prec_min() is equal to 2.

        sci_not -- (default: True) whether or not to display
                using scientific notation

        rnd -- (string) the rounding mode
                RNDN -- round to nearest  (default)
                RNDZ -- round towards zero
                RNDU -- round towards plus infinity
                RNDD -- round towards minus infinity

    EXAMPLES:
        sage: RealField(10)
        Real Field with 10 bits of precision
        sage: RealField()
        Real Field with 53 bits of precision
        sage: RealField(100000)
        Real Field with 100000 bits of precision

    NOTE: The default precision is 53, since according to the GMP
       manual: "mpfr should be able to exactly reproduce all
       computations with double-precision machine floating-point
       numbers (double type in C), except the default exponent
       range is much wider and subnormal numbers are not
       implemented."
    """

    def __init__(self, int prec=53, int sci_not=0, rnd="RNDN"):
        if prec < MPFR_PREC_MIN or prec > MPFR_PREC_MAX:
            raise ValueError, "prec (=%s) must be >= %s and <= %s."%(
                prec, MPFR_PREC_MIN, MPFR_PREC_MAX)
        self.__prec = prec
        if not isinstance(rnd, str):
            raise TypeError, "rnd must be a string"
        self.sci_not = sci_not
        n = _rounding_modes.index(rnd)
        if n == -1:
            raise ValueError, "rnd (=%s) must be one of RNDN, RNDZ, RNDU, or RNDD"%rnd
        self.rnd = n
        self.rnd_str = rnd

    def _repr_(self):
        s = "Real Field with %s bits of precision"%self.__prec
        if self.rnd != GMP_RNDN:
            s = s + " and rounding %s"%(self.rnd_str)
        return s

    def _latex_(self):
        return "\\R"

    def __call__(self, x, base=10):
        """
        Coerce x into this real field.

        EXAMPLES:
            sage: R = RealField(10)
            sage: R('1.234')
            1.2344
            sage: R('2', base=2)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert x (='2') to real number.
            sage: a = R('1.1001', base=2); a
            1.5625
            sage: a.str(2)
            '1.100100000'
        """
        if hasattr(x, '_mpfr_'):
            return x._mpfr_(self)
        return RealNumber(self, x, base)

    def _coerce_(self, x):
        cdef RealField K
        if isinstance(x, RealNumber):
            K = x.parent()
            if K is self:
                return x
            elif K.__prec >= self.__prec:
                return self(x)
            else:
                raise TypeError
        if isinstance(x, (int, long, sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            return self(x)
        try:
            return x._mpfr_(self)
        except AttributeError:
            pass
        raise TypeError

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: RealField(10) == RealField(11)
            False
            sage: RealField(10) == RealField(10)
            True
            sage: RealField(10,rnd='RNDN') == RealField(10,rnd='RNDZ')
            False
            sage: RealField(10,sci_not=True) == RealField(10,sci_not=False)
            False
            sage: RealField(10) == IntegerRing()
            False
        """
        if not isinstance(other, RealField):
            return -1
        cdef RealField _other
        _other = other  # to access C structure
        if self.__prec == _other.__prec and self.rnd == _other.rnd \
               and self.sci_not == _other.sci_not:
            return 0
        return 1

    def __reduce__(self):
        """
        EXAMPLES:
            sage: R = RealField(sci_not=1, prec=200, rnd='RNDU')
            sage: loads(dumps(R)) == R
            True
        """
        return sage.rings.real_field.__reduce__RealField, \
                (self.__prec, self.sci_not, self.rnd_str)

    def gen(self, i=0):
        if i == 0:
            return self(1)
        else:
            raise IndexError

    def complex_field(self):
        """
        Return complex field of the same precision.
        """
        return sage.rings.complex_field.ComplexField(self.prec())

    def ngens(self):
        return 1

    def gens(self):
        return [self.gen()]

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            s = codomain._coerce_(self(1))
        except TypeError:
            return False
        return s == im_gens[0]

    def is_atomic_repr(self):
        """
        Returns True, to signify that elements of this field print
        without sums, so parenthesis aren't required, e.g., in
        coefficients of polynomials.

        EXAMPLES:
            sage: RealField(10).is_atomic_repr()
            True
        """
        return True

    def is_finite(self):
        """
        Returns False, since the field of real numbers is not finite.

        EXAMPLES:
            sage: RealField(10).is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Returns 0, since the field of real numbers has characteristic 0.

        EXAMPLES:
            sage: RealField(10).characteristic()
            0
        """
        return 0

    def name(self):
        return "RealField%s_%s"%(self.__prec,self.rnd)

    def __hash__(self):
        return hash(self.name())

    def precision(self):
        return self.__prec

    def prec(self):
        return self.__prec

    # int mpfr_const_pi (mpfr_t rop, mp_rnd_t rnd)
    def pi(self):
        """
        Returns pi to the precision of this field.

        EXAMPLES:
            sage: R = RealField(100)
            sage: R.pi()
            3.1415926535897932384626433832793
            sage: R.pi().sqrt()/2
            0.88622692545275801364908374167063
        """
        cdef RealNumber x
        x = RealNumber(self)
        mpfr_const_pi(x.value, self.rnd)
        return x


    # int mpfr_const_euler (mpfr_t rop, mp_rnd_t rnd)
    def euler_constant(self):
        """
        Returns Euler's gamma constant to the precision of this field.

        EXAMPLES:
            sage: RealField(100).euler_constant()
            0.57721566490153286060651209008234
        """
        cdef RealNumber x
        x = RealNumber(self)
        mpfr_const_euler(x.value, self.rnd)
        return x

    # int mpfr_const_catalan (mpfr_t rop, mp_rnd_t rnd)
    def catalan_constant(self):
        """
        Returns Catalan's constant to the precision of this field.

        EXAMPLES:
            sage: RealField(100).catalan_constant()
            0.91596559417721901505460351493252
        """
        cdef RealNumber x
        x = RealNumber(self)
        mpfr_const_catalan(x.value, self.rnd)
        return x

    # int mpfr_const_log2 (mpfr_t rop, mp_rnd_t rnd)
    def log2(self):
        """
        Returns log(2) to the precision of this field.

        EXAMPLES:
            sage: R=RealField(100)
            sage: R.log2()
            0.69314718055994530941723212145798
            sage: R(2).log()
            0.69314718055994530941723212145798
        """
        cdef RealNumber x
        x = RealNumber(self)
        mpfr_const_log2(x.value, self.rnd)
        return x

    def factorial(self, int n):
        """
        Return the factorial of the integer n as a real number.
        """
        cdef RealNumber x
        if n < 0:
            raise ArithmeticError, "n must be nonnegative"
        x = RealNumber(self)
        mpfr_fac_ui(x.value, n, self.rnd)
        return x

    def rounding_mode(self):
        return _rounding_modes[self.rnd]

    def scientific_notation(self, status=None):
        """
        Set or return the scientific notation printing flag.  If this flag
        is True then real numbers with this space as parent print using
        scientific notation.

        INPUT:
            status -- (bool --) optional flag
        """
        if status is None:
            return bool(self.sci_not)
        else:
            self.sci_not = status

    def zeta(self, n=2):
        """
        Return an $n$-th root of unity in the real field,
        if one exists, or raise a ValueError otherwise.

        EXAMPLES:
            sage: R = RealField()
            sage: R.zeta()
            -1.0000000000000000
            sage: R.zeta(1)
            1.0000000000000000
            sage: R.zeta(5)
            Traceback (most recent call last):
            ...
            ValueError: No 5th root of unity in self
        """
        if n == 1:
            return self(1)
        elif n == 2:
            return self(-1)
        raise ValueError, "No %sth root of unity in self"%n

R = RealField()

#*****************************************************************************
#
#     RealNumber -- element of Real Field
#
#
#
#*****************************************************************************
cdef class RealNumber(sage.structure.element.RingElement):
    """
    A real number.
    """
    def __new__(self, RealField parent, x=0, int base=10):
        self.init = 0

    def __init__(self, RealField parent, x=0, int base=10, special=None):
        """
        Create a real number.  Should be called by first creating
        a RealField, as illustrated in the examples.

        EXAMPLES:
            sage: R = RealField()
            sage: R('1.2456')
            1.2456000000000000
            sage: R = RealField(3)
            sage: R('1.2456')
            1.2

        EXAMPLE: Rounding Modes
            sage: w = RealField(3)(5/2)
            sage: RealField(2, rnd="RNDZ")(w).str(2)
            '10'
            sage: RealField(2, rnd="RNDD")(w).str(2)
            '10'
            sage: RealField(2, rnd="RNDU")(w).str(2)
            '11'
            sage: RealField(2, rnd="RNDN")(w).str(2)
            '10'

        NOTES: A real number is an arbitrary precision mantissa with a
        limited precision exponent.  A real number can have three
        special values: Not-a-Number (NaN) or plus or minus
        Infinity. NaN represents an uninitialized object, the result
        of an invalid operation (like 0 divided by 0), or a value that
        cannot be determined (like +Infinity minus
        +Infinity). Moreover, like in the IEEE 754-1985 standard, zero
        is signed, i.e. there are both +0 and ?0; the behavior is the
        same as in the IEEE 754-1985 standard and it is generalized to
        the other functions supported by MPFR.

        """
        if parent is None:
            raise TypeError
        self._parent = parent
        mpfr_init2(self.value, parent.__prec)
        self.init = 1
        if x is None: return
        cdef RealNumber _x, n, d
        if isinstance(x, RealNumber):
            _x = x  # so we can get at x.value
            mpfr_set(self.value, _x.value, parent.rnd)
        elif isinstance(x, sage.rings.rational.Rational):
            n = parent(x.numerator())
            d = parent(x.denominator())
            mpfr_div(self.value, n.value, d.value, parent.rnd)
        else:
            s = str(x).replace(' ','')
            if mpfr_set_str(self.value, s, base, parent.rnd):
                if s == 'NaN' or s == '@NaN@':
                    mpfr_set_nan(self.value)
                elif s == '+infinity':
                    mpfr_set_inf(self.value, 1)
                elif s == '-infinity':
                    mpfr_set_inf(self.value, -1)
                else:
                    raise TypeError, "Unable to convert x (='%s') to real number."%s

    cdef RealField c_parent(RealNumber self):
        """
        Returns the parent of self as a RealField.

        TODO: we *really* need to work out a way to make this insanely
        fast. This is called extremely frequently. All this overhead is
        probably already unacceptable.
        """
        return <RealField> self._parent


    def __reduce__(self):
        """
        EXAMPLES:
            sage: R = RealField(sci_not=1, prec=200, rnd='RNDU')
            sage: b = R('393.39203845902384098234098230948209384028340')
            sage: loads(dumps(b)) == b
            True
            sage: b = R(1)/R(0); b
            +infinity
            sage: loads(dumps(b)) == b
            True
            sage: b = R(-1)/R(0); b
            -infinity
            sage: loads(dumps(b)) == b
            True
            sage: b = R(-1).sqrt(); b
            1.0000000000000000000000000000000000000000000000000000000000000*I
            sage: loads(dumps(b)) == b
            True
        """
        s = self.str(32, no_sci=False, e='@')
        return (sage.rings.real_field.__reduce__RealNumber, (self.c_parent(), s, 32))

    def  __dealloc__(self):
        if self.init:
            mpfr_clear(self.value)

    def __repr__(self):
        return self.str(10)

    def _latex_(self):
        return str(self)

    def _interface_init_(self):
        """
        Return string representation of self in base 10 with
        no scientific notation.

        This is most likely to make sense in other computer algebra
        systems (this function is the default for exporting to other
        computer algebra systems).

        EXAMPLES:
            sage: n = 1.3939494594
            sage: n._interface_init_()
            '1.3939494593999999'
        """
        return self.str(10, no_sci=True)

    def __hash__(self):
        return hash(self.str(16))

    def _im_gens_(self, codomain, im_gens):
        return codomain(self) # since 1 |--> 1

    def real(self):
        """
        Return the real part of self.

        (Since self is a real number, this simply returns self.)
        """
        return self

    def parent(self):
        """
        EXAMPLES:
            sage: R = RealField()
            sage: a = R('1.2456')
            sage: a.parent()
            Real Field with 53 bits of precision
        """
        return self.c_parent()

    def str(self, int base=10, no_sci=None, e='e'):
        if base < 2 or base > 36:
            raise ValueError, "the base (=%s) must be between 2 and 36"%base
        if mpfr_nan_p(self.value):
            if base >= 24:
                return "@NaN@"
            else:
                return "NaN"
        elif mpfr_inf_p(self.value):
            if mpfr_sgn(self.value) > 0:
                return "+infinity"
            else:
                return "-infinity"

        cdef char *s
        cdef mp_exp_t exponent

        _sig_on
        s = mpfr_get_str(<char*>0, &exponent, base, 0,
                         self.value, self.c_parent().rnd)
        _sig_off
        if s == <char*> 0:
            raise RuntimeError, "Unable to convert an mpfr number to a string."
        t = str(s)
        free(s)

        if no_sci==False or (self.c_parent().sci_not and not (no_sci==True)):
            if t[0] == "-":
                return "-%s.%s%s%s"%(t[1:2], t[2:], e, exponent-1)
            return "%s.%s%s%s"%(t[0], t[1:], e, exponent-1)

        n = abs(exponent)
        lpad = ''
        if exponent <= 0:
            n = len(t)
            lpad = '0.' + '0'*abs(exponent)
        if t[0] == '-':
            lpad = '-' + lpad
            t = t[1:]
        z = lpad + str(t[:n])
        w = t[n:]
        if len(w) > 0:
            z = z + ".%s"%w
        elif lpad == '':
            z = z + '0'*(n-len(t))
        return z

    def copy(self):
        cdef RealNumber z
        z = RealNumber(self.c_parent())
        mpfr_set(z.value, self.value, self.c_parent().rnd)
        return z

    def integer_part(self):
        """
        If in decimal this number is written n.defg, returns n.
        The output is a SAGE integer.
        """
        s = self.str(base=32, no_sci=True)
        i = s.find(".")
        return sage.rings.integer.Integer(s[:i], base=32)

    ########################
    #   Basic Arithmetic
    ########################

    def _add_(RealNumber self, RealNumber other):
        """
        Add two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealField()
            sage: R(-1.5) + R(2.5)
            1.0000000000000000
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        mpfr_add(x.value, self.value, other.value, self.c_parent().rnd)
        return x

    def __invert__(self):
        return self.c_parent()(1) / self

    def _sub_(RealNumber self, RealNumber other):
        """
        Subtract two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealField()
            sage: R(-1.5) - R(2.5)
            -4.0000000000000000
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        mpfr_sub(x.value, self.value, other.value, self.c_parent().rnd)
        return x

    def _mul_(RealNumber self, RealNumber other):
        """
        Multiply two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealField()
            sage: R(-1.5) * R(2.5)
            -3.7500000000000000

        If two elements have different precision, arithmetic
        operations are performed after coercing to the lower
        precision.

            sage: R10 = RealField(10)
            sage: R100 = RealField(100)
            sage: a = R10('393.3902834028345')
            sage: b = R100('393.3902834028345')
            sage: a
            393.50
            sage: b
            393.39028340283450000000000000011
            sage: a*b
            154880
            sage: b*a
            154880
            sage: parent(b*a)
            Real Field with 10 bits of precision
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        mpfr_mul(x.value, self.value, other.value, self.c_parent().rnd)
        return x

    def __div_(RealNumber self, RealNumber other):
        if not other:
            raise ZeroDivisionError, "RealNumber division by zero"
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        mpfr_div(x.value, self.value, other.value, self.c_parent().rnd)
        return x

    def __div__(x, y):
        """
        EXAMPLES:
            sage: R = RealField()
            sage: R(-1.5) / R(2.5)
            -0.59999999999999998
        """
        if isinstance(x, RealNumber) and isinstance(y, RealNumber):
            return x.__div_(y)
        return sage.rings.coerce.bin_op(x, y, operator.div)

    def __neg__(self):
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        mpfr_neg(x.value, self.value, self.c_parent().rnd)
        return x

    def __pos__(self):
        return self

    def __abs__(self):
        return self.abs()

    cdef RealNumber abs(RealNumber self):
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        mpfr_abs(x.value, self.value, self.c_parent().rnd)
        return x

    # Bit shifting
    def _lshift_(RealNumber self, n):
        cdef RealNumber x
        if n > sys.maxint:
            raise OverflowError, "n (=%s) must be <= %s"%(n, sys.maxint)
        x = RealNumber(self.c_parent(), None)
        mpfr_mul_2exp(x.value, self.value, n, self.c_parent().rnd)
        return x

    def __lshift__(x, y):
        """
        EXAMPLES:
            sage: 1.0 << 32
            4294967296.0000000
        """
        if isinstance(x, RealNumber) and isinstance(y, (int,long,sage.rings.integer.Integer)):
            return x._lshift_(y)
        return sage.rings.coerce.bin_op(x, y, operator.lshift)

    def _rshift_(RealNumber self, n):
        if n > sys.maxint:
            raise OverflowError, "n (=%s) must be <= %s"%(n, sys.maxint)
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        mpfr_div_2exp(x.value, self.value, n, self.c_parent().rnd)
        return x

    def __rshift__(x, y):
        """
        EXAMPLES:
            sage: 1024.0 >> 7
            8.0000000000000000
        """
        if isinstance(x, RealNumber) and isinstance(y, (int,long,sage.rings.integer.Integer)):
            return x._rshift_(y)
        return sage.rings.coerce.bin_op(x, y, operator.rshift)

    def multiplicative_order(self):
        if self == 1:
            return 1
        elif self == -1:
            return -1
        return sage.rings.infinity.infinity

    def sign(self):
        return mpfr_sgn(self.value)

    def prec(self):
        return self.c_parent().__prec

    ###################
    # Rounding etc
    ###################

    def round(self):
        """
        Rounds self to the nearest real number. There are 4
        rounding modes. They are

        EXAMPLES:
            RNDN -- round to nearest:

            sage: R = RealField(10,False,'RNDN')
            sage: R(22.454)
            22.469
            sage: R(22.491)
            22.500

            RNDZ -- round towards zero:
            sage: R = RealField(10,False,'RNDZ')
            sage: R(22.454)
            22.437
            sage: R(22.491)
            22.468

            RNDU -- round towards plus infinity:
            sage: R10 = RealField(10,False,'RNDU')
            sage: R10(22.454)
            22.469
            sage: R10(22.491)
            22.500

            RNDU -- round towards minus infinity:
            sage: R10 = RealField(10,False,'RNDD')
            sage: R10(22.454)
            22.437
            sage: R10(22.491)
            22.468

        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        mpfr_round(x.value, self.value)
        return x

    def floor(self):
        """
        Returns the floor of this number

        EXAMPLES:
            sage: R = RealField()
            sage: (2.99).floor()
            2
            sage: (2.00).floor()
            2
            sage: floor(RR(-5/2))
            -3
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        mpfr_floor(x.value, self.value)
        return x.integer_part()

    def ceil(self):
        """
        Returns the ceiling of this number

        OUTPUT:
            integer

        EXAMPLES:
            sage: (2.99).ceil()
            3
            sage: (2.00).ceil()
            2
            sage: (2.01).ceil()
            3
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        mpfr_ceil(x.value, self.value)
        return x.integer_part()

    def ceiling(self):
        return self.ceil()

    def trunc(self):
        """
        Truncates this number

        EXAMPLES:
            sage: (2.99).trunc()
            2.0000000000000000
            sage: (-0.00).trunc()
            -0.00000000000000000
            sage: (0.00).trunc()
            0.00000000000000000
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        mpfr_trunc(x.value, self.value)
        return x

    def frac(self):
        """
        frac returns a real number > -1 and < 1. it satisfies the
        relation:
            x = x.trunc() + x.frac()

        EXAMPLES:
            sage: (2.99).frac()
            0.99000000000000021
            sage: (2.50).frac()
            0.50000000000000000
            sage: (-2.79).frac()
            -0.79000000000000004
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        mpfr_frac(x.value, self.value, self.c_parent().rnd)
        return x

    ###########################################
    # Conversions
    ###########################################

    def __float__(self):
        return mpfr_get_d(self.value, self.c_parent().rnd)

    def __int__(self):
        """
        Returns integer truncation of this real number.
        """
        s = self.str(32)
        i = s.find('.')
        return int(s[:i], 32)

    def __long__(self):
        """
        Returns long integer truncation of this real number.
        """
        s = self.str(32)
        i = s.find('.')
        return long(s[:i], 32)

    def __complex__(self):
        return complex(float(self))

    def _complex_number_(self):
        return sage.rings.complex_field.ComplexField(self.prec())(self)

    def _pari_(self):
        return sage.libs.pari.all.pari.new_with_bits_prec(str(self), self.c_parent().__prec)


    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    def is_NaN(self):
        return bool(mpfr_nan_p(self.value))

    cdef int cmp(RealNumber self, RealNumber x):
        cdef int a,b
        a = mpfr_nan_p(self.value)
        b = mpfr_nan_p(x.value)
        if a != b:
            return -1    # nothing is equal to Nan
        cdef int i
        i = mpfr_cmp(self.value, x.value)
        if i < 0:
            return -1
        elif i == 0:
            return 0
        else:
            return 1

    def __cmp__(RealNumber self, RealNumber x):
        return self.cmp(x)

    def __richcmp__(RealNumber self, x, int op):
        cdef int n
        if not isinstance(x, RealNumber):
            try:
                x = RealNumber(self.c_parent(), x)
            except TypeError:
                n = sage.rings.coerce.cmp(self, x)
            else:
                n = self.cmp(x)
        else:
            n = self.cmp(x)
        if op == 0:
            return bool(n < 0)
        elif op == 1:
            return bool(n <= 0)
        elif op == 2:
            return bool(n == 0)
        elif op == 3:
            return bool(n != 0)
        elif op == 4:
            return bool(n > 0)
        elif op == 5:
            return bool(n >= 0)


    ############################
    # Special Functions
    ############################

    def sqrt(self):
        """
        Return a square root of self.

        If self is negative a complex number is returned.

        If you use self.square_root() then a real number will always
        be returned (though it will be NaN if self is negative).

        EXAMPLES:
            sage: r = 4.0
            sage: r.sqrt()
            2.0000000000000000
            sage: r.sqrt()^2 == r
            True

            sage: r = 4344
            sage: r.sqrt()
            65.909028213136324
            sage: r.sqrt()^2 == r
            True

            sage: r = -2.0
            sage: r.sqrt()
            1.4142135623730951*I
            """
        if self >= 0:
            return self.square_root()
        return self._complex_number_().sqrt()


    def square_root(self):
        """
        Return a square root of self.  A real number will always be
        returned (though it will be NaN if self is negative).

        Use self.sqrt() to get a complex number if self is negative.

        EXAMPLES:
            sage: r = -2.0
            sage: r.square_root()
            NaN
            sage: r.sqrt()
            1.4142135623730951*I
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_sqrt(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def cube_root(self):
        """
        Return the cubic root (defined over the real numbers) of self.

        EXAMPLES:
            sage: r = 125.0; r.cube_root()
            5.0000000000000000
            sage: r = -119.0
            sage: r.cube_root()^3 - r       # illustrates precision loss
            -0.000000000000014210854715202004
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_cbrt(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def __pow(self, RealNumber exponent):
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_pow(x.value, self.value, exponent.value, self.c_parent().rnd)
        _sig_off
        if mpfr_nan_p(x.value):
            return self._complex_number_()**exponent._complex_number_()
        return x

    def __pow__(self, exponent, modulus):
        """
        Compute self raised to the power of exponent, rounded in
        the direction specified by the parent of self.

        If the result is not a real number, self and the exponent are
        both coerced to complex numbers (with sufficient precision),
        then the exponentiation is computed in the complex numbers.
        Thus this function can return either a real or complex number.

        EXAMPLES:
            sage: R = RealField(30)
            sage: a = R('1.23456')
            sage: a^20
            67.646297455
            sage: a^a
            1.2971114814
            sage: b = R(-1)
            sage: b^(1/2)
            -0.00000000000000000010842021725 + 1.0000000000*I
        """
        cdef RealNumber x
        if not isinstance(self, RealNumber):
            return self.__pow__(float(exponent))
        if not isinstance(exponent, RealNumber):
            x = self
            exponent = x.c_parent()(exponent)
        return self.__pow(exponent)

    def log(self, base='e'):
        """
        EXAMPLES:
            sage: R = RealField()
            sage: R(2).log()
            0.69314718055994529
        """
        cdef RealNumber x
        if base == 'e':
            x = RealNumber(self.c_parent(), None)
            _sig_on
            mpfr_log(x.value, self.value, self.c_parent().rnd)
            _sig_off
            return x
        elif base == 10:
            return self.log10()
        elif base == 2:
            return self.log2()
        else:
            return self.log() / (self.parent()(base)).log()

    def log2(self):
        """
        Returns log to the base 2 of self

        EXAMPLES:
            sage: r = 16.0
            sage: r.log2()
            4.0000000000000000

            sage: r = 31.9; r.log2()
            4.9954845188775066

            sage: r = 0.0
            sage: r.log2()
            -infinity
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_log2(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def log10(self):
        """
        Returns log to the base 10 of self

        EXAMPLES:
            sage: r = 16.0; r.log10()
            1.2041199826559248
            sage: r.log() / log(10)
            1.2041199826559246

            sage: r = 39.9; r.log10()
            1.6009728956867482

            sage: r = 0.0
            sage: r.log10()
            -infinity

            sage: r = -1.0
            sage: r.log10()
            NaN

        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_log10(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def exp(self):
        r"""
        Returns $e^\code{self}$

        EXAMPLES:
            sage: r = 0.0
            sage: r.exp()
            1.0000000000000000

            sage: r = 32.3
            sage: a = r.exp(); a
            106588847274864.47
            sage: a.log()
            32.299999999999997

            sage: r = -32.3
            sage: r.exp()
            0.0000000000000093818445884986851
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_exp(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def exp2(self):
        """
        Returns $2^\code{self}$

        EXAMPLES:
            sage: r = 0.0
            sage: r.exp2()
            1.0000000000000000

            sage: r = 32.0
            sage: r.exp2()
            4294967296.0000000

            sage: r = -32.3
            sage: r.exp2()
            0.00000000018911724825302072

        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_exp2(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def exp10(self):
        r"""
        Returns $10^\code{self}$

        EXAMPLES:
            sage: r = 0.0
            sage: r.exp10()
            1.0000000000000000

            sage: r = 32.0
            sage: r.exp10()
            100000000000000010000000000000000

            sage: r = -32.3
            sage: r.exp10()
            0.0000000000000000000000000000000050118723362727556
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_exp10(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def cos(self):
        """
        Returns the cosine of this number

        EXAMPLES:
            sage: t=RR.pi()/2
            sage: t.cos()
            0.000000000000000061232339957367660
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_cos(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    ##########################################################
    # it would be nice to get zero back here:
    # sage: R(-1).acos().sin()
    # _57 = -0.50165576126683320234e-19
    # i think this could be "fixed" by using MPFI. (put on to-do list.)
    #
    # this seems to work ok:
    # sage: R(-1).acos().cos()
    # _58 = -0.10000000000000000000e1
    def sin(self):
        """
        Returns the sine of this number

        EXAMPLES:
            sage: R = RealField(100)
            sage: R(2).sin()
            0.90929742682568169539601986591150
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_sin(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def tan(self):
        """
        Returns the tangent of this number

        EXAMPLES:
            sage: q = RR.pi()/3
            sage: q.tan()
            1.7320508075688767
            sage: q = RR.pi()/6
            sage: q.tan()
            0.57735026918962573
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_tan(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def sincos(self):
        """
        Returns a pair consisting of the sine and cosine.

        EXAMPLES:
            sage: R = RealField()
            sage: t = R.pi()/6
            sage: t.sincos()
            (0.49999999999999994, 0.86602540378443871)
        """
        cdef RealNumber x,y
        x = RealNumber(self.c_parent(), None)
        y = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_sin_cos(x.value, y.value, self.value, self.c_parent().rnd)
        _sig_off
        return x,y


    # int mpfr_sin_cos (mpfr_t rop, mpfr_t op, mpfr_t, mp_rnd_t rnd)

    def acos(self):
        """
        Returns the inverse cosine of this number

        EXAMPLES:
            sage: q = RR.pi()/3
            sage: i = q.cos()
            sage: i.acos() == q
            True
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_acos(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def asin(self):
        """
        Returns the inverse sine of this number

        EXAMPLES:
            sage: q = RR.pi()/5
            sage: i = q.sin()
            sage: i.asin() == q
            True
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_asin(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def atan(self):
        """
        Returns the inverse tangent of this number

        EXAMPLES:
            sage: q = RR.pi()/5
            sage: i = q.tan()
            sage: i.atan() == q
            True
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_atan(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    #int mpfr_acos _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
    #int mpfr_asin _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));
    #int mpfr_atan _PROTO ((mpfr_ptr, mpfr_srcptr, mp_rnd_t));

    def cosh(self):
        """
        Returns the hyperbolic cosine of this number

        EXAMPLES:
            sage: q = RR.pi()/12
            sage: q.cosh()
            1.0344656400955106
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_cosh(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def sinh(self):
        """
        Returns the hyperbolic sine of this number

        EXAMPLES:
            sage: q = RR.pi()/12
            sage: q.sinh()
            0.26480022760227073

        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_sinh(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def tanh(self):
        """
        Returns the hyperbolic tangent of this number

        EXAMPLES:
            sage: q = RR.pi()/11
            sage: q.tanh()
            0.27807942929585028
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_tanh(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def acosh(self):
        """
        Returns the hyperbolic inverse cosine of this number

        EXAMPLES:
            sage: q = RR.pi()/2
            sage: i = q.cosh() ; i
            2.5091784786580567
            sage: i.acosh() == q
            True
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_acosh(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def asinh(self):
        """
        Returns the hyperbolic inverse sine of this number

        EXAMPLES:
            sage: q = RR.pi()/7
            sage: i = q.sinh() ; i
            0.46401763049299094
            sage: i.asinh() == q
            True
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_asinh(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def atanh(self):
        """
        Returns the hyperbolic inverse tangent of this number

        EXAMPLES:
            sage: q = RR.pi()/7
            sage: i = q.tanh() ; i
            0.42091124104853489
            sage: i.atanh() == q
            True
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_atanh(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def agm(self, other):
        """
        Return the arithmetic-geometric mean of self and other. The
        arithmetic-geometric mean is the common limit of the sequences
        $u_n$ and $v_n$, where $u_0$ is self, $v_0$ is other,
        $u_{n+1}$ is the arithmetic mean of $u_n$ and $v_n$, and
        $v_{n+1}$ is the geometric mean of u_n and v_n. If any operand
        is negative, the return value is \code{NaN}.
        """
        cdef RealNumber x, _other
        if not isinstance(other, RealNumber) or other.parent() != self.c_parent():
            _other = self.c_parent()(other)
        else:
            _other = other
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_agm(x.value, self.value, _other.value, self.c_parent().rnd)
        _sig_off
        return x


    def erf(self):
        """
        Returns the value of the error function on self.

        EXAMPLES:
           sage: R = RealField()
           sage: R(6).erf()
           1.0000000000000000
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_erf(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x


    def gamma(self):
        """
        The Euler gamma function. Return gamma of self.

        EXAMPLES:
           sage: R = RealField()
           sage: R(6).gamma()
           120.00000000000000
           sage: R(1.5).gamma()
           0.88622692545275805
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_gamma(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def zeta(self):
        r"""
        Return the Riemann zeta function evaluated at this real number.

        \note{PARI is vastly more efficient at computing the Riemann zeta
        function.   See the example below for how to use it.}

        EXAMPLES:
            sage: R = RealField()
            sage: R(2).zeta()
            1.6449340668482264
            sage: R.pi()^2/6
            1.6449340668482264
            sage: R(-2).zeta()
            0.00000000000000000
            sage: R(1).zeta()
            +infinity

        Computing zeta using PARI is much more efficient in difficult cases.
        Here's how to compute zeta with at least a given precision:

             sage: z = pari.new_with_bits_prec(2, 53).zeta(); z
             1.644934066848226436472415167              # 32-bit
             1.6449340668482264364724151666460251892    # 64-bit

        Note that the number of bits of precision in the constructor only
        effects the internel precision of the pari number, not the number
        of digits that gets displayed.  To increase that you must
        use \code{pari.set_real_precision}.

             sage: type(z)
             <type 'gen.gen'>
             sage: R(z)
             1.6449340668482264
        """
        cdef RealNumber x
        x = RealNumber(self.c_parent(), None)
        _sig_on
        mpfr_zeta(x.value, self.value, self.c_parent().rnd)
        _sig_off
        return x

    def algdep(self, n):
        """
        Returns a polynomial of degree at most $n$ which is approximately
        satisfied by this number.  Note that the returned polynomial
        need not be irreducible, and indeed usually won't be if this number
        is a good approximation to an algebraic number of degree less than $n$.

        ALGORITHM: Uses the PARI C-library algdep command.

        EXAMPLE:
             sage: r = sqrt(2); r
             1.4142135623730951
             sage: r.algdep(5)
             x^2 - 2                        # 32-bit
             x^5 - x^4 - 2*x^3 + x^2 + 2    # 64-bit
        """
        return sage.rings.arith.algdep(self,n)

    def algebraic_dependency(self, n):
        """
         Returns a polynomial of degree at most $n$ which is approximately
         satisfied by this number.  Note that the returned polynomial
         need not be irreducible, and indeed usually won't be if this number
         is a good approximation to an algebraic number of degree less than $n$.

         ALGORITHM: Uses the PARI C-library algdep command.

         EXAMPLE:
              sage: r = sqrt(2); r
              1.4142135623730951
              sage: r.algdep(5)
              x^2 - 2                              # 32-bit
              x^5 - x^4 - 2*x^3 + x^2 + 2          # 64-bit
        """
        return sage.rings.arith.algdep(self,n)

RR = RealField()


def create_RealNumber(s, int base=10, int pad=0, rnd="RNDN", min_prec=53):
    r"""
    Return the real number defined by the string s as an element
    of \code{RealField(prec=n)}, where n has slightly more (controlled
    by pad) bits than given by s.

    INPUT:
        s -- a string that defines a real number
        base -- an integer between 2 and 36
        pad -- an integer >= 1.
        rnd -- rounding mode: RNDN, RNDZ, RNDU, RNDD
        min_prec -- number will have at least this many bits of precision, no matter what.
    """
    if base == 10:
        bits = int(3.32192*len(s))
    else:
        bits = int(math.log(base,2)*len(s))
    R = RealField(prec=max(bits+pad, min_prec), rnd=rnd)
    return RealNumber(R, s, base)

