"""
Field of Arbitrary Precision Real Numbers

AUTHORS:
   -- Kyle Schalm <kschalm@math.utexas.edu> (2005-09)
   -- William Stein <wstein@gmail.com>: bug fixes, examples, maintenance
   -- Didier Deshommes <dfdeshom@gmail.com> (2006-03-19): examples
   -- David Harvey (2006-09-20): compatibility with Element._parent
   -- William Stein (2006-10): default printing truncates to avoid base-2
              rounding confusing (fix suggested by Bill Hart)

EXAMPLES:

A difficult conversion:

    sage: RR(sys.maxint)
    9223372036854780000     # 64-bit
    2147483647.00000        # 32-bit
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

import math # for log
import sys

include '../ext/interrupt.pxi'
include "../ext/stdsage.pxi"

cimport sage.rings.ring
import  sage.rings.ring

cimport sage.structure.element
from sage.structure.element cimport RingElement, Element, ModuleElement
import  sage.structure.element

import sage.misc.misc as misc

import sage.structure.coerce
import operator

from sage.libs.pari.gen import PariInstance, gen
from sage.libs.pari.gen cimport PariInstance, gen

from integer import Integer
from integer cimport Integer
from rational import Rational
from rational cimport Rational

from real_double cimport RealDoubleElement
from real_double import is_RealDoubleElement

import sage.rings.complex_field

import sage.rings.infinity

from sage.structure.parent_gens cimport ParentWithGens

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

        sci_not -- (default: False) if True, always display
                using scientific notation; if False, display
                using scientific notation only for very large or
                very small numbers

        rnd -- (string) the rounding mode
                RNDN -- (default) round to nearest: Knuth says this is
                        the best choice to prevent ``floating point
                        drift''.
                RNDD -- round towards minus infinity
                RNDZ -- round towards zero
                RNDU -- round towards plus infinity

    EXAMPLES:
        sage: RealField(10)
        Real Field with 10 bits of precision
        sage: RealField()
        Real Field with 53 bits of precision
        sage: RealField(100000)
        Real Field with 100000 bits of precision

    NOTE: The default precision is 53, since according to the MPFR
       manual: 'mpfr should be able to exactly reproduce all
       computations with double-precision machine floating-point
       numbers (double type in C), except the default exponent
       range is much wider and subnormal numbers are not
       implemented.'
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
        ParentWithGens.__init__(self, self, tuple([]), False)

    cdef RealNumber _new(self):
        """
        Return a new real number with parent self.
        """
        cdef RealNumber x
        x = PY_NEW(RealNumber)
        x._parent = self
        mpfr_init2(x.value, self.__prec)
        x.init = 1
        return x

    def _repr_(self):
        s = "Real Field with %s bits of precision"%self.__prec
        if self.rnd != GMP_RNDN:
            s = s + " and rounding %s"%(self.rnd_str)
        return s

    def _latex_(self):
        return "\\R"

    def is_exact(self):
        return False

    def __call__(self, x, base=10):
        """
        Coerce x into this real field.

        EXAMPLES:
            sage: R = RealField(20)
            sage: R('1.234')
            1.2340
            sage: R('2', base=2)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert x (='2') to real number.
            sage: a = R('1.1001', base=2); a
            1.5625
            sage: a.str(2)
            '1.1001000000000000000'
        """
        if hasattr(x, '_mpfr_'):
            # This design with the hasattr is very annoying.
            # The only thing that uses it right now is symbolic constants
            # and symbolic function evaluation.
            # Getting rid of this would speed things up.
            return x._mpfr_(self)
        cdef RealNumber z
        z = self._new()
        z._set(x, base)
        return z

    cdef _coerce_c_impl(self, x):
        """
        Canonical coercion of x to this mpfr real field.

        The rings that canonically coerce to this mpfr real field are:
             * this real field itself
             * any other mpfr real field with precision that is as large as this one
             * int, long, integer, and rational rings.
             * real mathematical constants
        """
        if isinstance(x, RealNumber):
            P = x.parent()
            if (<RealField> P).__prec >= self.__prec:
                return self(x)
            else:
                raise TypeError, "Canonical coercion from lower to higher precision not defined"
        elif isinstance(x, (Integer, Rational)):
            return self(x)
        elif self.__prec <= 53 and is_RealDoubleElement(x):
            return self(x)
        import sage.functions.constants
        return self._coerce_try(x, [sage.functions.constants.ConstantRing])

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
        return __create__RealField_version0, (self.__prec, self.sci_not, self.rnd_str)

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
            3.1415926535897932384626433833
            sage: R.pi().sqrt()/2
            0.88622692545275801364908374167
            sage: R = RealField(150)
            sage: R.pi().sqrt()/2
            0.88622692545275801364908374167057259139877473
        """
        cdef RealNumber x
        x = self._new()
        mpfr_const_pi(x.value, self.rnd)
        return x


    # int mpfr_const_euler (mpfr_t rop, mp_rnd_t rnd)
    def euler_constant(self):
        """
        Returns Euler's gamma constant to the precision of this field.

        EXAMPLES:
            sage: RealField(100).euler_constant()
            0.57721566490153286060651209008
        """
        cdef RealNumber x
        x = self._new()
        mpfr_const_euler(x.value, self.rnd)
        return x

    # int mpfr_const_catalan (mpfr_t rop, mp_rnd_t rnd)
    def catalan_constant(self):
        """
        Returns Catalan's constant to the precision of this field.

        EXAMPLES:
            sage: RealField(100).catalan_constant()
            0.91596559417721901505460351493
        """
        cdef RealNumber x
        x = self._new()
        mpfr_const_catalan(x.value, self.rnd)
        return x

    # int mpfr_const_log2 (mpfr_t rop, mp_rnd_t rnd)
    def log2(self):
        """
        Returns log(2) to the precision of this field.

        EXAMPLES:
            sage: R=RealField(100)
            sage: R.log2()
            0.69314718055994530941723212146
            sage: R(2).log()
            0.69314718055994530941723212146
        """
        cdef RealNumber x
        x = self._new()
        mpfr_const_log2(x.value, self.rnd)
        return x

    def factorial(self, int n):
        """
        Return the factorial of the integer n as a real number.
        """
        cdef RealNumber x
        if n < 0:
            raise ArithmeticError, "n must be nonnegative"
        x = self._new()
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
            -1.00000000000000
            sage: R.zeta(1)
            1.00000000000000
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

    Real numbers are printed to slightly less digits than their
    internal precision, in order to avoid confusing roundoff issues
    that occur because numbers are stored internally in binary.
    """
    cdef RealNumber _new(self):
        """
        Return a new real number with same parent as self.
        """
        cdef RealNumber x
        x = PY_NEW(RealNumber)
        x._parent = self._parent
        mpfr_init2(x.value, (<RealField>self._parent).__prec)
        x.init = 1
        return x

    def __init__(self, RealField parent, x=0, int base=10):
        """
        Create a real number.  Should be called by first creating
        a RealField, as illustrated in the examples.

        EXAMPLES:
            sage: R = RealField()
            sage: R('1.2456')
            1.24560000000000
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
        is signed, i.e. there are both +0 and -0; the behavior is the
        same as in the IEEE 754-1985 standard and it is generalized to
        the other functions supported by MPFR.

        """
        self.init = 0
        if parent is None:
            raise TypeError
        self._parent = parent
        mpfr_init2(self.value, parent.__prec)
        self.init = 1
        if x is None: return
        self._set(x, base)

    cdef _set(self, x, int base):
        # This should not be called except when the number is being created.
        # Real Numbers are supposed to be immutable.
        cdef RealNumber _x, n, d
        cdef Integer _ix
        cdef RealField parent
        cdef gen _gen
        cdef RealDoubleElement rd
        parent = self._parent
        if PY_TYPE_CHECK(x, RealNumber):
            _x = x  # so we can get at x.value
            mpfr_set(self.value, _x.value, parent.rnd)
        elif PY_TYPE_CHECK(x, Integer):
            mpfr_set_z(self.value, (<Integer>x).value, parent.rnd)
        elif PY_TYPE_CHECK(x, Rational):
            mpfr_set_q(self.value, (<Rational>x).value, parent.rnd)
        elif PY_TYPE_CHECK(x, gen) and x.type() == "t_REAL":
            _gen = x
            self._set_from_GEN_REAL(_gen.g)
        elif isinstance(x, (int, long)):
            _ix = Integer(x)
            mpfr_set_z(self.value, _ix.value, parent.rnd)
        elif isinstance(x, float):
            mpfr_set_d(self.value, x, parent.rnd)
        elif PY_TYPE_CHECK(x, RealDoubleElement):
            rd = x
            mpfr_set_d(self.value, rd._value, parent.rnd)
        #elif hasattr(x, '_mpfr_'):
        #    return x._mpfr_(self)
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

    cdef _set_from_GEN_REAL(self, GEN g):
        """
        EXAMPLES:
            sage: rt2 = sqrt(pari('2.0'))
            sage: rt2
            1.414213562373095048801688724              # 32-bit
            1.4142135623730950488016887242096980786    # 64-bit
            sage: rt2.python()
            1.414213562373095048801688724              # 32-bit
            1.4142135623730950488016887242096980786    # 64-bit
            sage: rt2.python().prec()
            96                                         # 32-bit
            128                                        # 64-bit
            sage: pari(rt2.python()) == rt2
            True
            sage: for i in xrange(1, 1000):
            ...       assert(sqrt(pari(i)) == pari(sqrt(pari(i)).python()))
            sage: (-3.1415)._pari_().python()
            -3.14150000000000018
        """
        cdef int sgn
        sgn = signe(g)

        if sgn == 0:
            mpfr_set_ui(self.value, 0, GMP_RNDN)
            return

        cdef int wordsize

        if sage.misc.misc.is_64_bit:
            wordsize = 64
        else:
            wordsize = 32

        cdef mpz_t mantissa
        mpz_init(mantissa)

        mpz_import(mantissa, lg(g) - 2, 1, wordsize/8, 0, 0, &g[2])

        cdef int exponent
        exponent = expo(g)

        # Round to nearest for best results when setting a low-precision
        # MPFR from a high-precision GEN
        mpfr_set_z(self.value, mantissa, GMP_RNDN)
        mpfr_mul_2si(self.value, self.value, exponent - wordsize * (lg(g) - 2) + 1, GMP_RNDN)

        if sgn < 0:
            mpfr_neg(self.value, self.value, GMP_RNDN)

        mpz_clear(mantissa)

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
            1.0000000000000000000000000000000000000000000000000000000000*I
            sage: loads(dumps(b)) == b
            True
        """
        s = self.str(32, no_sci=False, e='@')
        return (__create__RealNumber_version0, (self._parent, s, 32))

    def  __dealloc__(self):
        if self.init:
            mpfr_clear(self.value)

    def __repr__(self):
        return self.str(10)

    def _latex_(self):
        return str(self)

    def _interface_init_(self):
        """
        Return string representation of self in base 10, avoiding
        scientific notation except for very large or very small numbers.

        This is most likely to make sense in other computer algebra
        systems (this function is the default for exporting to other
        computer algebra systems).

        EXAMPLES:
            sage: n = 1.3939494594
            sage: n._interface_init_()
            '1.39394945940000'
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
        return self._parent

    def str(self, int base=10, no_sci=None, e=None, int truncate=1):
        """
        INPUT:
             base -- base for output
             no_sci -- if 2, never print using scientific notation;
                       if 1 or True, print using scientific notation only
                       for very large or very small numbers;
                       if 0 or False always print with scientific notation;
                       if None (the default), print how the parent prints.
             e - symbol used in scientific notation; defaults to 'e' for
                       base<=10, and '@' otherwise
             truncate -- if True, round off the last digits in printing to
                       lessen confusing base-2 roundoff issues.

        EXAMPLES:
            sage: a = 61/3.0; a
            20.3333333333333
            sage: a.str(truncate=False)
            '20.333333333333332'
            sage: a.str(2)
            '10100.010101010101010101010101010101010101010101010101'
            sage: a.str(no_sci=False)
            '2.03333333333333e1'
            sage: a.str(16, no_sci=False)
            '1.4555555555555@1'
            sage: b = 2.0^99
            sage: b.str()
            '633825300114115000000000000000'
            sage: b.str(no_sci=False)
            '6.33825300114115e29'
            sage: b.str(no_sci=True)
            '633825300114115000000000000000'
            sage: c = 2.0^100
            sage: c.str()
            '1.26765060022823e30'
            sage: c.str(no_sci=False)
            '1.26765060022823e30'
            sage: c.str(no_sci=True)
            '1.26765060022823e30'
            sage: c.str(no_sci=2)
            '1267650600228230000000000000000'
            sage: 0.5^53
            0.000000000000000111022302462516
            sage: 0.5^54
            5.55111512312578e-17
        """
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

        if e is None:
            if base > 10:
                e = '@'
            else:
                e = 'e'

        cdef char *s
        cdef mp_exp_t exponent

        cdef int reqdigits

        reqdigits = 0

        if base == 10 and truncate:

            # This computes reqdigits == floor(log_{10}(2^(b-1))),
            # which is the number of *decimal* digits that are
            # "right", given that the last binary bit of the binary
            # number can be off.  That is, if this real is within a
            # relative error of 2^(-b) of an exact decimal with
            # reqdigits digits, that decimal will be returned.
            # This is equivalent to saying that exact decimals with
            # reqdigits digits differ by at least 2*2^(-b) (relative).

            # (Depending on the precision and the exact number involved,
            # adjacent exact decimals can differ by far more than 2*2^(-b)
            # (relative).)

            # This avoids the confusion a lot of people have with the last
            # 1-2 binary digits being wrong due to rounding coming from
            # representating numbers in binary.

            reqdigits = ((<RealField>self._parent).__prec - 1) * 0.3010299956
            if reqdigits <= 1: reqdigits = 2

        _sig_on
        s = mpfr_get_str(<char*>0, &exponent, base, reqdigits,
                         self.value, (<RealField>self._parent).rnd)
        _sig_off
        if s == <char*> 0:
            raise RuntimeError, "Unable to convert an mpfr number to a string."
        t = str(s)
        mpfr_free_str(s)


        cdef int digits
        digits = len(t)
        if t[0] == "-":
            digits = digits - 1

        if no_sci is None:
            no_sci = not (<RealField>self._parent).sci_not

        if no_sci==True and (-exponent > digits or exponent > 2*digits):
            no_sci = False

        if no_sci==False:
            if t[0] == "-":
                return "-%s.%s%s%s"%(t[1:2], t[2:], e, exponent-1)
            return "%s.%s%s%s"%(t[0], t[1:], e, exponent-1)

        lpad = ''
        if exponent <= 0:
            n = len(t)
            lpad = '0.' + '0'*abs(exponent)
        else:
            n = exponent
        if t[0] == '-':
            lpad = '-' + lpad
            t = t[1:]
        z = lpad + str(t[:n])
        w = t[n:]
        if len(w) > 0:
            z = z + ".%s"%w
        elif exponent > 0:
            z = z + '0'*(n-len(t))
        return z

    def __copy__(self):
        """
        Return copy of self -- since self is immutable, we just return self again.

        EXAMPLES:
            sage: a = 3.5
            sage: copy(a) is  a
            True
        """
        return self    # since object is immutable.

    def integer_part(self):
        """
        If in decimal this number is written n.defg, returns n.

        OUTPUT:
            -- a SAGE Integer

        EXAMPLE:
            sage: a = 119.41212
            sage: a.integer_part()
            119

        A big number with no decimal point:
            sage: a = RR(10^17); a
            100000000000000000
            sage: a.integer_part()
            100000000000000000
        """
        s = self.str(base=32, no_sci=True)
        i = s.find(".")
        if i != -1:
            return Integer(s[:i], base=32)
        else:
            return Integer(s, base=32)

    ########################
    #   Basic Arithmetic
    ########################

    cdef ModuleElement _add_c_impl(self, ModuleElement other):
        """
        Add two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealField()
            sage: R(-1.5) + R(2.5)
            1.00000000000000
        """
        cdef RealNumber x
        x = self._new()
        mpfr_add(x.value, self.value, (<RealNumber>other).value, (<RealField>self._parent).rnd)
        return x

    def __invert__(self):
        return self._parent(1) / self

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        Subtract two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealField()
            sage: R(-1.5) - R(2.5)
            -4.00000000000000
        """
        cdef RealNumber x
        x = self._new()
        mpfr_sub(x.value, self.value, (<RealNumber>right).value, (<RealField> self._parent).rnd)
        return x

    cdef RingElement _mul_c_impl(self, RingElement right):
        """
        Multiply two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealField()
            sage: R(-1.5) * R(2.5)
            -3.75000000000000

        If two elements have different precision, arithmetic
        operations are performed after coercing to the lower
        precision.

            sage: R20 = RealField(20)
            sage: R100 = RealField(100)
            sage: a = R20('393.3902834028345')
            sage: b = R100('393.3902834028345')
            sage: a
            393.39
            sage: b
            393.39028340283450000000000000
            sage: a*b
            154760
            sage: b*a
            154760
            sage: parent(b*a)
            Real Field with 20 bits of precision
        """
        cdef RealNumber x
        x = self._new()
        mpfr_mul(x.value, self.value, (<RealNumber>right).value, (<RealField>self._parent).rnd)
        return x


    cdef RingElement _div_c_impl(self, RingElement right):
        """
        Divide self by other, where both are real numbers with the same parent.

        EXAMPLES:
            sage: RR(1)/RR(3)
            0.333333333333333
            sage: RR(1)/RR(0)
            +infinity

            sage: R = RealField()
            sage: R(-1.5) / R(2.5)
            -0.600000000000000
        """
        cdef RealNumber x
        x = self._new()
        mpfr_div((<RealNumber>x).value, self.value,
                 (<RealNumber>right).value, (<RealField>self._parent).rnd)
        return x

    cdef ModuleElement _neg_c_impl(self):
        cdef RealNumber x
        x = self._new()
        mpfr_neg(x.value, self.value, (<RealField>self._parent).rnd)
        return x

    def __abs__(self):
        return self.abs()

    cdef RealNumber abs(RealNumber self):
        cdef RealNumber x
        x = self._new()
        mpfr_abs(x.value, self.value, (<RealField>self._parent).rnd)
        return x

    # Bit shifting
    def _lshift_(RealNumber self, n):
        cdef RealNumber x
        if n > sys.maxint:
            raise OverflowError, "n (=%s) must be <= %s"%(n, sys.maxint)
        x = self._new()
        mpfr_mul_2exp(x.value, self.value, n, (<RealField>self._parent).rnd)
        return x

    def __lshift__(x, y):
        """
        EXAMPLES:
            sage: 1.0 << 32
            4294967296.00000
        """
        if isinstance(x, RealNumber) and isinstance(y, (int,long, Integer)):
            return x._lshift_(y)
        return sage.structure.coerce.bin_op(x, y, operator.lshift)

    def _rshift_(RealNumber self, n):
        if n > sys.maxint:
            raise OverflowError, "n (=%s) must be <= %s"%(n, sys.maxint)
        cdef RealNumber x
        x = self._new()
        mpfr_div_2exp(x.value, self.value, n, (<RealField>self._parent).rnd)
        return x

    def __rshift__(x, y):
        """
        EXAMPLES:
            sage: 1024.0 >> 7
            8.00000000000000
        """
        if isinstance(x, RealNumber) and isinstance(y, (int,long,Integer)):
            return x._rshift_(y)
        return sage.structure.coerce.bin_op(x, y, operator.rshift)

    def multiplicative_order(self):
        if self == 1:
            return 1
        elif self == -1:
            return -1
        return sage.rings.infinity.infinity

    def sign(self):
        return mpfr_sgn(self.value)

    def prec(self):
        return (<RealField>self._parent).__prec

    ###################
    # Rounding etc
    ###################

    def round(self):
        """
        Rounds self to the nearest real number. There are 4
        rounding modes. They are

        EXAMPLES:
            RNDN -- round to nearest:

            sage: R = RealField(20,False,'RNDN')
            sage: R(22.454)
            22.454
            sage: R(22.491)
            22.491

            RNDZ -- round towards zero:
            sage: R = RealField(20,False,'RNDZ')
            sage: R(22.454)
            22.453
            sage: R(22.491)
            22.490

            RNDU -- round towards plus infinity:
            sage: R = RealField(20,False,'RNDU')
            sage: R(22.454)
            22.455
            sage: R(22.491)
            22.492

            RNDU -- round towards minus infinity:
            sage: R = RealField(20,False,'RNDD')
            sage: R(22.454)
            22.453
            sage: R(22.491)
            22.490
        """
        cdef RealNumber x
        x = self._new()
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
        x = self._new()
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

            sage: ceil(10^16 * 1.0)
            10000000000000000
            sage: ceil(10^17 * 1.0)
            100000000000000000
        """
        cdef RealNumber x
        x = self._new()
        mpfr_ceil(x.value, self.value)
        return x.integer_part()

    def ceiling(self):
        return self.ceil()

    def trunc(self):
        """
        Truncates this number

        EXAMPLES:
            sage: (2.99).trunc()
            2.00000000000000
            sage: (-0.00).trunc()
            -0.000000000000000
            sage: (0.00).trunc()
            0.000000000000000
        """
        cdef RealNumber x
        x = self._new()
        mpfr_trunc(x.value, self.value)
        return x

    def frac(self):
        """
        frac returns a real number > -1 and < 1. it satisfies the
        relation:
            x = x.trunc() + x.frac()

        EXAMPLES:
            sage: (2.99).frac()
            0.990000000000000
            sage: (2.50).frac()
            0.500000000000000
            sage: (-2.79).frac()
            -0.790000000000000
        """
        cdef RealNumber x
        x = self._new()
        mpfr_frac(x.value, self.value, (<RealField>self._parent).rnd)
        return x

    ###########################################
    # Conversions
    ###########################################

    def __float__(self):
        return mpfr_get_d(self.value, (<RealField>self._parent).rnd)

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
        """
        Returns self as a Pari floating-point number.

        EXAMPLES:
            sage: RR(2.0)._pari_()
            2.000000000000000000

        The current Pari precision affects the printing of this number, but
        Pari does maintain the same 250-bit number on both 32-bit and 64-bit
        platforms.

            sage: RealField(250).pi()._pari_()
            3.141592653589793238462643383               # 32-bit
            3.1415926535897932384626433832795028842     # 64-bit
            sage: RR(0.0)._pari_()
            0.E-19
            sage: RR(-1.234567)._pari_()
            -1.2345669999999999700
            sage: RR(2.0).sqrt()._pari_()
            1.4142135623730951455
            sage: RR(2.0).sqrt()._pari_().python()
            1.41421356237309515
            sage: RR(2.0).sqrt()._pari_().python().prec()
            64
            sage: RealField(70)(pi)._pari_().python().prec()
            96                                         # 32-bit
            128                                        # 64-bit
            sage: for i in xrange(1, 1000):
            ...       assert(RR(i).sqrt() == RR(i).sqrt()._pari_().python())
        """
        # return sage.libs.pari.all.pari.new_with_bits_prec(str(self), (<RealField>self._parent).__prec)

        # This uses interfaces of MPFR and Pari which are documented
        # (and not marked subject-to-change).  It could be faster
        # by using internal interfaces of MPFR, which are documented
        # as subject-to-change.

        if mpfr_nan_p(self.value) or mpfr_inf_p(self.value):
            raise ValueError, 'Cannot convert NaN or infinity to Pari float'

        cdef int wordsize

        if sage.misc.misc.is_64_bit:
            wordsize = 64
        else:
            wordsize = 32

        cdef int prec
        prec = (<RealField>self._parent).__prec

        # We round up the precision to the nearest multiple of wordsize.
        cdef int rounded_prec
        rounded_prec = (self.prec() + wordsize - 1) & ~(wordsize - 1)

        # Yes, assigning to self works fine, even in Pyrex.
        if rounded_prec > prec:
            self = RealField(rounded_prec)(self)

        # Now we can extract the mantissa, and it will be normalized
        # (the most significant bit of the most significant word will be 1).
        cdef mpz_t mantissa
        cdef mp_exp_t exponent
        mpz_init(mantissa)

        exponent = mpfr_get_z_exp(mantissa, self.value)

        cdef GEN pari_float
        pari_float = cgetr(2 + rounded_prec / wordsize)

        mpz_export(&pari_float[2], NULL, 1, wordsize/8, 0, 0, mantissa)

        if mpfr_zero_p(self.value):
            setexpo(pari_float, -rounded_prec)
        else:
            setexpo(pari_float, exponent + rounded_prec - 1)
        setsigne(pari_float, mpfr_sgn(self.value))

        cdef PariInstance P
        P = sage.libs.pari.all.pari

        gen = P.new_gen(pari_float)

        mpz_clear(mantissa)

        return gen

    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    def is_NaN(self):
        return bool(mpfr_nan_p(self.value))

    def is_positive_infinity(self):
        """
        EXAMPLES:
            sage: a = RR('1.494') / RR(0); a
            +infinity
            sage: a.is_positive_infinity()
            True
            sage: a = -RR('1.494') / RR(0); a
            -infinity
            sage: RR(1.5).is_positive_infinity()
            False
            sage: a.is_positive_infinity()
            False
        """
        return bool(mpfr_inf_p(self.value) and mpfr_sgn(self.value) > 0)

    def is_negative_infinity(self):
        """
        EXAMPLES:
            sage: a = RR('1.494') / RR(0); a
            +infinity
            sage: a.is_negative_infinity()
            False
            sage: a = -RR('1.494') / RR(0); a
            -infinity
            sage: RR(1.5).is_negative_infinity()
            False
            sage: a.is_negative_infinity()
            True
        """
        return bool(mpfr_inf_p(self.value) and mpfr_sgn(self.value) < 0)

    def is_infinity(self):
        """
        EXAMPLES:
            sage: a = RR('1.494') / RR(0); a
            +infinity
            sage: a.is_infinity()
            True
            sage: a = -RR('1.494') / RR(0); a
            -infinity
            sage: a.is_infinity()
            True
            sage: RR(1.5).is_infinity()
            False
        """
        return bool(mpfr_inf_p(self.value))

    def __richcmp__(left, right, int op):
        return (<RingElement>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        cdef RealNumber self, x
        self = left
        x = right

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
            2.00000000000000
            sage: r.sqrt()^2 == r
            True

            sage: r = 4344
            sage: r.sqrt()
            65.9090282131363
            sage: r.sqrt()^2 == r
            True
            sage: r.sqrt()^2 - r
             0.000000000000000

            sage: r = -2.0
            sage: r.sqrt()
            1.41421356237310*I
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
            1.41421356237310*I
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_sqrt(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def cube_root(self):
        """
        Return the cubic root (defined over the real numbers) of self.

        EXAMPLES:
            sage: r = 125.0; r.cube_root()
            5.00000000000000
            sage: r = -119.0
            sage: r.cube_root()^3 - r       # illustrates precision loss
            -0.0000000000000142108547152020
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_cbrt(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def __pow(self, RealNumber exponent):
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_pow(x.value, self.value, exponent.value, (<RealField>self._parent).rnd)
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
            67.646297
            sage: a^a
            1.2971115
            sage: b = R(-1)
            sage: b^(1/2)
            1.0000000*I                    # 32-bit
            -1.0842022e-19 + 1.0000000*I   # 64-bit
        """
        cdef RealNumber x
        if not PY_TYPE_CHECK(self, RealNumber):
            return self.__pow__(float(exponent))
        if not PY_TYPE_CHECK(exponent, RealNumber):
            x = self
            exponent = x._parent(exponent)
        return self.__pow(exponent)

    def log(self, base='e'):
        """
        EXAMPLES:
            sage: R = RealField()
            sage: R(2).log()
            0.693147180559945
        """
        cdef RealNumber x
        if base == 'e':
            x = self._new()
            _sig_on
            mpfr_log(x.value, self.value, (<RealField>self._parent).rnd)
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
            4.00000000000000

            sage: r = 31.9; r.log2()
            4.99548451887751

            sage: r = 0.0
            sage: r.log2()
            -infinity
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_log2(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def log10(self):
        """
        Returns log to the base 10 of self

        EXAMPLES:
            sage: r = 16.0; r.log10()
            1.20411998265592
            sage: r.log() / log(10.0)
            1.20411998265592

            sage: r = 39.9; r.log10()
            1.60097289568675

            sage: r = 0.0
            sage: r.log10()
            -infinity

            sage: r = -1.0
            sage: r.log10()
            NaN

        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_log10(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def exp(self):
        r"""
        Returns $e^\code{self}$

        EXAMPLES:
            sage: r = 0.0
            sage: r.exp()
            1.00000000000000

            sage: r = 32.3
            sage: a = r.exp(); a
            106588847274864
            sage: a.log()
            32.3000000000000

            sage: r = -32.3
            sage: r.exp()
            0.00000000000000938184458849869
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_exp(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def exp2(self):
        """
        Returns $2^\code{self}$

        EXAMPLES:
            sage: r = 0.0
            sage: r.exp2()
            1.00000000000000

            sage: r = 32.0
            sage: r.exp2()
            4294967296.00000

            sage: r = -32.3
            sage: r.exp2()
            0.000000000189117248253021

        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_exp2(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def exp10(self):
        r"""
        Returns $10^\code{self}$

        EXAMPLES:
            sage: r = 0.0
            sage: r.exp10()
            1.00000000000000

            sage: r = 32.0
            sage: r.exp10()
            1.00000000000000e32

            sage: r = -32.3
            sage: r.exp10()
            5.01187233627276e-33
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_exp10(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def cos(self):
        """
        Returns the cosine of this number

        EXAMPLES:
            sage: t=RR.pi()/2
            sage: t.cos()
            6.12323399573677e-17
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_cos(x.value, self.value, (<RealField>self._parent).rnd)
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
            0.90929742682568169539601986591
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_sin(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def tan(self):
        """
        Returns the tangent of this number

        EXAMPLES:
            sage: q = RR.pi()/3
            sage: q.tan()
            1.73205080756888
            sage: q = RR.pi()/6
            sage: q.tan()
            0.577350269189626
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_tan(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def sincos(self):
        """
        Returns a pair consisting of the sine and cosine.

        EXAMPLES:
            sage: R = RealField()
            sage: t = R.pi()/6
            sage: t.sincos()
            (0.500000000000000, 0.866025403784439)
        """
        cdef RealNumber x,y
        x = self._new()
        y = self._new()
        _sig_on
        mpfr_sin_cos(x.value, y.value, self.value, (<RealField>self._parent).rnd)
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
        x = self._new()
        _sig_on
        mpfr_acos(x.value, self.value, (<RealField>self._parent).rnd)
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
            sage: i.asin() - q
            0.000000000000000
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_asin(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def atan(self):
        """
        Returns the inverse tangent of this number

        EXAMPLES:
            sage: q = RR.pi()/5
            sage: i = q.tan()
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_atan(x.value, self.value, (<RealField>self._parent).rnd)
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
            1.03446564009551
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_cosh(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def sinh(self):
        """
        Returns the hyperbolic sine of this number

        EXAMPLES:
            sage: q = RR.pi()/12
            sage: q.sinh()
            0.264800227602271
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_sinh(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def tanh(self):
        """
        Returns the hyperbolic tangent of this number

        EXAMPLES:
            sage: q = RR.pi()/11
            sage: q.tanh()
            0.278079429295850
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_tanh(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def acosh(self):
        """
        Returns the hyperbolic inverse cosine of this number

        EXAMPLES:
            sage: q = RR.pi()/2
            sage: i = q.cosh() ; i
            2.50917847865806
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_acosh(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def asinh(self):
        """
        Returns the hyperbolic inverse sine of this number

        EXAMPLES:
            sage: q = RR.pi()/7
            sage: i = q.sinh() ; i
            0.464017630492991
            sage: i.asinh() - q
            0.000000000000000
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_asinh(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x

    def atanh(self):
        """
        Returns the hyperbolic inverse tangent of this number

        EXAMPLES:
            sage: q = RR.pi()/7
            sage: i = q.tanh() ; i
            0.420911241048535
            sage: i.atanh() - q
            0.000000000000000
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_atanh(x.value, self.value, (<RealField>self._parent).rnd)
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
        if not isinstance(other, RealNumber) or other.parent() != self._parent:
            _other = self._parent(other)
        else:
            _other = other
        x = self._new()
        _sig_on
        mpfr_agm(x.value, self.value, _other.value, (<RealField>self._parent).rnd)
        _sig_off
        return x


    def erf(self):
        """
        Returns the value of the error function on self.

        EXAMPLES:
           sage: R = RealField()
           sage: R(6).erf()
           1.00000000000000
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_erf(x.value, self.value, (<RealField>self._parent).rnd)
        _sig_off
        return x


    def gamma(self):
        """
        The Euler gamma function. Return gamma of self.

        EXAMPLES:
           sage: R = RealField()
           sage: R(6).gamma()
           120.000000000000
           sage: R(1.5).gamma()
           0.886226925452758
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_gamma(x.value, self.value, (<RealField>self._parent).rnd)
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
            1.64493406684823
            sage: R.pi()^2/6
            1.64493406684823
            sage: R(-2).zeta()
            0.000000000000000
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
             <type 'sage.libs.pari.gen.gen'>
             sage: R(z)
             1.64493406684823
        """
        cdef RealNumber x
        x = self._new()
        _sig_on
        mpfr_zeta(x.value, self.value, (<RealField>self._parent).rnd)
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
             1.41421356237310
             sage: r.algdep(5)
             x^2 - 2                          # 32-bit
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
              1.41421356237310
              sage: r.algdep(5)
              x^2 - 2                        # 32-bit
              x^5 - x^4 - 2*x^3 + x^2 + 2    # 64-bit
        """
        return sage.rings.arith.algdep(self,n)

    def nth_root(self, int n):
        r"""
        Returns an $n^{th}$ root of self.

        INPUT:
            n -- A positive number, rounded down to the nearest integer.
                 Note that $n$ should be less than $\code{sys.maxint}$.

        EXAMPLES:
            sage: R = RealField()
            sage: R(8).nth_root(3)
            2.00000000000000
            sage: R(8).nth_root(3.7)    # illustrate rounding down
            2.00000000000000
            sage: R(-8).nth_root(3)
            -2.00000000000000
            sage: R(0).nth_root(3)
            0.000000000000000
            sage: R(32).nth_root(-1)
            Traceback (most recent call last):
            ...
            ValueError: n must be nonnegative
            sage: R(32).nth_root(1.0)
            32.0000000000000

        Note that for negative numbers, any even root returns NaN
            sage: R(-2).nth_root(6)
            NaN

        The $n^{th}$ root of 0 is defined to be 0, for any $n$
            sage: R(0).nth_root(6)
            0.000000000000000

            sage: R(0).nth_root(7)
            0.000000000000000

        AUTHOR: Didier Deshommes (2007-02)
        REFEREE: David Harvey

        TODO: (trac \#294) the underlying mpfr_root function is unforgivably
        slow when n is large. e.g. RealNumber(8).nth_root(100000) is very slow.
        This should be investigated further and possibly discussed with the
        mpfr developers.
        """
        cdef RealNumber x

        if n < 0:
            raise ValueError, "n must be nonnegative"

        x = self._new()
        _sig_on
        mpfr_root(x.value, self.value, n, (<RealField>self._parent).rnd)
        _sig_off
        return x


RR = RealField()


def create_RealNumber(s, int base=10, int pad=0, rnd="RNDN", min_prec=53):
    r"""
    Return the real number defined by the string s as an element of
    \code{RealField(prec=n)}, where n potentially has slightly more
    (controlled by pad) bits than given by s.

    INPUT:
        s -- a string that defines a real number (or something whose
             string representation defines a number)
        base -- an integer between 2 and 36
        pad -- an integer >= 0.
        rnd -- rounding mode: RNDN, RNDZ, RNDU, RNDD
        min_prec -- number will have at least this many bits of precision, no matter what.

    EXAMPLES:
        sage: RealNumber('2.3')
        2.30000000000000
        sage: RealNumber(10)
        10.0000000000000
        sage: RealNumber('1.0000000000000000000000000000000000')
        1.0000000000000000000000000000000000
    """
    if not isinstance(s, str):
        s = str(s)
    if base == 10:
        bits = int(3.32192*len(s))
    else:
        bits = int(math.log(base,2)*len(s))
    R = RealField(prec=max(bits+pad, min_prec), rnd=rnd)
    return RealNumber(R, s, base)



def is_RealField(x):
    return PY_TYPE_CHECK(x, RealField)

def is_RealNumber(x):
    return PY_TYPE_CHECK(x, RealNumber)

def __create__RealField_version0(prec, sci_not, rnd):
    return RealField(prec, sci_not, rnd)

def __create__RealNumber_version0(parent, x, base=10):
    return RealNumber(parent, x, base=base)
