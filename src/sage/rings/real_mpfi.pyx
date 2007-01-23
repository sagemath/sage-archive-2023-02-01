"""
Aribtrary Precision Real Numbers -- Interval Arithmetic

AUTHORS: Kyle Schalm <kschalm@math.utexas.edu> (2005-09)
         William Stein <wstein@gmail.com>: bug fixes, examples, maintenance
         Didier Deshommes <dfdeshom@gmail.com> (2006-03-19): examples
         David Harvey (2006-09-20): compatibility with Element._parent
         William Stein (2006-10): default printing truncates to avoid base-2
              rounding confusing (fix suggested by Bill Hart)
         Carl Witty (2007-01-21): COPIED THIS FILE FROM real_mpfr.pyx;
               changing it to use mpfi rather than mpfr
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
# TypeError: Unable to convert x (='1/2') to mpfi.
#
# sage: 1 + R(42)
# _49 = 1
#*****************************************************************************

import math # for log
import sys

include '../ext/interrupt.pxi'
include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"

cimport sage.rings.ring
import  sage.rings.ring

cimport sage.structure.element
from sage.structure.element cimport RingElement, Element, ModuleElement
import  sage.structure.element

cimport real_mpfr
from real_mpfr cimport RealField, RealNumber
import real_mpfr

import sage.structure.coerce
import operator

from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational

import sage.rings.complex_field

import sage.rings.infinity

from sage.structure.parent_gens cimport ParentWithGens

cdef class RealIntervalFieldElement(sage.structure.element.RingElement)

#*****************************************************************************
#
#       Implementation
#
#*****************************************************************************

#*****************************************************************************
#
#       Real Field
#
#*****************************************************************************
# The real field is in Pyrex, so mpfi elements will have access to
# their parent via direct C calls, which will be faster.

cdef class RealIntervalField(sage.rings.ring.Field):
    """
    RealIntervalField(prec, sci_not, rnd):

    INPUT:
        prec -- (integer) precision; default = 53
                prec is the number of bits used to represent the
                mantissa of a floating-point number.  The
                precision can be any integer between mpfr_prec_min()
                and mpfr_prec_max(). In the current implementation,
                mpfr_prec_min() is equal to 2.

        sci_not -- (default: False) whether or not to display
                using scientific notation

    EXAMPLES:
        sage: RealIntervalField(10)
        Real Interval Field with 10 bits of precision
        sage: RealIntervalField()
        Real Interval Field with 53 bits of precision
        sage: RealIntervalField(100000)
        Real Interval Field with 100000 bits of precision

    NOTE: The default precision is 53, since according to the GMP
       manual: 'mpfr should be able to exactly reproduce all
       computations with double-precision machine floating-point
       numbers (double type in C), except the default exponent
       range is much wider and subnormal numbers are not
       implemented.'
    """

    def __init__(self, int prec=53, int sci_not=0):
        if prec < MPFR_PREC_MIN or prec > MPFR_PREC_MAX:
            raise ValueError, "prec (=%s) must be >= %s and <= %s."%(
                prec, MPFR_PREC_MIN, MPFR_PREC_MAX)
        self.__prec = prec
        self.sci_not = sci_not
        self.__lower_field = RealField(prec, sci_not, "RNDD")
        self.__upper_field = RealField(prec, sci_not, "RNDU")
        ParentWithGens.__init__(self, self, tuple([]), False)

    cdef RealIntervalFieldElement _new(self):
        """
        Return a new real number with parent self.
        """
        cdef RealIntervalFieldElement x
        x = PY_NEW(RealIntervalFieldElement)
        x._parent = self
        mpfi_init2(x.value, self.__prec)
        x.init = 1
        return x

    def _repr_(self):
        s = "Real Interval Field with %s bits of precision"%self.__prec
        return s

    def _latex_(self):
        # Is there a standard notation for this?
        return "\\R_I"

    def is_exact(self):
        return False

    def __call__(self, x, base=10):
        """
        Coerce x into this real field.

        EXAMPLES:
            sage: R = RealIntervalField(20)
            sage: R('1.234')
            [1.2339992 ... 1.2340012]
            sage: R('2', base=2)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert x (='2') to real number.
            sage: a = R('1.1001', base=2); a
            [1.5625000 ... 1.5625000]
            sage: a.str(2)
            '[1.1001000000000000000 ... 1.1001000000000000000]'
        """
        if hasattr(x, '_mpfi_'):
            return x._mpfi_(self)
        return RealIntervalFieldElement(self, x, base)

    cdef _coerce_c_impl(self, x):
        """
        Canonical coercion of x to this mpfi real field.

        The rings that canonically coerce to this mpfi real field are:
             * this real field itself
             * any mpfr real field with precision that is as large as this one
             * any other mpfi real field with precision that is as large as this one
             * int, long, integer, and rational rings.
             * real mathematical constants
        """
        if isinstance(x, real_mpfr.RealNumber):
            P = x.parent()
            if (<RealField> P).__prec >= self.__prec:
                return self(x)
            else:
                raise TypeError, "Canonical coercion from lower to higher precision not defined"
        if isinstance(x, RealIntervalFieldElement):
            P = x.parent()
            if (<RealIntervalField> P).__prec >= self.__prec:
                return self(x)
            else:
                raise TypeError, "Canonical coercion from lower to higher precision not defined"
        if isinstance(x, (Integer, Rational)):
            return self(x)
        import sage.functions.constants
        return self._coerce_try(x, [sage.functions.constants.ConstantRing])

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: RealIntervalField(10) == RealIntervalField(11)
            False
            sage: RealIntervalField(10) == RealIntervalField(10)
            True
            sage: RealIntervalField(10,sci_not=True) == RealIntervalField(10,sci_not=False)
            False
            sage: RealIntervalField(10) == IntegerRing()
            False
        """
        if not isinstance(other, RealIntervalField):
            return -1
        cdef RealIntervalField _other
        _other = other  # to access C structure
        if self.__prec == _other.__prec and self.sci_not == _other.sci_not:
            return 0
        return 1

    def __reduce__(self):
        """
        EXAMPLES:
            sage: R = RealIntervalField(sci_not=1, prec=200)
            sage: loads(dumps(R)) == R
            True
        """
        return sage.rings.real_interval_field.__reduce__RealIntervalField, \
                (self.__prec, self.sci_not)

    def gen(self, i=0):
        if i == 0:
            return self(1)
        else:
            raise IndexError

#    def complex_field(self):
#        """
#        Return complex field of the same precision.
#        """
#        return sage.rings.complex_field.ComplexField(self.prec())

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
            sage: RealIntervalField(10).is_atomic_repr()
            True
        """
        return True

    def is_finite(self):
        """
        Returns False, since the field of real numbers is not finite.

        EXAMPLES:
            sage: RealIntervalField(10).is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Returns 0, since the field of real numbers has characteristic 0.

        EXAMPLES:
            sage: RealIntervalField(10).characteristic()
            0
        """
        return 0

    def name(self):
        return "IntervalRealIntervalField%s"%(self.__prec)

    def __hash__(self):
        return hash(self.name())

    def precision(self):
        return self.__prec

    def prec(self):
        return self.__prec

    # int mpfi_const_pi (mpfi_ptr)
    def pi(self):
        """
        Returns pi to the precision of this field.

        EXAMPLES:
            sage: R = RealIntervalField(100)
            sage: R.pi()
            [3.1415926535897932384626433832793 ... 3.1415926535897932384626433832825]
            sage: R.pi().sqrt()/2
            [0.88622692545275801364908374166983 ... 0.88622692545275801364908374167142]
            sage: R = RealIntervalField(150)
            sage: R.pi().sqrt()/2
            [0.88622692545275801364908374167057259139877472759 ... 0.88622692545275801364908374167057259139877472830]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_const_pi(x.value)
        return x


    # int mpfi_const_euler (mpfi_ptr)
    def euler_constant(self):
        """
        Returns Euler's gamma constant to the precision of this field.

        EXAMPLES:
            sage: RealIntervalField(100).euler_constant()
            [0.57721566490153286060651209008233 ... 0.57721566490153286060651209008313]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_const_euler(x.value)
        return x

#     def catalan_constant(self):
#         """
#         Returns Catalan's constant to the precision of this field.

#         EXAMPLES:
#             sage: RealIntervalField(100).catalan_constant()
#             0.91596559417721901505460351493
#         """
#         cdef RealIntervalFieldElement x
#         x = self._new()
#         mpfr_const_catalan(x.value, self.rnd)
#         return x

    def log2(self):
        """
        Returns log(2) to the precision of this field.

        EXAMPLES:
            sage: R=RealIntervalField(100)
            sage: R.log2()
            [0.69314718055994530941723212145798 ... 0.69314718055994530941723212145878]
            sage: R(2).log()
            [0.69314718055994530941723212145798 ... 0.69314718055994530941723212145878]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_const_log2(x.value)
        return x

# MPFI does not have factorial
#     def factorial(self, int n):
#         """
#         Return the factorial of the integer n as a real number.
#         """
#         cdef RealIntervalFieldElement x
#         if n < 0:
#             raise ArithmeticError, "n must be nonnegative"
#         x = self._new()
#         mpfr_fac_ui(x.value, n, self.rnd)
#         return x

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
            sage: R = RealIntervalField()
            sage: R.zeta()
            [-1.0000000000000000 ... -1.0000000000000000]
            sage: R.zeta(1)
            [1.0000000000000000 ... 1.0000000000000000]
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

R = RealIntervalField()

#*****************************************************************************
#
#     RealIntervalFieldElement -- element of Real Field
#
#
#
#*****************************************************************************
cdef class RealIntervalFieldElement(sage.structure.element.RingElement):
    """
    A real number interval.
    """
    cdef RealIntervalFieldElement _new(self):
        """
        Return a new real interval with same parent as self.
        """
        cdef RealIntervalFieldElement x
        x = PY_NEW(RealIntervalFieldElement)
        x._parent = self._parent
        mpfi_init2(x.value, (<RealIntervalField>self._parent).__prec)
        x.init = 1
        return x

    def __init__(self, RealIntervalField parent, x=0, int base=10, special=None):
        """
        Create a real number.  Should be called by first creating
        a RealIntervalField, as illustrated in the examples.

        EXAMPLES:
            sage: R = RealIntervalField()
            sage: R('1.2456')
            [1.2455999999999998 ... 1.2456000000000001]
            sage: R = RealIntervalField(3)
            sage: R('1.2456')
            [1.0 ... 1.3]

        EXAMPLE: Rounding
            sage: w = RealIntervalField(3)(5/2)
            sage: RealIntervalField(2)(w).str(2)
            '[10 ... 11]'

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
        self.init = 0
        if parent is None:
            raise TypeError
        self._parent = parent
        mpfi_init2(self.value, parent.__prec)
        self.init = 1
        if x is None: return
        cdef RealIntervalFieldElement _x, n, d
        cdef RealNumber rn
        cdef Rational rat
        cdef int ix
        if PY_TYPE_CHECK(x, RealIntervalFieldElement):
            _x = x  # so we can get at x.value
            mpfi_set(self.value, _x.value)
        elif PY_TYPE_CHECK(x, RealNumber):
            rn = x
            mpfi_set_fr(self.value, <mpfr_t> rn.value)
        elif PY_TYPE_CHECK(x, sage.rings.rational.Rational):
            rat = x
            mpfi_set_q(self.value, <mpq_t> rat.value)
        elif PY_TYPE_CHECK(x, int):
            ix = x
            mpfi_set_si(self.value, ix)
        else:
            s = str(x).replace('...', ',').replace(' ','').replace('+infinity', '@inf@').replace('-infinity','-@inf@')
            if mpfi_set_str(self.value, s, base):
#                 if s == 'NaN' or s == '@NaN@' or s == '[... NaN ...]' or s == '[... @NaN@ ...]':
#                     mpfi_set_nan(self.value)
#                 elif s == '+infinity':
#                     mpfr_set_inf(self.value, 1)
#                 elif s == '-infinity':
#                     mpfr_set_inf(self.value, -1)
#                 else:
                    raise TypeError, "Unable to convert x (='%s') to real number."%s


    def __reduce__(self):
        """
        EXAMPLES:
            sage: R = RealIntervalField(sci_not=1, prec=200)
            sage: b = R('393.39203845902384098234098230948209384028340')
            sage: loads(dumps(b)) == b
            True
            sage: b = R(1)/R(0); b
            [+infinity ... +infinity]
            sage: loads(dumps(b)) == b
            True
            sage: b = R(-1)/R(0); b
            [-infinity ... -infinity]
            sage: loads(dumps(b)) == b
            True
            sage: b = R('[2 ... 3]'); b
            [2.0000000000000000000000000000000000000000000000000000000000000e0 ... 3.0000000000000000000000000000000000000000000000000000000000000e0]
            sage: loads(dumps(b)) == b
            True
        """
        s = self.str(32, no_sci=False, e='@')
        return (sage.rings.real_interval_field.__reduce__RealIntervalFieldElement, (self._parent, s, 32))

    def  __dealloc__(self):
        if self.init:
            mpfi_clear(self.value)

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
            '1.39394945939999'
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
            sage: R = RealIntervalField()
            sage: a = R('1.2456')
            sage: a.parent()
            Real Interval Field with 53 bits of precision
        """
        return self._parent

    # MPFR had an elaborate "truncation" scheme to avoid printing
    # inaccurate-looking results; this has been removed for MPFI,
    # because I think it's less confusing to get such results than to
    # see [0.333 ... 0.333] for an interval with unequal left and right
    # sides.
    def str(self, int base=10, no_sci=None, e='e'):
        """
        INPUT:
             base -- base for output
             no_sci -- if True do not print using scientific notation; if False
                       print with scientific notation; if None (the default), print how the parent prints.
             e - symbol used in scientific notation

        EXAMPLES:
            sage: a = 61/3.0; a
            20.3333333333333
            sage: a.str(truncate=False)
            '20.333333333333332'
            sage: a.str(2)
            '10100.010101010101010101010101010101010101010101010101'
            sage: a.str(no_sci=False)
            '2.03333333333333e1'
        """
        if base < 2 or base > 36:
            raise ValueError, "the base (=%s) must be between 2 and 36"%base
        if mpfi_nan_p(self.value):
            if base >= 24:
                return "[... @NaN@ ...]"
            else:
                return "[... NaN ...]"

        t1 = self.lower().str(base=base, no_sci=no_sci, e=e, truncate=False)
        t2 = self.upper().str(base=base, no_sci=no_sci, e=e, truncate=False)

        return "[%s ... %s]"%(t1, t2)

    def __copy__(self):
        """
        Return copy of self -- since self is immutable, we just return self again.

        EXAMPLES:
            sage: a = 3.5
            sage: copy(a) is  a
            True
        """
        return self

    # Interval-specific functions
    def lower(self):
        """
        Returns the lower bound of this interval

        EXAMPLES:
            sage: R = RealIntervalField(13)
            sage: R.pi().lower().str(truncate=False)
            '3.1411'
        """
        cdef RealNumber x
        x = (<RealIntervalField>self._parent).__lower_field._new()
        mpfi_get_left(<mpfr_t> x.value, self.value)
        return x

    def upper(self):
        """
        Returns the upper bound of this interval

        EXAMPLES:
            sage: R = RealIntervalField(13)
            sage: R.pi().upper().str(truncate=False)
            '3.1417'
        """
        cdef RealNumber x
        x = (<RealIntervalField>self._parent).__upper_field._new()
        mpfi_get_right(<mpfr_t> x.value, self.value)
        return x

#     def integer_part(self):
#         """
#         If in decimal this number is written n.defg, returns n.

#         OUTPUT:
#             -- a SAGE Integer

#         EXAMPLE:
#             sage: a = 119.41212
#             sage: a.integer_part()
#             119
#         """
#         s = self.str(base=32, no_sci=True)
#         i = s.find(".")
#         return Integer(s[:i], base=32)

    ########################
    #   Basic Arithmetic
    ########################

    cdef ModuleElement _add_c_impl(self, ModuleElement other):
        """
        Add two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealIntervalField()
            sage: R(-1.5) + R(2.5)
            [1.0000000000000000 ... 1.0000000000000000]
            sage: R('-1.3') + R('2.3')
            [0.99999999999999977 ... 1.0000000000000005]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_add(x.value, self.value, (<RealIntervalFieldElement>other).value)
        return x

    def __invert__(self):
        return self._parent(1) / self

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        Subtract two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealIntervalField()
            sage: R(-1.5) - R(2.5)
            [-4.0000000000000000 ... -4.0000000000000000]
            sage: R('-1.3') - R('2.7')
            [-4.0000000000000009 ... -3.9999999999999995]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_sub(x.value, self.value, (<RealIntervalFieldElement>right).value)
        return x

    cdef RingElement _mul_c_impl(self, RingElement right):
        """
        Multiply two real numbers with the same parent.

        EXAMPLES:
            sage: R = RealIntervalField()
            sage: R(-1.5) * R(2.5)
            [-3.7500000000000000 ... -3.7500000000000000]
            sage: R('-1.3') * R('2.3')
            [-2.9900000000000007 ... -2.9899999999999993]

        If two elements have different precision, arithmetic
        operations are performed after coercing to the lower
        precision.

            sage: R20 = RealIntervalField(20)
            sage: R100 = RealIntervalField(100)
            sage: a = R20('393.3902834028345')
            sage: b = R100('393.3902834028345')
            sage: a
            [393.39013 ... 393.39063]
            sage: b
            [393.39028340283449999999999999970 ... 393.39028340283450000000000000011]
            sage: a*b
            [154755.75 ... 154756.25]
            sage: b*a
            [154755.75 ... 154756.25]
            sage: parent(b*a)
            Real Interval Field with 20 bits of precision
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_mul(x.value, self.value, (<RealIntervalFieldElement>right).value)
        return x


    cdef RingElement _div_c_impl(self, RingElement right):
        """
        Divide self by other, where both are real numbers with the same parent.

        EXAMPLES:
            sage: R = RealIntervalField()
            sage: R(1)/R(3)
            [0.33333333333333331 ... 0.33333333333333338]
            sage: R(1)/R(0)
            [+infinity ... +infinity]

            sage: R(-1.5) / R(2.5)
            [-0.60000000000000009 ... -0.59999999999999997]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_div((<RealIntervalFieldElement>x).value, self.value,
                 (<RealIntervalFieldElement>right).value)
        return x

    cdef ModuleElement _neg_c_impl(self):
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_neg(x.value, self.value)
        return x

    def __abs__(self):
        return self.abs()

    cdef RealIntervalFieldElement abs(RealIntervalFieldElement self):
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_abs(x.value, self.value)
        return x

    # Bit shifting
    def _lshift_(RealIntervalFieldElement self, n):
        cdef RealIntervalFieldElement x
        if n > sys.maxint:
            raise OverflowError, "n (=%s) must be <= %s"%(n, sys.maxint)
        x = self._new()
        mpfi_mul_2exp(x.value, self.value, n)
        return x

    def __lshift__(x, y):
        """
        EXAMPLES:
            sage: 1.0 << 32
            4294967296.00000
        """
        if isinstance(x, RealIntervalFieldElement) and isinstance(y, (int,long, Integer)):
            return x._lshift_(y)
        return sage.structure.coerce.bin_op(x, y, operator.lshift)

    def _rshift_(RealIntervalFieldElement self, n):
        if n > sys.maxint:
            raise OverflowError, "n (=%s) must be <= %s"%(n, sys.maxint)
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_div_2exp(x.value, self.value, n)
        return x

    def __rshift__(x, y):
        """
        EXAMPLES:
            sage: 1024.0 >> 7
            8.00000000000000
        """
        if isinstance(x, RealIntervalFieldElement) and isinstance(y, (int,long,Integer)):
            return x._rshift_(y)
        return sage.structure.coerce.bin_op(x, y, operator.rshift)

    def multiplicative_order(self):
        if self == 1:
            return 1
        elif self == -1:
            return -1
        return sage.rings.infinity.infinity

#     def sign(self):
#         return mpfr_sgn(self.value)

    def prec(self):
        return (<RealIntervalField>self._parent).__prec

    ###################
    # Rounding etc
    ###################

    # Not implemented on intervals (for good reason!)
#     def round(self):
#         """
#         Rounds self to the nearest real number. There are 4
#         rounding modes. They are

#         EXAMPLES:
#             RNDN -- round to nearest:

#             sage: R = RealIntervalField(20,False,'RNDN')
#             sage: R(22.454)
#             22.454
#             sage: R(22.491)
#             22.490

#             RNDZ -- round towards zero:
#             sage: R = RealIntervalField(20,False,'RNDZ')
#             sage: R(22.454)
#             22.453
#             sage: R(22.491)
#             22.490

#             RNDU -- round towards plus infinity:
#             sage: R = RealIntervalField(20,False,'RNDU')
#             sage: R(22.454)
#             22.454
#             sage: R(22.491)
#             22.491

#             RNDU -- round towards minus infinity:
#             sage: R = RealIntervalField(20,False,'RNDD')
#             sage: R(22.454)
#             22.453
#             sage: R(22.491)
#             22.490
#         """
#         cdef RealIntervalFieldElement x
#         x = self._new()
#         mpfr_round(x.value, self.value)
#         return x

#     def floor(self):
#         """
#         Returns the floor of this number

#         EXAMPLES:
#             sage: R = RealIntervalField()
#             sage: (2.99).floor()
#             2
#             sage: (2.00).floor()
#             2
#             sage: floor(RR(-5/2))
#             -3
#         """
#         cdef RealIntervalFieldElement x
#         x = self._new()
#         mpfr_floor(x.value, self.value)
#         return x.integer_part()

#     def ceil(self):
#         """
#         Returns the ceiling of this number

#         OUTPUT:
#             integer

#         EXAMPLES:
#             sage: (2.99).ceil()
#             3
#             sage: (2.00).ceil()
#             2
#             sage: (2.01).ceil()
#             3
#         """
#         cdef RealIntervalFieldElement x
#         x = self._new()
#         mpfr_ceil(x.value, self.value)
#         return x.integer_part()

#     def ceiling(self):
#         return self.ceil()

#     def trunc(self):
#         """
#         Truncates this number

#         EXAMPLES:
#             sage: (2.99).trunc()
#             2.00000000000000
#             sage: (-0.00).trunc()
#             -0.000000000000000
#             sage: (0.00).trunc()
#             0.000000000000000
#         """
#         cdef RealIntervalFieldElement x
#         x = self._new()
#         mpfr_trunc(x.value, self.value)
#         return x

#     def frac(self):
#         """
#         frac returns a real number > -1 and < 1. it satisfies the
#         relation:
#             x = x.trunc() + x.frac()

#         EXAMPLES:
#             sage: (2.99).frac()
#             0.990000000000000
#             sage: (2.50).frac()
#             0.500000000000000
#             sage: (-2.79).frac()
#             -0.790000000000000
#         """
#         cdef RealIntervalFieldElement x
#         x = self._new()
#         mpfr_frac(x.value, self.value, (<RealIntervalField>self._parent).rnd)
#         return x

    ###########################################
    # Conversions
    ###########################################

#     def __float__(self):
#         return mpfr_get_d(self.value, (<RealIntervalField>self._parent).rnd)

#     def __int__(self):
#         """
#         Returns integer truncation of this real number.
#         """
#         s = self.str(32)
#         i = s.find('.')
#         return int(s[:i], 32)

#     def __long__(self):
#         """
#         Returns long integer truncation of this real number.
#         """
#         s = self.str(32)
#         i = s.find('.')
#         return long(s[:i], 32)

#     def __complex__(self):
#         return complex(float(self))

#     def _complex_number_(self):
#         return sage.rings.complex_field.ComplexField(self.prec())(self)

#     def _pari_(self):
#         return sage.libs.pari.all.pari.new_with_bits_prec(str(self), (<RealIntervalField>self._parent).__prec)


    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    def is_NaN(self):
        return bool(mpfi_nan_p(self.value))

    def __richcmp__(left, right, int op):
        return (<RingElement>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        cdef RealIntervalFieldElement self, x
        self = left
        x = right

        cdef int i
        i = mpfi_cmp(self.value, x.value)
        if i < 0:
            return -1
        elif i == 0:
            return 0
        else:
            return 1


    ############################
    # Special Functions
    ############################

    # XXX This is wrong...
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
            False
            sage: r.sqrt()^2 - r
             -0.000000000000909494701772928

            sage: r = -2.0
            sage: r.sqrt()
            1.41421356237309*I
            """
        # XXX Probably does the wrong thing with intervals spanning 0...
        if self >= 0:
            return self.square_root()
        raise ValueError, "self (=%s) is not >= 0"%self


    def square_root(self):
        """
        Return a square root of self.  An interval will always be
        returned (though it will be NaN if self is negative).

        EXAMPLES:
            sage: r = -2.0
            sage: r.square_root()
            NaN
            sage: r.sqrt()
            1.41421356237309*I
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_sqrt(x.value, self.value)
        _sig_off
        return x

# MPFI does not have cbrt.
#     def cube_root(self):
#         """
#         Return the cubic root (defined over the real numbers) of self.

#         EXAMPLES:
#             sage: r = 125.0; r.cube_root()
#             5.00000000000000
#             sage: r = -119.0
#             sage: r.cube_root()^3 - r       # illustrates precision loss
#             -0.0000000000000142108547152020
#         """
#         cdef RealIntervalFieldElement x
#         x = self._new()
#         _sig_on
#         mpfr_cbrt(x.value, self.value, (<RealIntervalField>self._parent).rnd)
#         _sig_off
#         return x

# MPFI does not have pow.
#     def __pow(self, RealIntervalFieldElement exponent):
#         cdef RealIntervalFieldElement x
#         x = self._new()
#         _sig_on
#         mpfr_pow(x.value, self.value, exponent.value, (<RealIntervalField>self._parent).rnd)
#         _sig_off
#         if mpfr_nan_p(x.value):
#             return self._complex_number_()**exponent._complex_number_()
#         return x

#     def __pow__(self, exponent, modulus):
#         """
#         Compute self raised to the power of exponent, rounded in
#         the direction specified by the parent of self.

#         If the result is not a real number, self and the exponent are
#         both coerced to complex numbers (with sufficient precision),
#         then the exponentiation is computed in the complex numbers.
#         Thus this function can return either a real or complex number.

#         EXAMPLES:
#             sage: R = RealIntervalField(30)
#             sage: a = R('1.23456')
#             sage: a^20
#             67.646297
#             sage: a^a
#             1.2971114
#             sage: b = R(-1)
#             sage: b^(1/2)
#             1.0000000*I                   # 32-bit
#             -0.00000000000000000010842021 + 0.99999999*I   # 64-bit
#         """
#         cdef RealIntervalFieldElement x
#         if not PY_TYPE_CHECK(self, RealIntervalFieldElement):
#             return self.__pow__(float(exponent))
#         if not PY_TYPE_CHECK(exponent, RealIntervalFieldElement):
#             x = self
#             exponent = x._parent(exponent)
#         return self.__pow(exponent)

    def log(self, base='e'):
        """
        EXAMPLES:
            sage: R = RealIntervalField()
            sage: R(2).log()
            [0.69314718055994528 ... 0.69314718055994540]
        """
        cdef RealIntervalFieldElement x
        if base == 'e':
            x = self._new()
            _sig_on
            mpfi_log(x.value, self.value)
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
            4.99548451887750

            sage: r = 0.0
            sage: r.log2()
            -infinity
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_log2(x.value, self.value)
        _sig_off
        return x

    def log10(self):
        """
        Returns log to the base 10 of self

        EXAMPLES:
            sage: r = 16.0; r.log10()
            1.20411998265592
            sage: r.log() / log(10)
            1.20411998265592

            sage: r = 39.9; r.log10()
            1.60097289568674

            sage: r = 0.0
            sage: r.log10()
            -infinity

            sage: r = -1.0
            sage: r.log10()
            NaN

        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_log10(x.value, self.value)
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
            32.2999999999999

            sage: r = -32.3
            sage: r.exp()
            0.00000000000000938184458849868
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_exp(x.value, self.value)
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
            0.000000000189117248253020

        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_exp2(x.value, self.value)
        _sig_off
        return x

# MPFI does not have exp10.  (Could easily be synthesized if anybody cares.)
#     def exp10(self):
#         r"""
#         Returns $10^\code{self}$

#         EXAMPLES:
#             sage: r = 0.0
#             sage: r.exp10()
#             1.00000000000000

#             sage: r = 32.0
#             sage: r.exp10()
#             100000000000000000000000000000000

#             sage: r = -32.3
#             sage: r.exp10()
#             0.00000000000000000000000000000000501187233627275
#         """
#         cdef RealIntervalFieldElement x
#         x = self._new()
#         _sig_on
#         mpfr_exp10(x.value, self.value, (<RealIntervalField>self._parent).rnd)
#         _sig_off
#         return x

    def cos(self):
        """
        Returns the cosine of this number

        EXAMPLES:
            sage: t=RR.pi()/2
            sage: t.cos()
            0.0000000000000000612323399573676
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_cos(x.value, self.value)
        _sig_off
        return x

    def sin(self):
        """
        Returns the sine of this number

        EXAMPLES:
            sage: R = RealIntervalField(100)
            sage: R(2).sin()
            [0.90929742682568169539601986591150 ... 0.90929742682568169539601986591230]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_sin(x.value, self.value)
        _sig_off
        return x

    def tan(self):
        """
        Returns the tangent of this number

        EXAMPLES:
            sage: q = RR.pi()/3
            sage: q.tan()
            1.73205080756887
            sage: q = RR.pi()/6
            sage: q.tan()
            0.577350269189625
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_tan(x.value, self.value)
        _sig_off
        return x

# MPFI does not have sincos
#     def sincos(self):
#         """
#         Returns a pair consisting of the sine and cosine.

#         EXAMPLES:
#             sage: R = RealIntervalField()
#             sage: t = R.pi()/6
#             sage: t.sincos()
#             (0.499999999999999, 0.866025403784438)
#         """
#         cdef RealIntervalFieldElement x,y
#         x = self._new()
#         y = self._new()
#         _sig_on
#         mpfi_sin_cos(x.value, y.value, self.value)
#         _sig_off
#         return x,y

    def acos(self):
        """
        Returns the inverse cosine of this number

        EXAMPLES:
            sage: q = RR.pi()/3
            sage: i = q.cos()
            sage: i.acos() == q
            True
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_acos(x.value, self.value)
        _sig_off
        return x

    def asin(self):
        """
        Returns the inverse sine of this number

        EXAMPLES:
            sage: q = RR.pi()/5
            sage: i = q.sin()
            sage: i.asin() == q
            False
            sage: i.asin() - q
            -0.000000000000000111022302462515
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_asin(x.value, self.value)
        _sig_off
        return x

    def atan(self):
        """
        Returns the inverse tangent of this number

        EXAMPLES:
            sage: q = RR.pi()/5
            sage: i = q.tan()
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_atan(x.value, self.value)
        _sig_off
        return x

    def cosh(self):
        """
        Returns the hyperbolic cosine of this number

        EXAMPLES:
            sage: q = RR.pi()/12
            sage: q.cosh()
            1.03446564009551
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_cosh(x.value, self.value)
        _sig_off
        return x

    def sinh(self):
        """
        Returns the hyperbolic sine of this number

        EXAMPLES:
            sage: q = RR.pi()/12
            sage: q.sinh()
            0.264800227602270
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_sinh(x.value, self.value)
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
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_tanh(x.value, self.value)
        _sig_off
        return x

    def acosh(self):
        """
        Returns the hyperbolic inverse cosine of this number

        EXAMPLES:
            sage: q = RR.pi()/2
            sage: i = q.cosh() ; i
            2.50917847865805
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_acosh(x.value, self.value)
        _sig_off
        return x

    def asinh(self):
        """
        Returns the hyperbolic inverse sine of this number

        EXAMPLES:
            sage: q = RR.pi()/7
            sage: i = q.sinh() ; i
            0.464017630492990
            sage: i.asinh() - q
            -0.0000000000000000555111512312578
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_asinh(x.value, self.value)
        _sig_off
        return x

    def atanh(self):
        """
        Returns the hyperbolic inverse tangent of this number

        EXAMPLES:
            sage: q = RR.pi()/7
            sage: i = q.tanh() ; i
            0.420911241048534
            sage: i.atanh() - q
            -0.0000000000000000555111512312578
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_atanh(x.value, self.value)
        _sig_off
        return x

# MPFI does not have: agm, erf, gamma, zeta
#     def agm(self, other):
#         """
#         Return the arithmetic-geometric mean of self and other. The
#         arithmetic-geometric mean is the common limit of the sequences
#         $u_n$ and $v_n$, where $u_0$ is self, $v_0$ is other,
#         $u_{n+1}$ is the arithmetic mean of $u_n$ and $v_n$, and
#         $v_{n+1}$ is the geometric mean of u_n and v_n. If any operand
#         is negative, the return value is \code{NaN}.
#         """
#         cdef RealIntervalFieldElement x, _other
#         if not isinstance(other, RealIntervalFieldElement) or other.parent() != self._parent:
#             _other = self._parent(other)
#         else:
#             _other = other
#         x = self._new()
#         _sig_on
#         mpfi_agm(x.value, self.value, _other.value)
#         _sig_off
#         return x


#     def erf(self):
#         """
#         Returns the value of the error function on self.

#         EXAMPLES:
#            sage: R = RealIntervalField()
#            sage: R(6).erf()
#            0.999999999999999
#         """
#         cdef RealIntervalFieldElement x
#         x = self._new()
#         _sig_on
#         mpfi_erf(x.value, self.value)
#         _sig_off
#         return x


#     def gamma(self):
#         """
#         The Euler gamma function. Return gamma of self.

#         EXAMPLES:
#            sage: R = RealIntervalField()
#            sage: R(6).gamma()
#            120.000000000000
#            sage: R(1.5).gamma()
#            0.886226925452757
#         """
#         cdef RealIntervalFieldElement x
#         x = self._new()
#         _sig_on
#         mpfi_gamma(x.value, self.value)
#         _sig_off
#         return x

#     def zeta(self):
#         r"""
#         Return the Riemann zeta function evaluated at this real number.

#         \note{PARI is vastly more efficient at computing the Riemann zeta
#         function.   See the example below for how to use it.}

#         EXAMPLES:
#             sage: R = RealIntervalField()
#             sage: R(2).zeta()
#             1.64493406684822
#             sage: R.pi()^2/6
#             1.64493406684822
#             sage: R(-2).zeta()
#             0.000000000000000
#             sage: R(1).zeta()
#             +infinity

#         Computing zeta using PARI is much more efficient in difficult cases.
#         Here's how to compute zeta with at least a given precision:

#              sage: z = pari.new_with_bits_prec(2, 53).zeta(); z
#              1.644934066848226436472415167              # 32-bit
#              1.6449340668482264364724151666460251892    # 64-bit

#         Note that the number of bits of precision in the constructor only
#         effects the internel precision of the pari number, not the number
#         of digits that gets displayed.  To increase that you must
#         use \code{pari.set_real_precision}.

#              sage: type(z)
#              <type 'sage.libs.pari.gen.gen'>
#              sage: R(z)
#              1.64493406684822
#         """
#         cdef RealIntervalFieldElement x
#         x = self._new()
#         _sig_on
#         mpfi_zeta(x.value, self.value)
#         _sig_off
#         return x

    # Removed algdep/algebraic_dependency...these don't really make sense
    # for intervals

RR = RealIntervalField()


def create_RealIntervalFieldElement(s, int base=10, int pad=0, min_prec=53):
    r"""
    Return the real number defined by the string s as an element of
    \code{RealIntervalField(prec=n)}, where n potentially has slightly more
    (controlled by pad) bits than given by s.

    INPUT:
        s -- a string that defines a real number (or something whose
             string representation defines a number)
        base -- an integer between 2 and 36
        pad -- an integer >= 0.
        min_prec -- number will have at least this many bits of precision, no matter what.

    EXAMPLES:
        sage: RealIntervalFieldElement('2.3')
        [2.2999999999999998 ... 2.3000000000000003]
        sage: RealIntervalFieldElement(10)
        [10.000000000000000 ... 10.000000000000000]
        sage: RealIntervalFieldElement('1.0000000000000000000000000000000000')
        [1.000000000000000000000000000000000000 ... 1.000000000000000000000000000000000000]
    """
    if not isinstance(s, str):
        s = str(s)
    if base == 10:
        bits = int(3.32192*len(s))
    else:
        bits = int(math.log(base,2)*len(s))
    R = RealIntervalField(prec=max(bits+pad, min_prec))
    return RealIntervalFieldElement(R, s, base)

