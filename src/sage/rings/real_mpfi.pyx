"""
Field of Arbitrary Precision Real Intervals

AUTHORS:
    -- Carl Witty (2007-01-21): based on \code{real_mpfr.pyx};
               changed it to use mpfi rather than mpfr.
    -- William Stein (2007-01-24): modifications and clean up and docs, etc.

This is a straightforward binding to the MPFI library; it may be useful
to refer to its documentation for more details.

An interval is represented as a pair of floating-point numbers a and b
(where a <= b) and printed as [a .. b] (for instance,
[3.1415 .. 3.1416]).  These floating-point numbers implemented using
MPFR (the same as the RealNumber elements of RealField).

The interval represents the set { x : a <= x <= b } (so if a == b, then
the interval represents that particular floating-point number).  The
endpoints can include positive and negative infinity, with the obvious
meaning.  It is also possible to have a NaN (not-a-number) interval, which is
represented by having either endpoint be NaN.

PRINTING:

Intervals are printed with the left value rounded down and
the right rounded up, which is conservative, but in some
ways unsatisfying.

Consider a 3-bit interval containing exactly the floating-point number
1.25.  In round-to-nearest or round-down, this prints as 1.2; in
round-up, this prints as 1.3.  The straightforward options, then, are
to print this interval as [1.2 .. 1.2] (which does not even contain
the true value, 1.25), or to print it as [1.2 .. 1.3] (which gives
the impression that the upper and lower bounds are not equal, even
though they really are).  Neither of these is very satisfying, but I
have chosen the latter for now.

EXAMPLES:
    sage: R = RealIntervalField(3)
    sage: a = R(1.25)
    sage: a
    [1.2 .. 1.3]
    sage: a == 1.25
    True
    sage: a == 2
    False

COMPARISONS:

Comparison operations (==,!=,<,<=,>,>=) return true if every value in
the first interval has the given relation to every value in the second
interval.  The cmp(a, b) function works differently; it compares two
intervals lexicographically.  (However, the behavior is not specified
if given a non-interval and an interval.)

This convention for comparison operators has good and bad points.  The
good:
* Expected transitivity properties hold (if a > b and b == c, then a > c;
etc.)
* if a>b, then cmp(a, b) == 1; if a==b, then cmp(a,b) == 0; if a<b, then
cmp(a, b) == -1
* a==0 is true if the interval contains only the floating-point number 0;
similarly for a==1
* a>0 means something useful (that every value in the interval is greater
than 0)

The bad:
* Trichotomy fails to hold: there are values (a,b) such that none of a<b,
a==b, or a>b are true
* It is not the case that if cmp(a, b) == 0 then a==b, or that if
cmp(a, b) == 1 then a>b, or that if cmp(a, b) == -1 then a<b
* There are values a,b such that a<=b but neither a<b nor a==b hold

Note that intervals a and b overlap iff not(a != b).

EXAMPLES:
    sage: 0 < RIF(1, 2)
    True
    sage: 0 == RIF(0)
    True
    sage: not(0 == RIF(0, 1))
    True
    sage: not(0 != RIF(0, 1))
    True
    sage: 0 <= RIF(0, 1)
    True
    sage: not(0 < RIF(0, 1))
    True
    sage: cmp(RIF(0), RIF(0, 1))
    -1
    sage: cmp(RIF(0, 1), RIF(0))
    1
    sage: cmp(RIF(0, 1), RIF(1))
    -1
    sage: cmp(RIF(0, 1), RIF(0, 1))
    0
"""

############################################################################
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005-2006 William Stein <wstein@gmail.com>
#
#                  http://www.gnu.org/licenses/
############################################################################

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

import operator

from integer import Integer
from integer cimport Integer

from real_double import RealDoubleElement
from real_double cimport RealDoubleElement

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

    EXAMPLES: Creation of elements

    First with default precision.  First we coerce elements of various
    types, then we coerce intervals.

        sage: RIF = RealIntervalField(); RIF
        Real Interval Field with 53 bits of precision
        sage: RIF(3)
        [3.0000000000000000 .. 3.0000000000000000]
        sage: RIF(RIF(3))
        [3.0000000000000000 .. 3.0000000000000000]
        sage: RIF(pi)
        [3.1415926535897931 .. 3.1415926535897936]
        sage: RIF(RealField(53)('1.5'))
        [1.5000000000000000 .. 1.5000000000000000]
        sage: RIF(-2/19)
        [-0.10526315789473686 .. -0.10526315789473683]
        sage: RIF(-3939)
        [-3939.0000000000000 .. -3939.0000000000000]
        sage: RIF(-3939r)
        [-3939.0000000000000 .. -3939.0000000000000]
        sage: RIF('1.5')
        [1.5000000000000000 .. 1.5000000000000000]
        sage: RIF(RQDF.pi())
        [3.1415926535897926 .. 3.1415926535897941]
        sage: qdpi = RQDF(RealField(500).pi())
        sage: cmp(RealIntervalField(212)(qdpi), RealIntervalField(212).pi())
        0

    The base must be explicitly specified as a named parameter:
        sage: RIF('101101', base=2)
        [45.000000000000000 .. 45.000000000000000]
        sage: RIF('+infinity')
        [+infinity .. +infinity]
        sage: RIF('[1..3]')
        [1.0000000000000000 .. 3.0000000000000000]

    Next we coerce some 2-tuples, which define intervals.
        sage: RIF((-1.5, -1.3))
        [-1.5000000000000000 .. -1.3000000000000000]
        sage: RIF((RDF('-1.5'), RDF('-1.3')))
        [-1.5000000000000000 .. -1.3000000000000000]
        sage: RIF((1/3,2/3))
        [0.33333333333333331 .. 0.66666666666666675]

    The extra paranthesis aren't needed.
        sage: RIF(1/3,2/3)
        [0.33333333333333331 .. 0.66666666666666675]
        sage: RIF((1,2))
        [1.0000000000000000 .. 2.0000000000000000]
        sage: RIF((1r,2r))
        [1.0000000000000000 .. 2.0000000000000000]
        sage: RIF((pi, e))
        [2.7182818284590455 .. 3.1415926535897932]

    Values which can be represented as an exact floating-point number
    (of the precision of this RealIntervalField) result in a precise
    interval, where the lower bound is equal to the upper bound
    (even if they print differently).  Other values typically result
    in an interval where the lower and upper bounds are adjacent
    floating-point numbers.
        sage: def check(x):
        ...       return (x, x.lower() == x.upper())
        sage: check(RIF(pi))
        ([3.1415926535897931 .. 3.1415926535897936], False)
        sage: check(RIF(RR(pi)))
        ([3.1415926535897931 .. 3.1415926535897932], True)
        sage: check(RIF(1.5))
        ([1.5000000000000000 .. 1.5000000000000000], True)
        sage: check(RIF('1.5'))
        ([1.5000000000000000 .. 1.5000000000000000], True)
        sage: check(RIF(0.1))
        ([0.10000000000000000 .. 0.10000000000000001], True)
        sage: check(RIF(1/10))
        ([0.099999999999999991 .. 0.10000000000000001], False)
        sage: check(RIF('0.1'))
        ([0.099999999999999991 .. 0.10000000000000001], False)

    Similarly, when specifying both ends of an interval, the lower end
    is rounded down and the upper end is rounded up.
        sage: outward = RIF(1/10, 7/10); outward
        [0.099999999999999991 .. 0.70000000000000007]
        sage: nearest = RIF(RR(1/10), RR(7/10)); nearest
        [0.10000000000000000 .. 0.69999999999999996]
        sage: nearest.lower() - outward.lower()
        1.38777878078144e-17
        sage: outward.upper() - nearest.upper()
        1.11022302462516e-16

    Some examples with a real interval field of higher precision:
        sage: R = RealIntervalField(100)
        sage: R(3)
        [3.0000000000000000000000000000000 .. 3.0000000000000000000000000000000]
        sage: R(R(3))
        [3.0000000000000000000000000000000 .. 3.0000000000000000000000000000000]
        sage: R(pi)
        [3.1415926535897932384626433832793 .. 3.1415926535897932384626433832825]
        sage: R(-2/19)
        [-0.10526315789473684210526315789481 .. -0.10526315789473684210526315789470]
        sage: R(e,pi)
        [2.7182818284590452353602874713512 .. 3.1415926535897932384626433832825]
    """

    def __init__(self, int prec=53, int sci_not=0):
        if prec < MPFR_PREC_MIN or prec > MPFR_PREC_MAX:
            raise ValueError, "prec (=%s) must be >= %s and <= %s."%(
                prec, MPFR_PREC_MIN, MPFR_PREC_MAX)
        self.__prec = prec
        self.sci_not = sci_not
        self.__lower_field = RealField(prec, sci_not, "RNDD")
        self.__middle_field = RealField(prec, sci_not, "RNDN")
        self.__upper_field = RealField(prec, sci_not, "RNDU")
        ParentWithGens.__init__(self, self, tuple([]), False)

    def _lower_field(self):
        return self.__lower_field

    def _middle_field(self):
        return self.__middle_field

    def _upper_field(self):
        return self.__upper_field

    def _real_field(self, rnd):
        if rnd == "RNDD":
            return self._lower_field()
        elif rnd == "RNDN":
            return self._middle_field()
        elif rnd == "RNDU":
            return self._upper_field()
        else:
            return RealField(self.__prec, self.sci_not, "RNDZ")

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
        return "\\I \\R"

    def is_exact(self):
        return False

    def __call__(self, x, y=None, int base=10):
        """
        Create an element in this real interval field.

        INPUT:
            x -- a number, string, or 2-tuple.
            y -- (default: None); if given x is set to (x,y);
                 this is so you can write R(2,3) to make the interval from 2 to 3.
            base -- integer (default: 10) -- only used if x is a string.

        OUTPUT:
            an element of this real interval field.

        EXAMPLES:
            sage: R = RealIntervalField(20)
            sage: R('1.234')
            [1.2339992 .. 1.2340012]
            sage: R('2', base=2)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert number to real interval.
            sage: a = R('1.1001', base=2); a
            [1.5625000 .. 1.5625000]
            sage: a.str(2)
            '[1.1001000000000000000 .. 1.1001000000000000000]'

        Type:
            RealIntervalField?
        for more information.
        """
        if not y is None:
            x = (x,y)
        return RealIntervalFieldElement(self, x, base)

    def construction(self):
        """
        Returns the functorial construction of self, namely, completion of
        the rational numbers with respect to the prime at $\infinity$,
        and the note that this is an interval field.

        Also preserves other information that makes this field unique
        (e.g. precision, print mode).

        EXAMPLES:
            sage: R = RealIntervalField(123)
            sage: c, S = R.construction(); S
            Rational Field
            sage: R == c(S)
            True
        """
        from sage.categories.pushout import CompletionFunctor
        return (CompletionFunctor(sage.rings.infinity.Infinity,
                                  self.prec(),
                                  {'sci_not': self.scientific_notation(), 'type': 'Interval'}),
               sage.rings.rational_field.QQ)

    cdef _coerce_c_impl(self, x):
        """
        Canonical coercion of x to this mpfi real field.

        The rings that canonically coerce to this mpfi real field are:
             * this mpfi field itself
             * any mpfr real field with precision that is as large as this one
             * any other mpfi real field with precision that is as large as this one
             * anything that canonically coerces to the mpfr real
               field with same precision as self.

        Values which can be exactly represented as a floating-point
        number are coerced to a precise interval, with upper and lower
        bounds equal; otherwise, the upper and lower bounds will typically
        be adjacent floating-point numbers that surround the given value.
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
        cdef RealNumber lower, upper
        try:
            lower = self.__lower_field._coerce_(x)
            upper = self.__upper_field._coerce_(x)
            return self(lower, upper)
        except TypeError, msg:
            raise TypeError, "no canonical coercion of element into self"

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: RealIntervalField(10) == RealIntervalField(11)
            False
            sage: RealIntervalField(10) == RealIntervalField(10)
            True
            sage: RealIntervalField(10,sci_not=True) == RealIntervalField(10,sci_not=False)
            True
            sage: RealIntervalField(10) == IntegerRing()
            False
        """
        if not isinstance(other, RealIntervalField):
            return -1
        cdef RealIntervalField _other
        _other = other  # to access C structure
        if self.__prec == _other.__prec:
            return 0
        return 1

    def __reduce__(self):
        """
        EXAMPLES:
            sage: R = RealIntervalField(sci_not=1, prec=200)
            sage: loads(dumps(R)) == R
            True
        """
        return __create__RealIntervalField_version0, (self.__prec, self.sci_not)

    def gen(self, i=0):
        if i == 0:
            return self(1)
        else:
            raise IndexError

    def complex_field(self):
        """
        Return complex field of the same precision.

        EXAMPLES:
            sage: RIF.complex_field()
            Complex Interval Field with 53 bits of precision
        """
        return sage.rings.complex_interval_field.ComplexIntervalField(self.prec())

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

    def pi(self):
        """
        Returns pi to the precision of this field.

        EXAMPLES:
            sage: R = RealIntervalField(100)
            sage: R.pi()
            [3.1415926535897932384626433832793 .. 3.1415926535897932384626433832825]
            sage: R.pi().sqrt()/2
            [0.88622692545275801364908374166983 .. 0.88622692545275801364908374167142]
            sage: R = RealIntervalField(150)
            sage: R.pi().sqrt()/2
            [0.88622692545275801364908374167057259139877472759 .. 0.88622692545275801364908374167057259139877472830]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_const_pi(x.value)
        return x

    def euler_constant(self):
        """
        Returns Euler's gamma constant to the precision of this field.

        EXAMPLES:
            sage: RealIntervalField(100).euler_constant()
            [0.57721566490153286060651209008233 .. 0.57721566490153286060651209008313]
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
            [0.69314718055994530941723212145798 .. 0.69314718055994530941723212145878]
            sage: R(2).log()
            [0.69314718055994530941723212145798 .. 0.69314718055994530941723212145878]
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
            return self.sci_not
        else:
            self.sci_not = status

    def zeta(self, n=2):
        """
        Return an $n$-th root of unity in the real field,
        if one exists, or raise a ValueError otherwise.

        EXAMPLES:
            sage: R = RealIntervalField()
            sage: R.zeta()
            [-1.0000000000000000 .. -1.0000000000000000]
            sage: R.zeta(1)
            [1.0000000000000000 .. 1.0000000000000000]
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

    def __init__(self, RealIntervalField parent, x=0, int base=10):
        """
        Create a real interval element.  Should be called by first creating
        a RealIntervalField, as illustrated in the examples.

        EXAMPLES:
            sage: R = RealIntervalField()
            sage: R('1.2456')
            [1.2455999999999998 .. 1.2456000000000001]
            sage: R = RealIntervalField(3)
            sage: R('1.2456')
            [1.0 .. 1.3]

        EXAMPLE: Rounding
            sage: w = RealIntervalField(3)(5/2)
            sage: RealIntervalField(2)(w).str(2)
            '[10. .. 11.]'

        Type:
            RealIntervalField?
        for many more examples.
        """
        import sage.rings.qqbar

        self.init = 0
        if parent is None:
            raise TypeError
        self._parent = parent
        mpfi_init2(self.value, parent.__prec)
        self.init = 1
        if x is None: return
        cdef RealIntervalFieldElement _x, n, d
        cdef RealNumber rn, rn1
        cdef Rational rat, rat1
        cdef Integer integ, integ1
        cdef RealDoubleElement dx, dx1
        cdef int ix, ix1
        if PY_TYPE_CHECK(x, RealIntervalFieldElement):
            _x = x  # so we can get at x.value
            mpfi_set(self.value, _x.value)
        elif PY_TYPE_CHECK(x, RealNumber):
            rn = x
            mpfi_set_fr(self.value, <mpfr_t> rn.value)
        elif PY_TYPE_CHECK(x, Rational):
            rat = x
            mpfi_set_q(self.value, <mpq_t> rat.value)
        elif PY_TYPE_CHECK(x, Integer):
            integ = x
            mpfi_set_z(self.value, <mpz_t> integ.value)
        elif PY_TYPE_CHECK(x, int):
            ix = x
            mpfi_set_si(self.value, ix)
        elif isinstance(x, tuple):
            try:
                a, b = x
            except ValueError:
                raise TypeError, "tuple defining an interval must have length 2"
            if PY_TYPE_CHECK(a, RealNumber) and PY_TYPE_CHECK(b, RealNumber):
                rn = a
                rn1 = b
                mpfi_interv_fr(self.value, <mpfr_t> rn.value, <mpfr_t> rn1.value)
            elif PY_TYPE_CHECK(a, RealDoubleElement) and PY_TYPE_CHECK(b, RealDoubleElement):
                dx = a
                dx1 = b
                mpfi_interv_d(self.value, dx._value, dx1._value)
            elif PY_TYPE_CHECK(a, Rational) and PY_TYPE_CHECK(b, Rational):
                rat = a
                rat1 = b
                # todo: The <object> coerce is evidently to get around a weird bug in SageX (?)
                mpfi_interv_q(self.value, <object>rat.value, <object>rat1.value)
            elif PY_TYPE_CHECK(a, Integer) and PY_TYPE_CHECK(b, Integer):
                integ = a
                integ1 = b
                mpfi_interv_z(self.value, <object> integ.value, <object> integ1.value)
            elif PY_TYPE_CHECK(a, int) and PY_TYPE_CHECK(b, int):
                ix = a
                ix1 = b
                mpfi_interv_si(self.value, ix, ix1)
            else:  # generic fallback -- coerce both endpoints to reals.
                rn = self._parent._lower_field()(a)
                rn1 =self._parent._upper_field()(b)
                mpfi_interv_fr(self.value, <mpfr_t> rn.value, <mpfr_t> rn1.value)

        elif isinstance(x, sage.rings.qqbar.AlgebraicReal):
            d = x.interval(self._parent)
            mpfi_set(self.value, d.value)

        else:
            if isinstance(x, str):
                # string
                s = str(x).replace('..', ',').replace(' ','').replace('+infinity', '@inf@').replace('-infinity','-@inf@')
                if mpfi_set_str(self.value, s, base):
                    raise TypeError, "Unable to convert number to real interval."
            else:
                # try coercing to real
                try:
                    rn = self._parent._lower_field()(x)
                    rn1 = self._parent._upper_field()(x)
                except TypeError:
                    raise TypeError, "Unable to convert number to real interval."
                mpfi_interv_fr(self.value, <mpfr_t> rn.value, <mpfr_t> rn1.value)

    def __reduce__(self):
        """
        Pickling support.

        EXAMPLES:
            sage: a = RIF(5,5.5)
            sage: cmp(loads(dumps(a)), a)
            0
            sage: R = RealIntervalField(sci_not=1, prec=200)
            sage: b = R('393.39203845902384098234098230948209384028340')
            sage: cmp(loads(dumps(b)), b)
            0
            sage: b = R(1)/R(0); b
            [+infinity .. +infinity]
            sage: loads(dumps(b)) == b
            True
            sage: b = R(-1)/R(0); b
            [-infinity .. -infinity]
            sage: loads(dumps(b)) == b
            True
            sage: b = R('[2 .. 3]'); b
            [2.0000000000000000000000000000000000000000000000000000000000000e0 .. 3.0000000000000000000000000000000000000000000000000000000000000e0]
            sage: cmp(loads(dumps(b)), b)
            0
        """
        s = self.str(32, no_sci=False, e='@')
        return (__create__RealIntervalFieldElement_version0, (self._parent, s, 32))

    def  __dealloc__(self):
        if self.init:
            mpfi_clear(self.value)

    def __repr__(self):
        return self.str(10)

    def _latex_(self):
        return str(self)

    def _interface_init_(self):
        """
        Raise a TypeError.

        This function would return the string representation of self
        that makes sense as a default representation of a real
        interval in other computer algebra systems.  But, most
        other computer algebra systems do not support interval
        arithmetic, so instead we just raise a TypeError.

        Define the appropriate _cas_init_ function if there is a
        computer algebra system you would like to support.

        EXAMPLES:
            sage: n = RIF(1.3939494594)
            sage: n._interface_init_()
            Traceback (most recent call last):
            ...
            TypeError

        Here a conversion to Maxima happens implicitly, which
        results in a type error:
            sage: a = RealInterval('2.3')
            sage: erf(a)
            Traceback (most recent call last):
            ...
            TypeError
        """
        raise TypeError
        #return self.str(10, no_sci=True)

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
    # see [0.333 .. 0.333] for an interval with unequal left and right
    # sides.
    def str(self, int base=10, no_sci=None, e='e'):
        """
        INPUT:
             base -- base for output
             no_sci -- if True do not print using scientific notation; if False
                       print with scientific notation; if None (the default), print how the parent prints.
             e - symbol used in scientific notation

        EXAMPLES:
            sage: a = RIF(59/27); a
            [2.1851851851851851 .. 2.1851851851851856]
            sage: a.str(16)
            '[2.2f684bda12f68 .. 2.2f684bda12f6a]'
            sage: a.str(no_sci=False)
            '[2.1851851851851851e0 .. 2.1851851851851856e0]'
        """
        if base < 2 or base > 36:
            raise ValueError, "the base (=%s) must be between 2 and 36"%base
        if mpfi_nan_p(self.value):
            if base >= 24:
                return "[.. @NaN@ ..]"
            else:
                return "[.. NaN ..]"

        t1 = self.lower().str(base=base, no_sci=no_sci, e=e, truncate=False)
        t2 = self.upper().str(base=base, no_sci=no_sci, e=e, truncate=False)

        return "[%s .. %s]"%(t1, t2)

    def __copy__(self):
        """
        Return copy of self -- since self is immutable, we just return self again.

        EXAMPLES:
            sage: a = RIF(3.5)
            sage: copy(a) is  a
            True
        """
        return self

    # Interval-specific functions
    def lower(self, rnd=None):
        """
        Returns the lower bound of this interval

        rnd -- (string) the rounding mode
                RNDN -- round to nearest
                RNDD -- (default) round towards minus infinity
                RNDZ -- round towards zero
                RNDU -- round towards plus infinity

        The rounding mode does not affect the value returned as a
        floating-point number, but it does control which variety
        of RealField the returned number is in, which affects printing
        and subsequent operations.

        EXAMPLES:
            sage: R = RealIntervalField(13)
            sage: R.pi().lower().str(truncate=False)
            '3.1411'

            sage: x = R(1.2,1.3); x
            [1.1999 .. 1.3001]
            sage: x.lower()
            1.19
            sage: x.lower('RNDU')
            1.20
            sage: x.lower('RNDN')
            1.20
            sage: x.lower('RNDZ')
            1.19
            sage: x.lower().parent()
            Real Field with 13 bits of precision and rounding RNDD
            sage: x.lower('RNDU').parent()
            Real Field with 13 bits of precision and rounding RNDU
            sage: x.lower() == x.lower('RNDU')
            True
        """
        cdef RealNumber x
        if rnd is None:
            x = (<RealIntervalField>self._parent).__lower_field._new()
        else:
            x = (<RealField>(self._parent._real_field(rnd)))._new()
        mpfi_get_left(<mpfr_t> x.value, self.value)
        return x

    def upper(self, rnd =None):
        """
        Returns the upper bound of this interval

        rnd -- (string) the rounding mode
                RNDN -- round to nearest
                RNDD -- round towards minus infinity
                RNDZ -- round towards zero
                RNDU -- (default) round towards plus infinity

        The rounding mode does not affect the value returned as a
        floating-point number, but it does control which variety
        of RealField the returned number is in, which affects printing
        and subsequent operations.

        EXAMPLES:
            sage: R = RealIntervalField(13)
            sage: R.pi().upper().str(truncate=False)
            '3.1417'

            sage: R = RealIntervalField(13)
            sage: x = R(1.2,1.3); x
            [1.1999 .. 1.3001]
            sage: x.upper()
            1.31
            sage: x.upper('RNDU')
            1.31
            sage: x.upper('RNDN')
            1.30
            sage: x.upper('RNDD')
            1.30
            sage: x.upper('RNDZ')
            1.30
            sage: x.upper().parent()
            Real Field with 13 bits of precision and rounding RNDU
            sage: x.upper('RNDD').parent()
            Real Field with 13 bits of precision and rounding RNDD
            sage: x.upper() == x.upper('RNDD')
            True
        """
        cdef RealNumber x
        if rnd is None:
            x = (<RealIntervalField>self._parent).__upper_field._new()
        else:
            x = ((<RealField>self._parent._real_field(rnd)))._new()
        mpfi_get_right(<mpfr_t> x.value, self.value)
        return x

    def absolute_diameter(self):
        """
        The diameter of this interval (for [a .. b], this is b-a),
        rounded upward, as a RealNumber.

        EXAMPLES:
            sage: RIF(1, pi).absolute_diameter()
            2.14159265358979
        """
        cdef RealNumber x
        x = (<RealIntervalField>self._parent).__middle_field._new()
        mpfi_diam_abs(<mpfr_t> x.value, self.value)
        return x

    def relative_diameter(self):
        """
        The relative diameter of this interval (for [a .. b],
        this is (b-a)/((a+b)/2)), rounded upward, as a RealNumber.

        EXAMPLES:
            sage: RIF(1, pi).relative_diameter()
	    1.03418797197910
        """
        cdef RealNumber x
        x = (<RealIntervalField>self._parent).__middle_field._new()
        mpfi_diam_rel(<mpfr_t> x.value, self.value)
        return x

    def diameter(self):
        """
        If (0 in self), returns self.absolute_diameter(),
        otherwise self.relative_diameter().

        EXAMPLES:
            sage: RIF(1, 2).diameter()
            0.666666666666667
            sage: RIF(1, 2).absolute_diameter()
            1.00000000000000
            sage: RIF(1, 2).relative_diameter()
            0.666666666666667
            sage: RIF(pi).diameter()
            1.41357985842823e-16
            sage: RIF(pi).absolute_diameter()
            4.44089209850063e-16
            sage: RIF(pi).relative_diameter()
            1.41357985842823e-16
            sage: (RIF(pi) - RIF(3, 22/7)).diameter()
            0.142857142857144
            sage: (RIF(pi) - RIF(3, 22/7)).absolute_diameter()
            0.142857142857144
            sage: (RIF(pi) - RIF(3, 22/7)).relative_diameter()
            2.03604377705518
        """
        cdef RealNumber x
        x = (<RealIntervalField>self._parent).__middle_field._new()
        mpfi_diam(<mpfr_t> x.value, self.value)
        return x

    def is_exact(self):
        return mpfr_equal_p(&self.value.left, &self.value.right)

    def magnitude(self):
        """
        The largest absolute value of the elements of the interval.

        EXAMPLES:
            sage: RIF(-2, 1).magnitude()
            2.00000000000000
            sage: RIF(-1, 2).magnitude()
            2.00000000000000
        """
        cdef RealNumber x
        x = (<RealIntervalField>self._parent).__middle_field._new()
        mpfi_mag(<mpfr_t> x.value, self.value)
        return x

    def mignitude(self):
        """
        The smallest absolute value of the elements of the interval.

        EXAMPLES:
            sage: RIF(-2, 1).mignitude()
            0.000000000000000
            sage: RIF(-2, -1).mignitude()
            1.00000000000000
            sage: RIF(3, 4).mignitude()
            3.00000000000000
        """
        cdef RealNumber x
        x = (<RealIntervalField>self._parent).__middle_field._new()
        mpfi_mig(<mpfr_t> x.value, self.value)
        return x

    def center(self):
        """
        RIF(a, b).center() == (a+b)/2

        EXAMPLES:
            sage: RIF(1, 2).center()
            1.50000000000000
        """
        cdef RealNumber x
        x = (<RealIntervalField>self._parent).__middle_field._new()
        mpfi_mid(<mpfr_t> x.value, self.value)
        return x

    def alea(self):
        """
        RIF(a, b).alea() gives a floating-point number picked at
        random from the interval.

        EXAMPLES:
            sage: RIF(1, 2).alea() # random
            1.34696133696137
        """
        cdef RealNumber x
        x = (<RealIntervalField>self._parent).__middle_field._new()
        mpfi_alea(<mpfr_t> x.value, self.value)
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
        Add two real intervals with the same parent.

        EXAMPLES:
            sage: R = RealIntervalField()
            sage: R(-1.5) + R(2.5)
            [1.0000000000000000 .. 1.0000000000000000]
            sage: R('-1.3') + R('2.3')
            [0.99999999999999977 .. 1.0000000000000005]
            sage: R(1, 2) + R(3, 4)
            [4.0000000000000000 .. 6.0000000000000000]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_add(x.value, self.value, (<RealIntervalFieldElement>other).value)
        return x

    def __invert__(self):
        """
        Return the multiplicative "inverse" of this interval.
        (Technically, non-precise intervals don't have multiplicative
        inverses.)

        EXAMPLES:
            sage: v = RIF(2); v
            [2.0000000000000000 .. 2.0000000000000000]
            sage: ~v
            [0.50000000000000000 .. 0.50000000000000000]
            sage: v * ~v
            [1.0000000000000000 .. 1.0000000000000000]
            sage: v = RIF(1.5, 2.5); v
            [1.5000000000000000 .. 2.5000000000000000]
            sage: ~v
            [0.39999999999999996 .. 0.66666666666666675]
            sage: v * ~v
            [0.59999999999999986 .. 1.6666666666666670]
            sage: ~RIF(-1, 1)
            [-infinity .. +infinity]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_inv(x.value, self.value)
        return x

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        Subtract two real intervals with the same parent.

        EXAMPLES:
            sage: R = RealIntervalField()
            sage: R(-1.5) - R(2.5)
            [-4.0000000000000000 .. -4.0000000000000000]
            sage: R('-1.3') - R('2.7')
            [-4.0000000000000009 .. -3.9999999999999995]
            sage: R(1, 2) - R(3, 4)
            [-3.0000000000000000 .. -1.0000000000000000]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_sub(x.value, self.value, (<RealIntervalFieldElement>right).value)
        return x

    cdef RingElement _mul_c_impl(self, RingElement right):
        """
        Multiply two real intervals with the same parent.

        EXAMPLES:
            sage: R = RealIntervalField()
            sage: R(-1.5) * R(2.5)
            [-3.7500000000000000 .. -3.7500000000000000]
            sage: R('-1.3') * R('2.3')
            [-2.9900000000000007 .. -2.9899999999999993]
            sage: R(1, 2) * R(3, 4)
            [3.0000000000000000 .. 8.0000000000000000]

        If two elements have different precision, arithmetic
        operations are performed after coercing to the lower
        precision.

            sage: R20 = RealIntervalField(20)
            sage: R100 = RealIntervalField(100)
            sage: a = R20('393.3902834028345')
            sage: b = R100('393.3902834028345')
            sage: a
            [393.39013 .. 393.39063]
            sage: b
            [393.39028340283449999999999999970 .. 393.39028340283450000000000000011]
            sage: a*b
            [154755.75 .. 154756.25]
            sage: b*a
            [154755.75 .. 154756.25]
            sage: parent(b*a)
            Real Interval Field with 20 bits of precision
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_mul(x.value, self.value, (<RealIntervalFieldElement>right).value)
        return x


    cdef RingElement _div_c_impl(self, RingElement right):
        """
        Divide self by other, where both are real intervals with the same parent.

        EXAMPLES:
            sage: R = RealIntervalField()
            sage: R(1)/R(3)
            [0.33333333333333331 .. 0.33333333333333338]
            sage: R(1)/R(0)
            [+infinity .. +infinity]
            sage: R(1)/R(-1, 1)
            [-infinity .. +infinity]

            sage: R(-1.5) / R(2.5)
            [-0.60000000000000009 .. -0.59999999999999997]
            sage: R(1, 2) / R(3, 4)
            [0.25000000000000000 .. 0.66666666666666675]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_div((<RealIntervalFieldElement>x).value, self.value,
                 (<RealIntervalFieldElement>right).value)
        return x

    cdef ModuleElement _neg_c_impl(self):
        """
        Return the additive "inverse" of this interval.
        (Technically, non-precise intervals don't have additive inverses.)

        EXAMPLES:
            sage: v = RIF(2); v
            [2.0000000000000000 .. 2.0000000000000000]
            sage: -v
            [-2.0000000000000000 .. -2.0000000000000000]
            sage: v + -v
            [-0.00000000000000000 .. 0.00000000000000000]
            sage: v = RIF(1.5, 2.5); v
            [1.5000000000000000 .. 2.5000000000000000]
            sage: -v
            [-2.5000000000000000 .. -1.5000000000000000]
            sage: v + -v
            [-1.0000000000000000 .. 1.0000000000000000]
        """

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

    def square(self):
        """
        Return the square of the interval.  Note that squaring
        an interval is different than multiplying it by itself,
        because the square can never be negative.

        EXAMPLES:
            sage: RIF(1, 2).square()
            [1.0000000000000000 .. 4.0000000000000000]
            sage: RIF(-1, 1).square()
            [0.00000000000000000 .. 1.0000000000000000]
            sage: RIF(-1, 1) * RIF(-1, 1)
            [-1.0000000000000000 .. 1.0000000000000000]
        """

        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_sqr(x.value, self.value)
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
        Returns $x * 2^y$, for $y$ an integer.  Much faster than an ordinary
        multiplication.

        EXAMPLES:
            sage: RIF(1.0) << 32
            [4.2949672960000000e9 .. 4.2949672960000000e9]
        """
        if isinstance(x, RealIntervalFieldElement) and isinstance(y, (int,long, Integer)):
            return x._lshift_(y)
        return sage.structure.element.bin_op(x, y, operator.lshift)

    def _rshift_(RealIntervalFieldElement self, n):
        if n > sys.maxint:
            raise OverflowError, "n (=%s) must be <= %s"%(n, sys.maxint)
        cdef RealIntervalFieldElement x
        x = self._new()
        mpfi_div_2exp(x.value, self.value, n)
        return x

    def __rshift__(x, y):
        """
        Returns $x / 2^y$, for $y$ an integer.  Much faster than an ordinary
        division.

        EXAMPLES:
            sage: RIF(1024.0) >> 14
            [0.062500000000000000 .. 0.062500000000000000]
        """
        if isinstance(x, RealIntervalFieldElement) and \
               isinstance(y, (int,long,Integer)):
            return x._rshift_(y)
        return sage.structure.element.bin_op(x, y, operator.rshift)

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

    def floor(self):
        """
        Returns the floor of this number

        EXAMPLES:
            sage: R = RealIntervalField()
            sage: (2.99).floor()
            2
            sage: (2.00).floor()
            2
            sage: floor(RR(-5/2))
            -3
            sage: R = RealIntervalField(100)
            sage: a = R(9.5, 11.3); a
            [9.5000000000000000000000000000000 .. 11.300000000000000710542735760101]
            sage: floor(a)
            [9.0000000000000000000000000000000 .. 11.000000000000000000000000000000]
            sage: a.floor()
            [9.0000000000000000000000000000000 .. 11.000000000000000000000000000000]
            sage: ceil(a)
            [10.000000000000000000000000000000 .. 12.000000000000000000000000000000]
            sage: a.ceil()
            [10.000000000000000000000000000000 .. 12.000000000000000000000000000000]
        """
        return self.parent()(self.lower().floor(), self.upper().floor())

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
            sage: R = RealIntervalField(30)
            sage: a = R(-9.5, -11.3); a
            [-11.300000012 .. -9.5000000000]
            sage: a.floor()
            [-12.000000000 .. -10.000000000]
            sage: a.ceil()
            [-11.000000000 .. -9.0000000000]
            sage: ceil(a)
            [-11.000000000 .. -9.0000000000]
        """
        return self.parent()(self.lower().ceil(), self.upper().ceil())

    def ceiling(self):
         return self.ceil()

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

    def simplest_rational(self, low_open=False, high_open=False):
        """
        Returns the simplest rational in this interval.
        Given rationals a/b and c/d (both in lowest terms), the former
        is simpler if b<d or if b==d and |a|<|c|.

        If optional parameters low_open or high_open are true, then
        treat this as an open interval on that end.

        EXAMPLES:
            sage: RealIntervalField(10)(pi).simplest_rational()
            22/7
            sage: RealIntervalField(20)(pi).simplest_rational()
            355/113
            sage: RIF(0.123, 0.567).simplest_rational()
            1/2
            sage: RIF(RR(1/3).nextabove(), RR(3/7)).simplest_rational()
            2/5
            sage: RIF(1234/567).simplest_rational()
            1234/567
            sage: RIF(-8765/432).simplest_rational()
            -8765/432
            sage: RIF(-1.234, 0.003).simplest_rational()
            0
            sage: RIF(RR(1/3)).simplest_rational()
            6004799503160661/18014398509481984
            sage: RIF(RR(1/3)).simplest_rational(high_open=True)
            Traceback (most recent call last):
            ...
            ValueError: simplest_rational() on open, empty interval
            sage: RIF(1/3, 1/2).simplest_rational()
            1/2
            sage: RIF(1/3, 1/2).simplest_rational(high_open=True)
            1/3
            sage: phi = ((RealIntervalField(500)(5).sqrt() + 1)/2)
            sage: phi.simplest_rational() == fibonacci(362)/fibonacci(361)
            True

        """
        if mpfr_equal_p(&self.value.left, &self.value.right):
            if low_open or high_open:
                raise ValueError, 'simplest_rational() on open, empty interval'
            return self.lower().exact_rational()

        if mpfi_has_zero(self.value):
            return Rational(0)

        if mpfi_is_neg(self.value):
            return -(self._neg_c_impl().simplest_rational(low_open=high_open, high_open=low_open))

        low = self.lower()
        high = self.upper()

        # First, we try using approximate arithmetic of slightly higher
        # precision.
        cdef RealIntervalFieldElement highprec
        highprec = RealIntervalField(int(self.prec() * 1.2))(self)

        cdef Rational try1 = highprec._simplest_rational_helper()

        # Note that to compute "try1 >= low", SAGE converts try1 to a
        # floating-point number rounding down, and "try1 <= high"
        # rounds up (since "low" and "high" are in downward-rounding
        # and upward-rounding fields, respectively).
        if try1 >= low and try1 <= high:
            ok = True
            if low_open and (try1 == low.exact_rational()):
                ok = False
            if high_open and (try1 == high.exact_rational()):
                ok = False
            if ok:
                return try1

        # We could try again with higher precision; instead, we
        # go directly to using exact arithmetic.
        return _simplest_rational_exact(low.exact_rational(),
                                        high.exact_rational(),
                                        low_open,
                                        high_open)

    cdef Rational _simplest_rational_helper(self):
        """
        Returns the simplest rational in an interval which is
        either equal to or slightly larger than self.  We assume
        that both endpoints of self are nonnegative.
        """

        low = self.lower()

        cdef RealIntervalFieldElement new_elt

        if low <= 1:
            if low == 0:
                return Rational(0)
            if self.upper() >= 1:
                return Rational(1)
            new_elt = ~self
            return ~(new_elt._simplest_rational_helper())

        fl = low.floor()
        new_elt = self - fl
        return fl + new_elt._simplest_rational_helper()

    ###########################################
    # Comparisons: ==, !=, <, <=, >, >=
    ###########################################

    def is_NaN(self):
        return mpfi_nan_p(self.value)

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef _richcmp_c_impl(left, Element right, int op):
        """
        Implements comparisons between intervals.  (See the file header
        comment for more information on interval comparison.)

        EXAMPLES:
            sage: 0 < RIF(1, 3)
            True
            sage: 1 < RIF(1, 3)
            False
            sage: 2 < RIF(1, 3)
            False
            sage: 4 < RIF(1, 3)
            False
            sage: RIF(0, 1/2) < RIF(1, 3)
            True
            sage: RIF(0, 1) < RIF(1, 3)
            False
            sage: RIF(0, 2) < RIF(1, 3)
            False
            sage: RIF(1, 2) < RIF(1, 3)
            False
            sage: RIF(1, 3) < 4
            True
            sage: RIF(1, 3) < 3
            False
            sage: RIF(1, 3) < 2
            False
            sage: RIF(1, 3) < 0
            False
            sage: 0 <= RIF(1, 3)
            True
            sage: 1 <= RIF(1, 3)
            True
            sage: 2 <= RIF(1, 3)
            False
            sage: 4 <= RIF(1, 3)
            False
            sage: RIF(0, 1/2) <= RIF(1, 3)
            True
            sage: RIF(0, 1) <= RIF(1, 3)
            True
            sage: RIF(0, 2) <= RIF(1, 3)
            False
            sage: RIF(1, 2) <= RIF(1, 3)
            False
            sage: RIF(1, 3) <= 4
            True
            sage: RIF(1, 3) <= 3
            True
            sage: RIF(1, 3) <= 2
            False
            sage: RIF(1, 3) <= 0
            False
            sage: RIF(1, 3) > 0
            True
            sage: RIF(1, 3) > 1
            False
            sage: RIF(1, 3) > 2
            False
            sage: RIF(1, 3) > 4
            False
            sage: RIF(1, 3) > RIF(0, 1/2)
            True
            sage: RIF(1, 3) > RIF(0, 1)
            False
            sage: RIF(1, 3) > RIF(0, 2)
            False
            sage: RIF(1, 3) > RIF(1, 2)
            False
            sage: 4 > RIF(1, 3)
            True
            sage: 3 > RIF(1, 3)
            False
            sage: 2 > RIF(1, 3)
            False
            sage: 0 > RIF(1, 3)
            False
            sage: RIF(1, 3) >= 0
            True
            sage: RIF(1, 3) >= 1
            True
            sage: RIF(1, 3) >= 2
            False
            sage: RIF(1, 3) >= 4
            False
            sage: RIF(1, 3) >= RIF(0, 1/2)
            True
            sage: RIF(1, 3) >= RIF(0, 1)
            True
            sage: RIF(1, 3) >= RIF(0, 2)
            False
            sage: RIF(1, 3) >= RIF(1, 2)
            False
            sage: 4 >= RIF(1, 3)
            True
            sage: 3 >= RIF(1, 3)
            True
            sage: 2 >= RIF(1, 3)
            False
            sage: 0 >= RIF(1, 3)
            False
            sage: 0 == RIF(0)
            True
            sage: 0 == RIF(1)
            False
            sage: 1 == RIF(0)
            False
            sage: 0 == RIF(0, 1)
            False
            sage: 1 == RIF(0, 1)
            False
            sage: RIF(0, 1) == RIF(0, 1)
            False
            sage: RIF(1) == 0
            False
            sage: RIF(1) == 1
            True
            sage: RIF(0) == RIF(0)
            True
            sage: RIF(pi) == RIF(pi)
            False
            sage: RIF(0, 1) == RIF(1, 2)
            False
            sage: RIF(1, 2) == RIF(0, 1)
            False
            sage: 0 != RIF(0)
            False
            sage: 0 != RIF(1)
            True
            sage: 1 != RIF(0)
            True
            sage: 0 != RIF(0, 1)
            False
            sage: 1 != RIF(0, 1)
            False
            sage: RIF(0, 1) != RIF(0, 1)
            False
            sage: RIF(1) != 0
            True
            sage: RIF(1) != 1
            False
            sage: RIF(0) != RIF(0)
            False
            sage: RIF(pi) != RIF(pi)
            False
            sage: RIF(0, 1) != RIF(1, 2)
            False
            sage: RIF(1, 2) != RIF(0, 1)
            False
        """
        cdef RealIntervalFieldElement lt, rt

        lt = left
        rt = right

        if op == 0: #<
            return mpfr_less_p(&lt.value.right, &rt.value.left)
        elif op == 2: #==
            # a == b iff a<=b and b <= a
            # (this gives a result with two comparisons, where the
            # obvious approach would use three)
            return mpfr_lessequal_p(&lt.value.right, &rt.value.left) \
                and mpfr_lessequal_p(&rt.value.right, &lt.value.left)
        elif op == 4: #>
            return mpfr_less_p(&rt.value.right, &lt.value.left)
        elif op == 1: #<=
            return mpfr_lessequal_p(&lt.value.right, &rt.value.left)
        elif op == 3: #!=
            return mpfr_less_p(&lt.value.right, &rt.value.left) \
                or mpfr_less_p(&rt.value.right, &lt.value.left)
        elif op == 5: #>=
            return mpfr_lessequal_p(&rt.value.right, &lt.value.left)

    def __nonzero__(self):
        """
        Returns true if self is not known to be exactly zero.

        EXAMPLES:
            sage: RIF(0).__nonzero__()
            False
            sage: RIF(1).__nonzero__()
            True
            sage: RIF(1, 2).__nonzero__()
            True
            sage: RIF(0, 1).__nonzero__()
            True
            sage: RIF(-1, 1).__nonzero__()
            True
        """
        return not (mpfr_zero_p(&self.value.left) and mpfr_zero_p(&self.value.right))

    def __cmp__(left, right):
        """
        Compare two intervals lexicographically.  Returns 0 if they
        are the same interval, -1 if the second is larger, or 1 if
        the first is larger.

        EXAMPLES:
            sage: cmp(RIF(0), RIF(1))
            -1
            sage: cmp(RIF(0, 1), RIF(1))
            -1
            sage: cmp(RIF(0, 1), RIF(1, 2))
            -1
            sage: cmp(RIF(0, 0.99999), RIF(1, 2))
            -1
            sage: cmp(RIF(1, 2), RIF(0, 1))
            1
            sage: cmp(RIF(1, 2), RIF(0))
            1
            sage: cmp(RIF(0, 1), RIF(0, 2))
            -1
            sage: cmp(RIF(0, 1), RIF(0, 1))
            0
            sage: cmp(RIF(0, 1), RIF(0, 1/2))
            1
        """
        return (<Element>left)._cmp(right)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Implements the lexicographic total order on intervals.
        """
        cdef RealIntervalFieldElement lt, rt

        lt = left
        rt = right

        cdef int i
        i = mpfr_cmp(&lt.value.left, &rt.value.left)
        if i < 0:
            return -1
        elif i > 0:
            return 1
        i = mpfr_cmp(&lt.value.right, &rt.value.right)
        if i < 0:
            return -1
        elif i > 0:
            return 1
        else:
            return 0

    def __contains__(self, other):
        """
        Test whether one interval (or real number) is totally contained
        in another.

        EXAMPLES:
            sage: RIF(0, 2) in RIF(1, 3)
            False
            sage: RIF(0, 2) in RIF(0, 2)
            True
            sage: RIF(1, 2) in RIF(0, 3)
            True
            sage: 1.0 in RIF(0, 2)
            True
            sage: pi in RIF(3.1415, 3.1416)
            True
            sage: 22/7 in RIF(3.1415, 3.1416)
            False
        """

        cdef RealIntervalFieldElement other_intv
        cdef RealNumber other_rn
        if PY_TYPE_CHECK(other, RealIntervalFieldElement):
            other_intv = other
            return mpfi_is_inside(other_intv.value, self.value)
        elif PY_TYPE_CHECK(other, RealNumber):
            other_rn = other
            return mpfi_is_inside_fr(<mpfr_t> other_rn.value, self.value)
        try:
            other_intv = self._parent(other)
            return mpfi_is_inside(other_intv.value, self.value)
        except TypeError, msg:
            return False

    def contains_zero(self):
        """
        Returns True if self is an interval containing zero.

        EXAMPLES:
            sage: RIF(0).contains_zero()
            True
            sage: RIF(1, 2).contains_zero()
            False
            sage: RIF(-1, 1).contains_zero()
            True
            sage: RIF(-1, 0).contains_zero()
            True
        """
        return mpfi_has_zero(self.value)

    def overlaps(self, RealIntervalFieldElement other):
        """
        Returns true if self and other are intervals with at least
        one value in common.  For intervals a and b, we have
        a.overlaps(b) iff not(a!=b).

        EXAMPLES:
            sage: RIF(0, 1).overlaps(RIF(1, 2))
            True
            sage: RIF(1, 2).overlaps(RIF(0, 1))
            True
            sage: RIF(0, 1).overlaps(RIF(2, 3))
            False
            sage: RIF(2, 3).overlaps(RIF(0, 1))
            False
            sage: RIF(0, 3).overlaps(RIF(1, 2))
            True
            sage: RIF(0, 2).overlaps(RIF(1, 3))
            True
        """
        return mpfr_greaterequal_p(&self.value.right, &other.value.left) \
           and mpfr_greaterequal_p(&other.value.right, &self.value.left)

    def intersection(self, other):
        """
        Return the intersection of two intervals.  If the intervals
        do not overlap, raises a ValueError.

        EXAMPLES:
            sage: RIF(1, 2).intersection(RIF(1.5, 3))
            [1.5000000000000000 .. 2.0000000000000000]
            sage: RIF(1, 2).intersection(RIF(4/3, 5/3))
            [1.3333333333333332 .. 1.6666666666666668]
            sage: RIF(1, 2).intersection(RIF(3, 4))
            Traceback (most recent call last):
            ...
            ValueError: intersection of non-overlapping intervals
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        cdef RealIntervalFieldElement other_intv
        if PY_TYPE_CHECK(other, RealIntervalFieldElement):
            other_intv = other
        else:
            # Let type errors from _coerce_ propagate...
            other_intv = self._parent(other)

        mpfi_intersect(x.value, self.value, other_intv.value)
        if mpfr_less_p(&x.value.right, &x.value.left):
            raise ValueError, "intersection of non-overlapping intervals"
        return x

    def union(self, other):
        """
        Return the union of two intervals, or of an interval and
        a real number (more precisely, the convex hull).

        EXAMPLES:
            sage: RIF(1, 2).union(RIF(pi, 22/7))
            [1.0000000000000000 .. 3.1428571428571433]
            sage: RIF(1, 2).union(pi)
            [1.0000000000000000 .. 3.1415926535897936]
            sage: RIF(1).union(RIF(0, 2))
            [0.00000000000000000 .. 2.0000000000000000]
            sage: RIF(1).union(RIF(-1))
            [-1.0000000000000000 .. 1.0000000000000000]
        """

        cdef RealIntervalFieldElement x
        x = self._new()
        cdef RealIntervalFieldElement other_intv
        cdef RealNumber other_rn
        if PY_TYPE_CHECK(other, RealIntervalFieldElement):
            other_intv = other
            mpfi_union(x.value, self.value, other_intv.value)
        elif PY_TYPE_CHECK(other, RealNumber):
            other_rn = other
            mpfi_set(x.value, self.value)
            mpfi_put_fr(x.value, <mpfr_t> other_rn.value)
        else:
            # Let type errors from _coerce_ propagate...
            other_intv = self._parent(other)
            mpfi_union(x.value, self.value, other_intv.value)
        return x

    ############################
    # Special Functions
    ############################


    def sqrt(self):
        """
        Return a square root of self.  Raises an error if self is
        nonpositive.

        If you use self.square_root() then an interval will always
        be returned (though it will be NaN if self is nonpositive).

        EXAMPLES:
            sage: r = RIF(4.0)
            sage: r.sqrt()
            [2.0000000000000000 .. 2.0000000000000000]
            sage: r.sqrt()^2 == r
            True

            sage: r = RIF(4344)
            sage: r.sqrt()
            [65.909028213136323 .. 65.909028213136339]
            sage: r.sqrt()^2 == r
            False
            sage: r in r.sqrt()^2
            True
            sage: r.sqrt()^2 - r
            [-9.0949470177292824e-13 .. 1.8189894035458565e-12]

            sage: r = RIF(-2.0)
            sage: r.sqrt()
            Traceback (most recent call last):
            ...
            ValueError: self (=[-2.0000000000000000 .. -2.0000000000000000]) is not >= 0

            sage: r = RIF(-2, 2)
            sage: r.sqrt()
            Traceback (most recent call last):
            ...
            ValueError: self (=[-2.0000000000000000 .. 2.0000000000000000]) is not >= 0
            """
        if self.lower() < 0:
            raise ValueError, "self (=%s) is not >= 0"%self
        return self.square_root()


    def square_root(self):
        """
        Return a square root of self.  An interval will always be
        returned (though it will be NaN if self is nonpositive).

        EXAMPLES:
            sage: r = RIF(-2.0)
            sage: r.square_root()
            [.. NaN ..]
            sage: r.sqrt()
            Traceback (most recent call last):
            ...
            ValueError: self (=[-2.0000000000000000 .. -2.0000000000000000]) is not >= 0
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

    def __pow__(self, exponent, modulus):
        if isinstance(exponent, (int, long, Integer)):
            return sage.rings.ring_element.RingElement.__pow__(self, exponent)
        return (self.log() * exponent).exp()


    def log(self, base='e'):
        """
        EXAMPLES:
            sage: R = RealIntervalField()
            sage: R(2).log()
            [0.69314718055994528 .. 0.69314718055994540]
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
            sage: r = RIF(16.0)
            sage: r.log2()
            [4.0000000000000000 .. 4.0000000000000000]

            sage: r = RIF(31.9); r.log2()
            [4.9954845188775065 .. 4.9954845188775075]

            sage: r = RIF(0.0, 2.0)
            sage: r.log2()
            [-infinity .. 1.0000000000000000]
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
            sage: r = RIF(16.0); r.log10()
            [1.2041199826559245 ... 1.2041199826559248]
            sage: r.log() / log(10.0)
            [1.2041199826559245 ... 1.2041199826559251]

            sage: r = RIF(39.9); r.log10()
            [1.6009728956867481 .. 1.6009728956867484]

            sage: r = RIF(0.0)
            sage: r.log10()
            [-infinity .. -infinity]

            sage: r = RIF(-1.0)
            sage: r.log10()
            [.. NaN ..]

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
            sage: r = RIF(0.0)
            sage: r.exp()
            [1.0000000000000000 .. 1.0000000000000000]

            sage: r = RIF(32.3)
            sage: a = r.exp(); a
            [1.0658884727486446e14 .. 1.0658884727486449e14]
            sage: a.log()
            [32.299999999999990 .. 32.300000000000005]

            sage: r = RIF(-32.3)
            sage: r.exp()
            [9.3818445884986834e-15 .. 9.3818445884986851e-15]
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
            sage: r = RIF(0.0)
            sage: r.exp2()
            [1.0000000000000000 .. 1.0000000000000000]

            sage: r = RIF(32.0)
            sage: r.exp2()
            [4.2949672960000000e9 .. 4.2949672960000000e9]

            sage: r = RIF(-32.3)
            sage: r.exp2()
            [1.8911724825302069e-10 .. 1.8911724825302073e-10]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_exp2(x.value, self.value)
        _sig_off
        return x

    def is_int(self):
        """
        OUTPUT:
            bool -- True or False
            n -- an integer

        EXAMPLES:
            sage: a = RIF(0.8,1.5)
            sage: a.is_int()
            (True, 1)
            sage: a = RIF(1.1,1.5)
            sage: a.is_int()
            (False, None)
            sage: a = RIF(1,2)
            sage: a.is_int()
            (False, None)
            sage: a = RIF(-1.1, -0.9)
            sage: a.is_int()
            (True, -1)
            sage: a = RIF(0.1, 1.9)
            sage: a.is_int()
            (True, 1)
        """
        a = (self.lower()).ceil()
        b = (self.upper()).floor()
        if a == b:
            return True, a
        else:
            return False, None

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
        Returns the cosine of this number.

        EXAMPLES:
            sage: t=RIF(pi)/2
            sage: t.cos()
            [-1.6081226496766367e-16 .. 6.1232339957367661e-17]
            sage: t.cos().cos()
            [0.99999999999999988 .. 1.0000000000000000]

        TESTS:
        This looped forever with an earlier version of MPFI, but now
        it works.
            sage: RIF(-1, 1).cos()
            [0.54030230586813965 .. 1.0000000000000000]
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
            [0.90929742682568169539601986591150 .. 0.90929742682568169539601986591230]
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
            sage: q = RIF.pi()/3
            sage: q.tan()
            [1.7320508075688767 .. 1.7320508075688779]
            sage: q = RIF.pi()/6
            sage: q.tan()
            [0.57735026918962562 .. 0.57735026918962585]
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
            sage: q = RIF.pi()/3; q
            [1.0471975511965976 .. 1.0471975511965979]
            sage: i = q.cos(); i
            [0.49999999999999988 .. 0.50000000000000012]
            sage: q2 = i.acos(); q2
            [1.0471975511965974 .. 1.0471975511965981]
            sage: q == q2
            False
            sage: q != q2
            False
            sage: q2.lower() == q.lower()
            False
            sage: q - q2
            [-4.4408920985006262e-16 .. 4.4408920985006262e-16]
            sage: q in q2
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
            sage: q = RIF.pi()/5; q
            [0.62831853071795862 .. 0.62831853071795874]
            sage: i = q.sin(); i
            [0.58778525229247302 .. 0.58778525229247325]
            sage: q2 = i.asin(); q2
            [0.62831853071795851 .. 0.62831853071795885]
            sage: q == q2
            False
            sage: q != q2
            False
            sage: q2.lower() == q.lower()
            False
            sage: q - q2
            [-2.2204460492503131e-16 .. 2.2204460492503131e-16]
            sage: q in q2
            True
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
            sage: q = RIF.pi()/5; q
            [0.62831853071795862 .. 0.62831853071795874]
            sage: i = q.tan(); i
            [0.72654252800536078 .. 0.72654252800536113]
            sage: q2 = i.atan(); q2
            [0.62831853071795851 .. 0.62831853071795885]
            sage: q == q2
            False
            sage: q != q2
            False
            sage: q2.lower() == q.lower()
            False
            sage: q - q2
            [-2.2204460492503131e-16 .. 2.2204460492503131e-16]
            sage: q in q2
            True
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
            sage: q = RIF.pi()/12
            sage: q.cosh()
            [1.0344656400955103 .. 1.0344656400955108]
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
            sage: q = RIF.pi()/12
            sage: q.sinh()
            [0.26480022760227067 .. 0.26480022760227079]
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
            sage: q = RIF.pi()/11
            sage: q.tanh()
            [0.27807942929585022 .. 0.27807942929585039]
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
            sage: q = RIF.pi()/2
            sage: i = q.acosh() ; i
            [1.0232274785475503 .. 1.0232274785475509]
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
            sage: q = RIF.pi()/7
            sage: i = q.sinh() ; i
            [0.46401763049299088 .. 0.46401763049299106]
            sage: i.asinh() - q
            [-1.6653345369377349e-16 .. 1.6653345369377349e-16]
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
            sage: q = RIF.pi()/7
            sage: i = q.tanh() ; i
            [0.42091124104853488 .. 0.42091124104853501]
            sage: i.atanh() - q
            [-1.6653345369377349e-16 .. 1.6653345369377349e-16]
        """
        cdef RealIntervalFieldElement x
        x = self._new()
        _sig_on
        mpfi_atanh(x.value, self.value)
        _sig_off
        return x

    def algdep(self, n):
        """
        Returns a polynomial of degree at most $n$ which is approximately
        satisfied by this number.  Note that the returned polynomial
        need not be irreducible, and indeed usually won't be if this number
        is a good approximation to an algebraic number of degree less than $n$.

        Pari needs to know the number of "known good bits" in the number;
        we automatically get that from the interval width.

        ALGORITHM: Uses the PARI C-library algdep command.


        EXAMPLE:
            sage: r = sqrt(RIF(2)); r
            [1.4142135623730949 .. 1.4142135623730952]
            sage: r.algdep(5)
            x^2 - 2

        If we compute a wrong, but precise, interval, we get a wrong answer.
            sage: r = sqrt(RealIntervalField(200)(2)) + (1/2)^40; r
            [1.4142135623740045435034616524476131176321718753769480731766796 .. 1.4142135623740045435034616524476131176321718753769480731766809]
            sage: r.algdep(5)
            7266488*x^5 + 22441629*x^4 - 90470501*x^3 + 23297703*x^2 + 45778664*x + 13681026

        But if we compute an interval that includes the number we mean,
        we're much more likely to get the right answer, even if the interval
        is very imprecise.
            sage: r = r.union(sqrt(2.0))
            sage: r.algdep(5)
            x^2 - 2

        Even on this extremely imprecise interval we get an answer
        which is technically correct.
            sage: RIF(-1, 1).algdep(5)
            x
        """

        # If 0 is in the interval, then we have no known bits!  But
        # fortunately, there's a perfectly valid answer we can
        # return anyway.
        if 0 in self:
            #import sage.rings.polynomial.polynomial_ring
            return sage.rings.polynomial.polynomial_ring.polygen(
                sage.rings.integer_ring.IntegerRing())

        known_bits = -self.relative_diameter().log2()

        return sage.rings.arith.algdep(self.center(), n, known_bits=known_bits)

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


def _simplest_rational_test_helper(low, high, low_open=False, high_open=False):
    """
    Call _simplest_rational_exact().  Only used to allow doctests
    on that function.
    """
    return _simplest_rational_exact(low, high, low_open, high_open)

cdef _simplest_rational_exact(Rational low, Rational high, int low_open, int high_open):
    """
    Return the simplest rational between low and high.  May return low
    or high unless low_open or high_open (respectively) are true (non-zero).
    We assume that low and high are both nonnegative, and that high>low.

    This is a helper function for simplest_rational() on RealIntervalField,
    and should not be called directly.

    EXAMPLES:
        sage: test = sage.rings.real_mpfi._simplest_rational_test_helper
        sage: test(1/4, 1/3, 0, 0)
        1/3
        sage: test(1/4, 1/3, 0, 1)
        1/4
        sage: test(1/4, 1/3, 1, 1)
        2/7
        sage: test(QQ(0), QQ(2), 0, 0)
        0
        sage: test(QQ(0), QQ(2), 1, 0)
        1
        sage: test(QQ(0), QQ(1), 1, 0)
        1
        sage: test(QQ(0), QQ(1), 1, 1)
        1/2
        sage: test(1233/1234, QQ(1), 0, 0)
        1
        sage: test(1233/1234, QQ(1), 0, 1)
        1233/1234
        sage: test(10000/32007, 10001/32007, 0, 0)
        289/925
        sage: test(QQ(0), 1/3, 1, 0)
        1/3
        sage: test(QQ(0), 1/3, 1, 1)
        1/4
        sage: test(QQ(0), 2/5, 1, 0)
        1/3
        sage: test(QQ(0), 2/5, 1, 1)
        1/3
    """
    cdef Rational r

    if low < 1:
        if low == 0:
            if low_open:
                if high > 1:
                    return Rational(1)
                inv_high = ~high
                if high_open:
                    return ~Rational(inv_high.floor() + 1)
                else:
                    return ~Rational(inv_high.ceil())
            else:
                return Rational(0)

        if high > 1:
            return Rational(1)

        r = _simplest_rational_exact(~high, ~low, high_open, low_open)
        return ~r

    fl = low.floor()
    return fl + _simplest_rational_exact(low - fl, high - fl, low_open, high_open)


def RealInterval(s, upper=None, int base=10, int pad=0, min_prec=53):
    r"""
    Return the real number defined by the string s as an element of
    \code{RealIntervalField(prec=n)}, where n potentially has slightly more
    (controlled by pad) bits than given by s.

    INPUT:
        s -- a string that defines a real number (or something whose
             string representation defines a number)
        upper -- (default: None) -- upper endpoint of interval if given, in
             which case s is the lower endpoint.
        base -- an integer between 2 and 36
        pad -- an integer >= 0.
        min_prec -- number will have at least this many bits of
                    precision, no matter what.

    EXAMPLES:
        sage: RealInterval('2.3')
        [2.2999999999999998 .. 2.3000000000000003]
        sage: RealInterval(10)
        [10.000000000000000 .. 10.000000000000000]
        sage: RealInterval('1.0000000000000000000000000000000000')
        [1.000000000000000000000000000000000000 .. 1.000000000000000000000000000000000000]
        sage: RealInterval(29308290382930840239842390482, 3^20)
        [3.48678440100000000000000000000e9 .. 2.93082903829308402398423904820e28]
    """
    if not isinstance(s, str):
        s = str(s)
    if base == 10:
        bits = int(3.32192*len(s))
    else:
        bits = int(math.log(base,2)*len(s))
    R = RealIntervalField(prec=max(bits+pad, min_prec))
    return R(s, upper, base)

# The default real interval field, with precision 53 bits
RIF = RealIntervalField()

def is_RealIntervalField(x):
    return PY_TYPE_CHECK(x, RealIntervalField)

def is_RealIntervalFieldElement(x):
    return PY_TYPE_CHECK(x, RealIntervalFieldElement)


#### pickle functions
def __create__RealIntervalField_version0(prec, sci_not):
    return RealIntervalField(prec, sci_not)

## Keep all old versions!!!
def __create__RealIntervalFieldElement_version0(parent, x, base=10):
    return RealIntervalFieldElement(parent, x, base=base)
