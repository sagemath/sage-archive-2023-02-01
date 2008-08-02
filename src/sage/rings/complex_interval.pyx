"""
Arbitrary Precision Complex Intervals

This is a simple complex interval package, using intervals which are
axis-aligned rectangles in the complex plane.  It has very few special
functions, and it does not use any special tricks to keep the size of
the intervals down.

AUTHOR:
  These authors wrote complex_number.pyx.
    -- William Stein (2006-01-26): complete rewrite
    -- Joel B. Mohler (2006-12-16): naive rewrite into pyrex
    -- William Stein(2007-01): rewrite of Mohler's rewrite
  Then complex_number.pyx was copied to complex_interval.pyx and
  heavily modified:
    -- Carl Witty (2007-10-24): rewrite to become a complex interval package
"""

#################################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import math
import operator

from sage.structure.element cimport FieldElement, RingElement, Element, ModuleElement
from complex_number cimport ComplexNumber

import complex_interval_field
from complex_field import ComplexField
import sage.misc.misc
import integer
import infinity
import real_mpfi
import real_mpfr
cimport real_mpfr

include "../ext/stdsage.pxi"

def is_ComplexIntervalFieldElement(x):
    return isinstance(x, ComplexIntervalFieldElement)

cdef class ComplexIntervalFieldElement(sage.structure.element.FieldElement):
    """
    A complex interval.

    EXAMPLES:
        sage: I = CIF.gen()
        sage: b = 1.5 + 2.5*I
        sage: loads(b.dumps()) == b
        True
    """
    cdef ComplexIntervalFieldElement _new(self):
        """
        Quickly creates a new initialized complex interval with the
        same parent as self.
        """
        cdef ComplexIntervalFieldElement x
        x = PY_NEW(ComplexIntervalFieldElement)
        x._parent = self._parent
        x._prec = self._prec
        mpfi_init2(x.__re, self._prec)
        mpfi_init2(x.__im, self._prec)
        return x

    def __init__(self, parent, real, imag=None):
        """
        Initialize a complex interval.
        """
        cdef real_mpfi.RealIntervalFieldElement rr, ii
        self._parent = parent
        self._prec = self._parent._prec

        mpfi_init2(self.__re, self._prec)
        mpfi_init2(self.__im, self._prec)

        if imag is None:
            if PY_TYPE_CHECK(real, ComplexNumber):
                real, imag = (<ComplexNumber>real).real(), (<ComplexNumber>real).imag()
            elif PY_TYPE_CHECK(real, ComplexIntervalFieldElement):
                real, imag = (<ComplexIntervalFieldElement>real).real(), (<ComplexIntervalFieldElement>real).imag()
            elif isinstance(real, sage.libs.pari.all.pari_gen):
                real, imag = real.real(), real.imag()
            elif isinstance(real, list) or isinstance(real, tuple):
                re, imag = real
                real = re
            elif isinstance(real, complex):
                real, imag = real.real, real.imag
            else:
                imag = 0
        try:
            R = parent._real_field()
            rr = R(real)
            ii = R(imag)
            mpfi_set(self.__re, <mpfi_t> rr.value)
            mpfi_set(self.__im, <mpfi_t> ii.value)
        except TypeError:
            raise TypeError, "unable to coerce to a ComplexIntervalFieldElement"


    def  __dealloc__(self):
        mpfi_clear(self.__re)
        mpfi_clear(self.__im)

    def _repr_(self):
        return self.str(10)

    def __hash__(self):
        return hash(self.str())

    def __getitem__(self, i):
        if i == 0:
            return self.real()
        elif i == 1:
            return self.imag()
        raise IndexError, "i must be between 0 and 1."

    def __reduce__( self ):
        """
        Pickling support

        EXAMPLES:
            sage: a = CIF(1 + I)
            sage: loads(dumps(a)) == a
            True
        """
        # TODO: This is potentially slow -- make a 1 version that
        # is native and much faster -- doesn't use .real()/.imag()
        return (make_ComplexIntervalFieldElement0, (self._parent, self.real(), self.imag()))

    def str(self, base=10, style=None):
        s = ""
        if not self.real().is_zero():
            s = self.real().str(base=base, style=style)
        if not self.imag().is_zero():
            y  =  self.imag()
            if s!="":
                if y < 0:
                    s = s+" - "
                    y = -y
                else:
                    s = s+" + "
            s = s+"%s*I"%y.str(base=base, style=style)
        if len(s) == 0:
            s = "0"
        return s

    def _latex_(self):
        import re
        s = self.str().replace('*I', 'i')
        return re.sub(r"e(-?\d+)", r" \\times 10^{\1}", s)

    def is_exact(self):
        return mpfr_equal_p(&self.__re.left, &self.__re.right) and \
               mpfr_equal_p(&self.__im.left, &self.__im.right)

    def diameter(self):
        """
        Returns a somewhat-arbitrarily defined "diameter" for this
        interval: returns the maximum of the diameter of the real
        and imaginary components, where diameter on a real interval
        is defined as absolute diameter if the interval contains zero,
        and relative diameter otherwise.

        EXAMPLES:
            sage: CIF(RIF(-1, 1), RIF(13, 17)).diameter()
            2.00000000000000
            sage: CIF(RIF(-0.1, 0.1), RIF(13, 17)).diameter()
            0.266666666666667
            sage: CIF(RIF(-1, 1), 15).diameter()
            2.00000000000000
        """
        cdef real_mpfr.RealNumber diam
        diam = real_mpfr.RealNumber(self._parent._real_field()._middle_field(), None)
        cdef mpfr_t tmp
        mpfr_init2(tmp, self.prec())
        mpfi_diam(diam.value, self.__re)
        mpfi_diam(tmp, self.__im)
        mpfr_max(diam.value, diam.value, tmp, GMP_RNDU)
        mpfr_clear(tmp)
        return diam

    def overlaps(self, ComplexIntervalFieldElement other):
        """
        Returns True if self and other are intervals with at least
        one value in common.

        EXAMPLES:
            sage: CIF(0).overlaps(CIF(RIF(0, 1), RIF(-1, 0)))
            True
            sage: CIF(1).overlaps(CIF(1, 1))
            False
        """
        return mpfr_greaterequal_p(&self.__re.right, &other.__re.left) \
           and mpfr_greaterequal_p(&other.__re.right, &self.__re.left) \
           and mpfr_greaterequal_p(&self.__im.right, &other.__im.left) \
           and mpfr_greaterequal_p(&other.__im.right, &self.__im.left)

    def intersection(self, other):
        """
        Returns the intersection of two complex intervals.

        EXAMPLES:
            sage: CIF(RIF(1, 3), RIF(1, 3)).intersection(CIF(RIF(2, 4), RIF(2, 4))).str(style='brackets')
            '[2.0000000000000000 .. 3.0000000000000000] + [2.0000000000000000 .. 3.0000000000000000]*I'
            sage: CIF(RIF(1, 2), RIF(1, 3)).intersection(CIF(RIF(3, 4), RIF(2, 4)))
            Traceback (most recent call last):
            ...
            ValueError: intersection of non-overlapping intervals
        """

        cdef ComplexIntervalFieldElement x = self._new()
        cdef ComplexIntervalFieldElement other_intv
        if PY_TYPE_CHECK(other, ComplexIntervalFieldElement):
            other_intv = other
        else:
            # Let type errors from _coerce_ propagate...
            other_intv = self._parent(other)

        mpfi_intersect(x.__re, self.__re, other_intv.__re)
        mpfi_intersect(x.__im, self.__im, other_intv.__im)
        if mpfr_less_p(&x.__re.right, &x.__re.left) \
           or mpfr_less_p(&x.__im.right, &x.__im.left):
            raise ValueError, "intersection of non-overlapping intervals"

        return x

    def union(self, other):
        """
        Returns the smallest complex interval including the
        two complex intervals.

        EXAMPLES:
            sage: CIF(0).union(CIF(5, 5)).str(style='brackets')
            '[0.00000000000000000 .. 5.0000000000000000] + [0.00000000000000000 .. 5.0000000000000000]*I'
        """
        cdef ComplexIntervalFieldElement x = self._new()
        cdef ComplexIntervalFieldElement other_intv
        if PY_TYPE_CHECK(other, ComplexIntervalFieldElement):
            other_intv = other
        else:
            # Let type errors from _coerce_ propagate...
            other_intv = self._parent(other)

        mpfi_union(x.__re, self.__re, other_intv.__re)
        mpfi_union(x.__im, self.__im, other_intv.__im)
        return x

    def center(self):
        """
        Returns the closest floating-point approximation to the center
        of the interval.

        EXAMPLES:
            sage: CIF(RIF(1, 2), RIF(3, 4)).center()
            1.50000000000000 + 3.50000000000000*I
        """
        cdef complex_number.ComplexNumber center
        center = complex_number.ComplexNumber(self._parent._middle_field(), None)
        mpfi_mid(center.__re, self.__re)
        mpfi_mid(center.__im, self.__im)

        return center

    def __contains__(self, other):
        """
        Test whether one interval is totally contained in another.

        EXAMPLES:
            sage: CIF(1, 1) in CIF(RIF(1, 2), RIF(1, 2))
            True
        """
        # This could be more efficient (and support more types for "other").
        return (other.real() in self.real()) and (other.imag() in self.imag())

    def contains_zero(self):
        """
        Returns True if self is an interval containing zero.

        EXAMPLES:
            sage: CIF(0).contains_zero()
            True
            sage: CIF(RIF(-1, 1), 1).contains_zero()
            False
        """
        return mpfi_has_zero(self.__re) and mpfi_has_zero(self.__im)

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        cdef ComplexIntervalFieldElement x
        x = self._new()
        mpfi_add(x.__re, self.__re, (<ComplexIntervalFieldElement>right).__re)
        mpfi_add(x.__im, self.__im, (<ComplexIntervalFieldElement>right).__im)
        return x

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        cdef ComplexIntervalFieldElement x
        x = self._new()
        mpfi_sub(x.__re, self.__re, (<ComplexIntervalFieldElement>right).__re)
        mpfi_sub(x.__im, self.__im, (<ComplexIntervalFieldElement>right).__im)
        return x

    cdef RingElement _mul_c_impl(self, RingElement right):
        cdef ComplexIntervalFieldElement x
        x = self._new()
        cdef mpfi_t t0, t1
        mpfi_init2(t0, self._prec)
        mpfi_init2(t1, self._prec)
        mpfi_mul(t0, self.__re, (<ComplexIntervalFieldElement>right).__re)
        mpfi_mul(t1, self.__im, (<ComplexIntervalFieldElement>right).__im)
        mpfi_sub(x.__re, t0, t1)
        mpfi_mul(t0, self.__re, (<ComplexIntervalFieldElement>right).__im)
        mpfi_mul(t1, self.__im, (<ComplexIntervalFieldElement>right).__re)
        mpfi_add(x.__im, t0, t1)
        mpfi_clear(t0)
        mpfi_clear(t1)
        return x

    def norm(self):
        return self.norm_c()

    cdef real_mpfi.RealIntervalFieldElement norm_c(ComplexIntervalFieldElement self):
        cdef real_mpfi.RealIntervalFieldElement x
        x = real_mpfi.RealIntervalFieldElement(self._parent._real_field(), None)

        cdef mpfi_t t0, t1
        mpfi_init2(t0, self._prec)
        mpfi_init2(t1, self._prec)

        mpfi_sqr(t0, self.__re)
        mpfi_sqr(t1, self.__im)

        mpfi_add(<mpfi_t> x.value, t0, t1)

        mpfi_clear(t0)
        mpfi_clear(t1)
        return x

    cdef real_mpfi.RealIntervalFieldElement abs_c(ComplexIntervalFieldElement self):
        cdef real_mpfi.RealIntervalFieldElement x
        x = real_mpfi.RealIntervalFieldElement(self._parent._real_field(), None)

        cdef mpfi_t t0, t1
        mpfi_init2(t0, self._prec)
        mpfi_init2(t1, self._prec)

        mpfi_sqr(t0, self.__re)
        mpfi_sqr(t1, self.__im)

        mpfi_add(<mpfi_t> x.value, t0, t1)
        mpfi_sqrt(<mpfi_t> x.value, <mpfi_t> x.value)

        mpfi_clear(t0)
        mpfi_clear(t1)
        return x

    cdef RingElement _div_c_impl(self, RingElement right):
        cdef ComplexIntervalFieldElement x
        x = self._new()
        cdef mpfi_t a, b, t0, t1, right_nm
        mpfi_init2(t0, self._prec)
        mpfi_init2(t1, self._prec)
        mpfi_init2(a, self._prec)
        mpfi_init2(b, self._prec)
        mpfi_init2(right_nm, self._prec)

        mpfi_sqr(t0, (<ComplexIntervalFieldElement>right).__re)
        mpfi_sqr(t1, (<ComplexIntervalFieldElement>right).__im)
        mpfi_add(right_nm, t0, t1)

        mpfi_div(a, (<ComplexIntervalFieldElement>right).__re, right_nm)
        mpfi_div(b, (<ComplexIntervalFieldElement>right).__im, right_nm)

        ## Do this: x.__re =  a * self.__re + b * self.__im
        mpfi_mul(t0, a, self.__re)
        mpfi_mul(t1, b, self.__im)
        mpfi_add(x.__re, t0, t1)

        ## Do this: x.__im =  a * self.__im - b * self.__re
        mpfi_mul(t0, a, self.__im)
        mpfi_mul(t1, b, self.__re)
        mpfi_sub(x.__im, t0, t1)
        mpfi_clear(t0)
        mpfi_clear(t1)
        mpfi_clear(a)
        mpfi_clear(b)
        mpfi_clear(right_nm)
        return x

    def __rdiv__(self, left):
        return ComplexIntervalFieldElement(self._parent, left)/self

    def __pow__(self, right, modulus):
        """
        Compute $x^y$.  If y is an integer, uses multiplication;
        otherwise, uses the standard definition exp(log(x)*y).
        WARNING: If the interval x crosses the negative real axis,
        then we use a non-standard definition of log() (see
        the docstring for argument() for more details).  This means
        that we will not select the principal value of the power,
        for part of the input interval (and that we violate the
        interval guarantees).

        EXAMPLES:
            sage: C.<i> = ComplexIntervalField(20)
            sage: a = i^2; a
            -1.0000000?
            sage: a.parent()
            Complex Interval Field with 20 bits of precision
            sage: a = (1+i)^7; a
            8.0000000? - 8.0000000?*I
            sage: (1+i)^(1+i)
            0.27396? + 0.58370?*I
            sage: a.parent()
            Complex Interval Field with 20 bits of precision
            sage: (2+i)^(-39)
            1.688?e-14 + 1.628?e-14*I

        If the interval crosses the negative real axis, then we don't
        use the standard branch cut (and we violate the interval
        guarantees):
            sage: (CIF(-7, RIF(-1, 1)) ^ CIF(0.3)).str(style='brackets')
            '[0.99109735947126309 .. 1.1179269966896264] + [1.4042388462787560 .. 1.4984624123369835]*I'
            sage: CIF(-7, -1) ^ CIF(0.3)
            1.117926996689626? - 1.408500714575360?*I
        """
        if isinstance(right, (int, long, integer.Integer)):
            return sage.rings.ring_element.RingElement.__pow__(self, right)
        return (self.log() * right).exp()

    def prec(self):
        """
        Return precision of this complex number.

        EXAMPLES:
            sage: i = ComplexIntervalField(2000).0
            sage: i.prec()
            2000
        """
        return self._parent.prec()

    def real(self):
        """
        Return real part of self.

        EXAMPLES:
            sage: i = ComplexIntervalField(100).0
            sage: z = 2 + 3*i
            sage: x = z.real(); x
            2.0000000000000000000000000000000?
            sage: x.parent()
            Real Interval Field with 100 bits of precision
        """
        cdef real_mpfi.RealIntervalFieldElement x
        x = real_mpfi.RealIntervalFieldElement(self._parent._real_field(), None)
        mpfi_set(x.value, self.__re)
        return x

    def imag(self):
        """
        Return imaginary part of self.

        EXAMPLES:
            sage: i = ComplexIntervalField(100).0
            sage: z = 2 + 3*i
            sage: x = z.imag(); x
            3.0000000000000000000000000000000?
            sage: x.parent()
            Real Interval Field with 100 bits of precision
        """
        cdef real_mpfi.RealIntervalFieldElement x
        x = real_mpfi.RealIntervalFieldElement(self._parent._real_field(), None)
        mpfi_set(x.value, self.__im)
        return x

    def __neg__(self):
        cdef ComplexIntervalFieldElement x
        x = self._new()
        mpfi_neg(x.__re, self.__re)
        mpfi_neg(x.__im, self.__im)
        return x

    def __pos__(self):
        return self

    def __abs__(self):
        return self.abs_c()

    def __invert__(self):
        """
        Return the multiplicative inverse.

        EXAMPLES:
            sage: I = CIF.0
            sage: a = ~(5+I)
            sage: a * (5+I)
            1.000000000000000? + 0.?e-16*I
        """
        cdef ComplexIntervalFieldElement x
        x = self._new()

        cdef mpfi_t t0, t1
        mpfi_init2(t0, self._prec)
        mpfi_init2(t1, self._prec)

        mpfi_sqr(t0, self.__re)
        mpfi_sqr(t1, self.__im)

        mpfi_add(t0, t0, t1)         # now t0 is the norm
        mpfi_div(x.__re, self.__re, t0)   #     x.__re = self.__re/norm

        mpfi_neg(t1, self.__im)
        mpfi_div(x.__im, t1, t0)  #     x.__im = -self.__im/norm

        mpfi_clear(t0)
        mpfi_clear(t1)

        return x

    def __int__(self):
        raise TypeError, "can't convert complex interval to int"

    def __long__(self):
        raise TypeError, "can't convert complex interval to long"

    def __float__(self):
        raise TypeError, "can't convert complex interval to float"

    def __complex__(self):
        raise TypeError, "can't convert complex interval to complex"

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, sage.structure.element.Element right) except -2:
        """
        Intervals are compared lexicographically on the 4-tuple
        (x.real().lower(), x.real().upper(), x.imag().lower(), x.imag().upper())

        TESTS:
            sage: tests = []
            sage: for rl in (0, 1):
            ...       for ru in (rl, rl + 1):
            ...           for il in (0, 1):
            ...               for iu in (il, il + 1):
            ...                   tests.append((CIF(RIF(rl, ru), RIF(il, iu)), (rl, ru, il, iu)))
            sage: for (i1, t1) in tests:
            ...       for (i2, t2) in tests:
            ...           assert(cmp(i1, i2) == cmp(t1, t2))
        """
        cdef int a, b
        a = mpfi_nan_p(left.__re)
        b = mpfi_nan_p((<ComplexIntervalFieldElement>right).__re)
        if a != b:
            return -1

        cdef int i
        i = mpfr_cmp(&left.__re.left, &(<ComplexIntervalFieldElement>right).__re.left)
        if i < 0:
            return -1
        elif i > 0:
            return 1
        i = mpfr_cmp(&left.__re.right, &(<ComplexIntervalFieldElement>right).__re.right)
        if i < 0:
            return -1
        elif i > 0:
            return 1
        i = mpfr_cmp(&left.__im.left, &(<ComplexIntervalFieldElement>right).__im.left)
        if i < 0:
            return -1
        elif i > 0:
            return 1
        i = mpfr_cmp(&left.__im.right, &(<ComplexIntervalFieldElement>right).__im.right)
        if i < 0:
            return -1
        elif i > 0:
            return 1
        return 0

    ########################################################################
    # Transcendental (and other) functions
    ########################################################################

    def argument(self):
        r"""
        The argument (angle) of the complex number, normalized
        so that $-\pi < \theta.lower() \leq \pi$.

        We raise a ValueError if the interval strictly contains 0,
        or if the interval contains only 0.

        WARNING: We do not always use the standard branch cut for
        argument!  If the interval crosses the negative real axis,
        then the argument will be an interval whose lower bound is
        less than $\pi$ and whose upper bound is more than $\pi$; in
        effect, we move the branch cut away from the interval.

        EXAMPLES:
            sage: i = CIF.0
            sage: (i^2).argument()
            3.141592653589794?
            sage: (1+i).argument()
            0.785398163397449?
            sage: i.argument()
            1.5707963267948967?
            sage: (-i).argument()
            -1.570796326794897?
            sage: (RR('-0.001') - i).argument()
            -1.571796326461564?
            sage: CIF(2).argument()
            0.?e-17
            sage: CIF(-2).argument()
            3.141592653589794?

        Here we see that if the interval crosses the negative real
        axis, then the argument() can exceed $\pi$, and we
        we violate the standard interval guarantees in the process:
            sage: CIF(-2, RIF(-0.1, 0.1)).argument().str(style='brackets')
            '[3.0916342578678501 .. 3.1915510493117365]'
            sage: CIF(-2, -0.1).argument()
            -3.091634257867851?
        """
        if mpfi_has_zero(self.__re) and mpfi_has_zero(self.__im):

            if mpfi_is_zero(self.__re) and mpfi_is_zero(self.__im):
                raise ValueError, "Can't take the argument of complex zero"
            if not mpfi_is_nonpos(self.__re) and not mpfi_is_nonneg(self.__re) \
               and not mpfi_is_nonpos(self.__im) and not mpfi_is_nonneg(self.__im):
                raise ValueError, "Can't take the argument of interval strictly containing zero"

            # Now if we exclude zero from the interval, we know that the
            # argument of the remaining points is bounded.  Check which
            # axes the interval extends along (we can deduce information
            # about the quadrants from information about the axes).

            which_axes = [False, False, False, False]
            if not mpfi_is_nonpos(self.__re):
                which_axes[0] = True
            if not mpfi_is_nonpos(self.__im):
                which_axes[1] = True
            if not mpfi_is_nonneg(self.__re):
                which_axes[2] = True
            if not mpfi_is_nonneg(self.__im):
                which_axes[3] = True

            lower = None
            for i in range(-1, 3):
                if which_axes[i % 4] and not which_axes[(i - 1) % 4]:
                    if lower is not None:
                        raise ValueError, "Can't take the argument of line-segment interval strictly containing zero"
                    lower = i

            for i in range(lower, lower+4):
                if which_axes[i % 4] and not which_axes[(i + 1) % 4]:
                    upper = i
                    break

            fld = self.parent()._real_field()
            return fld.pi() * fld(lower, upper) * fld(0.5)

        else:

            # OK, we know that the interval is bounded away from zero
            # in either the real or the imaginary direction (or both).
            # We'll handle the "bounded away in the imaginary direction"
            # case first.

            fld = self.parent()._real_field()

            if mpfi_is_strictly_pos(self.__im):
                return (-self.real() / self.imag()).arctan() + fld.pi()/2
            if mpfi_is_strictly_neg(self.__im):
                return (-self.real() / self.imag()).arctan() - fld.pi()/2

            if mpfi_is_strictly_pos(self.__re):
                return (self.imag() / self.real()).arctan()

            # The only remaining case is that self.__re is strictly
            # negative and self.__im contains 0.  In that case, we
            # return an interval containing pi.

            return (self.imag() / self.real()).arctan() + fld.pi()

    def arg(self):
        """
        Same as argument.

        EXAMPLES:
            sage: i = CIF.0
            sage: (i^2).arg()
            3.141592653589794?
        """
        return self.argument()

    def crosses_log_branch_cut(self):
        """
        Returns true if this interval crosses the standard branch cut
        for log() (and hence for exponentiation) and for argument.
        (Recall that this branch cut is infinitesimally below the
        negative portion of the real axis.)
        """

        if mpfi_is_nonneg(self.__re):
            return False
        if mpfi_is_nonneg(self.__im):
            return False
        if mpfi_is_neg(self.__im):
            return False
        return True

    def conjugate(self):
        """
        Return the complex conjugate of this complex number.

        EXAMPLES:
            sage: i = CIF.0
            sage: (1+i).conjugate()
            1.0000000000000000? - 1.0000000000000000?*I
        """
        cdef ComplexIntervalFieldElement x
        x = self._new()

        mpfi_set(x.__re, self.__re)
        mpfi_neg(x.__im, self.__im)
        return x

    def exp(self):
        """
        Compute exp(z).

        EXAMPLES:
            sage: i = ComplexIntervalField(300).0
            sage: z = 1 + i
            sage: z.exp()
            1.46869393991588515713896759732660426132695673662900872279767567631093696585951213872272450? + 2.28735528717884239120817190670050180895558625666835568093865811410364716018934540926734485?*I
        """
        mag = self.real().exp()
        theta = self.imag()
        re = theta.cos() * mag
        im = theta.sin() * mag
        return ComplexIntervalFieldElement(self._parent, re, im)

    def log(self):
        """
        Complex logarithm of z.  WARNING: This does always not use the
        standard branch cut for complex log!  See the docstring for
        argument() to see what we do instead.

        EXAMPLES:
            sage: a = CIF(RIF(3, 4), RIF(13, 14))
            sage: a.log().str(style='brackets')
            '[2.5908917751460420 .. 2.6782931373360067] + [1.2722973952087170 .. 1.3597029935721503]*I'
            sage: a.log().exp().str(style='brackets')
            '[2.7954667135098274 .. 4.2819545928390213] + [12.751682453911920 .. 14.237018048974635]*I'
            sage: a in a.log().exp()
            True

        If the interval crosses the negative real axis, then we don't
        use the standard branch cut (and we violate the interval guarantees):
            sage: CIF(-3, RIF(-1/4, 1/4)).log().str(style='brackets')
            '[1.0986122886681095 .. 1.1020725100903968] + [3.0584514217013518 .. 3.2247338854782349]*I'
            sage: CIF(-3, -1/4).log()
            1.102072510090397? - 3.058451421701352?*I
        """
        theta = self.argument()
        rho = abs(self)
        return ComplexIntervalFieldElement(self._parent, rho.log(), theta)

    def sqrt(self, bint all=False, **kwds):
        """
        The square root function.

        WARNING: We approximate the standard branch cut along the
        negative real axis, with sqrt(-r^2) = i*r for positive real r;
        but if the interval crosses the negative real axis, we pick
        the root with positive imaginary component for the entire
        interval.

        INPUT:
            all -- bool (default: False); if True, return a list
                of all square roots.

        EXAMPLES:
            sage: a = CIF(-1).sqrt()^2; a
            -1.000000000000000? + 0.?e-15*I
            sage: sqrt(CIF(2))
            1.414213562373095?
            sage: sqrt(CIF(-1))
            0.?e-15 + 0.9999999999999999?*I
            sage: sqrt(CIF(2-I))^2
            2.00000000000000? - 1.00000000000000?*I
            sage: CC(-2-I).sqrt()^2
            -2.00000000000000 - 1.00000000000000*I

        Here, we select a non-principal root for part of the interval, and
        violate the standard interval guarantees:
            sage: CIF(-5, RIF(-1, 1)).sqrt().str(style='brackets')
            '[-0.22250788030178321 .. 0.22250788030178296] + [2.2251857651053086 .. 2.2581008643532262]*I'
            sage: CIF(-5, -1).sqrt()
            0.222507880301783? - 2.247111425095870?*I
        """
        if self.is_zero():
            return [self] if all else self
        theta = self.argument()/2
        rho = abs(self).sqrt()
        x = ComplexIntervalFieldElement(self._parent, rho*theta.cos(), rho*theta.sin())
        if all:
            return [x, -x]
        else:
            return x


    def is_square(self):
        """
        This function always returns true as $\C$ is algebraically closed.
        """
        return True

#     def algdep(self, n, **kwds):
#         """
#         Returns a polynomial of degree at most $n$ which is approximately
#         satisfied by this complex number.  Note that the returned polynomial
#         need not be irreducible, and indeed usually won't be if $z$ is a good
#         approximation to an algebraic number of degree less than $n$.

#         ALGORITHM: Uses the PARI C-library algdep command.

#         INPUT: Type algdep? at the top level prompt. All additional
#         parameters are passed onto the top-level algdep command.

#         EXAMPLE:
#             sage: C = ComplexIntervalField()
#             sage: z = (1/2)*(1 + sqrt(3.0) *C.0); z
#             0.500000000000000 + 0.866025403784439*I
#             sage: p = z.algdep(5); p
#             x^5 + x^2
#             sage: p.factor()
#             (x + 1) * x^2 * (x^2 - x + 1)
#             sage: z^2 - z + 1
#             0.000000000000000111022302462516
#         """
#         import sage.rings.arith
#         return sage.rings.arith.algdep(self,n, **kwds)

#     def algebraic_dependancy( self, n ):
#         return self.algdep( n )

def make_ComplexIntervalFieldElement0( fld, re, im ):
    x = ComplexIntervalFieldElement( fld, re, im )
    return x



def create_ComplexIntervalFieldElement(s_real, s_imag=None, int pad=0, min_prec=53):
    r"""
    Return the complex number defined by the strings s_real and s_imag as an element of
    \code{ComplexIntervalField(prec=n)}, where n potentially has slightly more
    (controlled by pad) bits than given by s.

    INPUT:
        s_real -- a string that defines a real number (or something whose
                  string representation defines a number)
        s_imag -- a string that defines a real number (or something whose
                  string representation defines a number)
        pad -- an integer >= 0.
        min_prec -- number will have at least this many bits of precision, no matter what.

    EXAMPLES:
        sage: ComplexIntervalFieldElement('2.3')
        2.300000000000000?
        sage: ComplexIntervalFieldElement('2.3','1.1')
        2.300000000000000? + 1.1000000000000000?*I
        sage: ComplexIntervalFieldElement(10)
        10.000000000000000?
        sage: ComplexIntervalFieldElement(10,10)
        10.000000000000000? + 10.000000000000000?*I
        sage: ComplexIntervalFieldElement(1.000000000000000000000000000,2)
        1.00000000000000000000000000000? + 2.00000000000000000000000000000?*I
        sage: ComplexIntervalFieldElement(1,2.000000000000000000000)
        1.00000000000000000000000? + 2.00000000000000000000000?*I
    """
    if s_imag is None:
        s_imag = 0

    if not isinstance(s_real, str):
        s_real = str(s_real).strip()
    if not isinstance(s_imag, str):
        s_imag = str(s_imag).strip()
    #if base == 10:
    bits = max(int(3.32192*len(s_real)),int(3.32192*len(s_imag)))
    #else:
    #    bits = max(int(math.log(base,2)*len(s_imag)),int(math.log(base,2)*len(s_imag)))

    C = complex_interval_field.ComplexIntervalField(prec=max(bits+pad, min_prec))

    return ComplexIntervalFieldElement(C, s_real, s_imag)
