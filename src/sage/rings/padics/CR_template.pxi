"""
Capped relative template for complete discrete valuation rings and their fraction fields.

In order to use this template you need to write a linkage file and gluing file.
For an example see mpz_linkage.pxi (linkage file) and padic_capped_relative_element.pyx (gluing file).

The linkage file implements a common API that is then used in the class CRElement defined here.
See the documentation of mpz_linkage.pxi for the functions needed.

The gluing file does the following:

- ctypedef's celement to be the appropriate type (e.g. mpz_t)
- includes the linkage file
- includes this template
- defines a concrete class inheriting from ``CRElement``, and implements
  any desired extra methods

AUTHORS:

- David Roe (2012-3-1) -- initial version
"""

#*****************************************************************************
#       Copyright (C) 2007-2012 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# This file implements common functionality among template elements
include "padic_template_element.pxi"

from sage.structure.element cimport Element
from sage.rings.padics.common_conversion cimport comb_prec, _process_args_and_kwds
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.categories.sets_cat import Sets
from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
from sage.categories.homset import Hom

cdef inline bint exactzero(long ordp):
    """
    Whether a given valuation represents an exact zero.
    """
    return ordp >= maxordp

cdef inline int check_ordp_mpz(mpz_t ordp) except -1:
    """
    Checks for overflow after addition or subtraction of valuations.

    There is another variant, :meth:`check_ordp`, for long input.

    If overflow is detected, raises an ``OverflowError``.
    """
    if mpz_fits_slong_p(ordp) == 0 or mpz_cmp_si(ordp, maxordp) > 0 or mpz_cmp_si(ordp, minusmaxordp) < 0:
        raise OverflowError("valuation overflow")

cdef inline int assert_nonzero(CRElement x) except -1:
    """
    Checks that ``x`` is distinguishable from zero.

    Used in division and floor division.
    """
    if exactzero(x.ordp):
        raise ZeroDivisionError("cannot divide by zero")
    if x.relprec == 0:
        raise PrecisionError("cannot divide by something indistinguishable from zero.")

cdef class CRElement(pAdicTemplateElement):
    cdef int _set(self, x, long val, long xprec, absprec, relprec) except -1:
        """
        Sets the value of this element from given defining data.

        This function is intended for use in conversion, and should
        not be called on an element created with :meth:`_new_c`.

        INPUT:

        - ``x`` -- data defining a `p`-adic element: int, long,
          Integer, Rational, other `p`-adic element...

        - ``val`` -- the valuation of the resulting element

        - ``xprec -- an inherent precision of ``x``

        - ``absprec`` -- an absolute precision cap for this element

        - ``relprec`` -- a relative precision cap for this element

        TESTS::

            sage: R = Zp(5)
            sage: R(15) #indirect doctest
            3*5 + O(5^21)
            sage: R(15, absprec=5)
            3*5 + O(5^5)
            sage: R(15, relprec=5)
            3*5 + O(5^6)
            sage: R(75, absprec = 10, relprec = 9) #indirect doctest
            3*5^2 + O(5^10)
            sage: R(25/9, relprec = 5) #indirect doctest
            4*5^2 + 2*5^3 + 5^5 + 2*5^6 + O(5^7)
            sage: R(25/9, relprec = 4, absprec = 5) #indirect doctest
            4*5^2 + 2*5^3 + O(5^5)

            sage: R = Zp(5,5)
            sage: R(25/9) #indirect doctest
            4*5^2 + 2*5^3 + 5^5 + 2*5^6 + O(5^7)
            sage: R(25/9, absprec = 5)
            4*5^2 + 2*5^3 + O(5^5)
            sage: R(25/9, relprec = 4)
            4*5^2 + 2*5^3 + 5^5 + O(5^6)

            sage: R = Zp(5); S = Zp(5, 6)
            sage: S(R(17)) # indirect doctest
            2 + 3*5 + O(5^6)
            sage: S(R(17),4) # indirect doctest
            2 + 3*5 + O(5^4)
            sage: T = Qp(5); a = T(1/5) - T(1/5)
            sage: R(a)
            O(5^19)
            sage: S(a)
            O(5^19)
            sage: S(a, 17)
            O(5^17)

            sage: R = Zp(5); S = ZpCA(5)
            sage: R(S(17, 5)) #indirect doctest
            2 + 3*5 + O(5^5)
        """
        cconstruct(self.unit, self.prime_pow)
        cdef long rprec = comb_prec(relprec, self.prime_pow.prec_cap)
        cdef long aprec = comb_prec(absprec, xprec)
        if aprec <= val: # this may also hit an exact zero, if aprec == val == maxordp
            self._set_inexact_zero(aprec)
        elif exactzero(val):
            self._set_exact_zero()
        else:
            self.relprec = min(rprec, aprec - val)
            self.ordp = val
            if isinstance(x,CRElement) and x.parent() is self.parent():
                cshift(self.unit, (<CRElement>x).unit, 0, self.relprec, self.prime_pow, True)
            else:
                cconv(self.unit, x, self.relprec, val, self.prime_pow)

    cdef int _set_exact_zero(self) except -1:
        """
        Sets ``self`` as an exact zero.

        TESTS::

            sage: R = Zp(5); R(0) #indirect doctest
            0
        """
        csetzero(self.unit, self.prime_pow)
        self.ordp = maxordp
        self.relprec = 0

    cdef int _set_inexact_zero(self, long absprec) except -1:
        """
        Sets ``self`` as an inexact zero with precision ``absprec``.

        TESTS::

            sage: R = Zp(5); R(0, 5) #indirect doctest
            O(5^5)
        """
        csetzero(self.unit, self.prime_pow)
        self.ordp = absprec
        self.relprec = 0

    cdef CRElement _new_c(self):
        """
        Creates a new element with the same basic info.

        TESTS::

            sage: R = Zp(5)
            sage: R(6,5) * R(7,8) #indirect doctest
            2 + 3*5 + 5^2 + O(5^5)
        """
        cdef type t = type(self)
        cdef CRElement ans = t.__new__(t)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        cconstruct(ans.unit, ans.prime_pow)
        return ans

    cdef int check_preccap(self) except -1:
        """
        Checks that this element doesn't have precision higher than
        allowed by the precision cap.

        TESTS::

            sage: Zp(5)(1).lift_to_precision(30)
            Traceback (most recent call last):
            ...
            PrecisionError: Precision higher than allowed by the precision cap.
        """
        if self.relprec > self.prime_pow.prec_cap:
            raise PrecisionError("Precision higher than allowed by the precision cap.")

    def __copy__(self):
        """
        Return a copy of this element.

        EXAMPLES::

            sage: a = Zp(5,6)(17); b = copy(a)
            sage: a == b
            True
            sage: a is b
            False
        """
        cdef CRElement ans = self._new_c()
        ans.relprec = self.relprec
        ans.ordp = self.ordp
        ccopy(ans.unit, self.unit, ans.prime_pow)
        return ans

    cdef int _normalize(self) except -1:
        """
        Normalizes this element, so that ``self.ordp`` is correct.

        TESTS::

            sage: R = Zp(5)
            sage: R(6) + R(4) #indirect doctest
            2*5 + O(5^20)
        """
        cdef long diff
        cdef bint is_zero
        if not exactzero(self.ordp):
            is_zero = creduce(self.unit, self.unit, self.relprec, self.prime_pow)
            if is_zero:
                self._set_inexact_zero(self.ordp + self.relprec)
            else:
                diff = cremove(self.unit, self.unit, self.relprec, self.prime_pow)
                # diff is less than self.relprec since the reduction didn't yield zero
                self.ordp += diff
                check_ordp(self.ordp)
                self.relprec -= diff

    def __dealloc__(self):
        """
        Deallocate the underlying data structure.

        TESTS::

            sage: R = Zp(5)
            sage: a = R(17)
            sage: del(a)
        """
        cdestruct(self.unit, self.prime_pow)

    def __reduce__(self):
        """
        Return a tuple of a function and data that can be used to unpickle this
        element.

        TESTS::

            sage: a = ZpCR(5)(-3)
            sage: type(a)
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>
            sage: loads(dumps(a)) == a  # indirect doctest
            True
        """
        return unpickle_cre_v2, (self.__class__, self.parent(), cpickle(self.unit, self.prime_pow), self.ordp, self.relprec)

    cpdef ModuleElement _neg_(self):
        """
        Return the additive inverse of this element.

        EXAMPLES::

            sage: R = Zp(5, 20, 'capped-rel', 'val-unit')
            sage: R(5) + (-R(5)) # indirect doctest
            O(5^21)
            sage: -R(1)
            95367431640624 + O(5^20)
            sage: -R(5)
            5 * 95367431640624 + O(5^21)
            sage: -R(0)
            0
        """
        cdef CRElement ans = self._new_c()
        ans.relprec = self.relprec
        ans.ordp = self.ordp
        if ans.relprec != 0:
            cneg(ans.unit, self.unit, ans.relprec, ans.prime_pow)
            creduce(ans.unit, ans.unit, ans.relprec, ans.prime_pow)
        return ans

    cpdef ModuleElement _add_(self, ModuleElement _right):
        """
        Return the sum of this element and ``_right``.

        EXAMPLES::

            sage: R = Zp(19, 5, 'capped-rel','series')
            sage: a = R(-1); a
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: b=R(-5/2); b
            7 + 9*19 + 9*19^2 + 9*19^3 + 9*19^4 + O(19^5)
            sage: a+b #indirect doctest
            6 + 9*19 + 9*19^2 + 9*19^3 + 9*19^4 + O(19^5)
        """
        cdef CRElement ans
        cdef CRElement right = _right
        cdef long tmpL
        if self.ordp == right.ordp:
            ans = self._new_c()
            # The relative precision of the sum is the minimum of the relative precisions in this case,
            # possibly decreasing if we got cancellation
            ans.ordp = self.ordp
            ans.relprec = min(self.relprec, right.relprec)
            if ans.relprec != 0:
                cadd(ans.unit, self.unit, right.unit, ans.relprec, ans.prime_pow)
                ans._normalize()
        else:
            if self.ordp > right.ordp:
                # Addition is commutative, swap so self.ordp < right.ordp
                ans = right; right = self; self = ans
            tmpL = right.ordp - self.ordp
            if tmpL > self.relprec:
                return self
            ans = self._new_c()
            ans.ordp = self.ordp
            ans.relprec = min(self.relprec, tmpL + right.relprec)
            if ans.relprec != 0:
                cshift(ans.unit, right.unit, tmpL, ans.relprec, ans.prime_pow, False)
                cadd(ans.unit, ans.unit, self.unit, ans.relprec, ans.prime_pow)
                creduce(ans.unit, ans.unit, ans.relprec, ans.prime_pow)
        return ans

    cpdef ModuleElement _sub_(self, ModuleElement _right):
        """
        Return the difference of this element and ``_right``.

        EXAMPLES::

            sage: R = Zp(13, 4)
            sage: R(10) - R(10) #indirect doctest
            O(13^4)
            sage: R(10) - R(11)
            12 + 12*13 + 12*13^2 + 12*13^3 + O(13^4)
        """
        cdef CRElement ans
        cdef CRElement right = _right
        cdef long tmpL
        if self.ordp == right.ordp:
            ans = self._new_c()
            # The relative precision of the difference is the minimum of the relative precisions in this case,
            # possibly decreasing if we got cancellation
            ans.ordp = self.ordp
            ans.relprec = min(self.relprec, right.relprec)
            if ans.relprec != 0:
                csub(ans.unit, self.unit, right.unit, ans.relprec, ans.prime_pow)
                ans._normalize()
        elif self.ordp < right.ordp:
            tmpL = right.ordp - self.ordp
            if tmpL > self.relprec:
                return self
            ans = self._new_c()
            ans.ordp = self.ordp
            ans.relprec = min(self.relprec, tmpL + right.relprec)
            if ans.relprec != 0:
                cshift(ans.unit, right.unit, tmpL, ans.relprec, ans.prime_pow, False)
                csub(ans.unit, self.unit, ans.unit, ans.relprec, ans.prime_pow)
                creduce(ans.unit, ans.unit, ans.relprec, ans.prime_pow)
        else:
            tmpL = self.ordp - right.ordp
            if tmpL > right.relprec:
                return right._neg_()
            ans = self._new_c()
            ans.ordp = right.ordp
            ans.relprec = min(right.relprec, tmpL + self.relprec)
            if ans.relprec != 0:
                cshift(ans.unit, self.unit, tmpL, ans.relprec, ans.prime_pow, False)
                csub(ans.unit, ans.unit, right.unit, ans.relprec, ans.prime_pow)
                creduce(ans.unit, ans.unit, ans.relprec, ans.prime_pow)
        return ans

    def __invert__(self):
        r"""
        Returns the multiplicative inverse of this element.

        .. NOTE::

            The result of inversion always lives in the fraction
            field, even if the element to be inverted is a unit.

        EXAMPLES::

            sage: R = Qp(7,4,'capped-rel','series'); a = R(3); a
            3 + O(7^4)
            sage: ~a   # indirect doctest
            5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)
        """
        assert_nonzero(self)
        cdef CRElement ans = self._new_c()
        if ans.prime_pow.in_field == 0:
            ans._parent = self._parent.fraction_field()
            ans.prime_pow = ans._parent.prime_pow
        ans.ordp = -self.ordp
        ans.relprec = self.relprec
        cinvert(ans.unit, self.unit, ans.relprec, ans.prime_pow)
        return ans

    cpdef RingElement _mul_(self, RingElement _right):
        r"""
        Return the product of this element and ``_right``.

        EXAMPLES::

            sage: R = Zp(5)
            sage: a = R(2385,11); a
            2*5 + 4*5^3 + 3*5^4 + O(5^11)
            sage: b = R(2387625, 16); b
            5^3 + 4*5^5 + 2*5^6 + 5^8 + 5^9 + O(5^16)
            sage: a * b # indirect doctest
            2*5^4 + 2*5^6 + 4*5^7 + 2*5^8 + 3*5^10 + 5^11 + 3*5^12 + 4*5^13 + O(5^14)
        """
        cdef CRElement ans
        cdef CRElement right = _right
        if exactzero(self.ordp):
            return self
        if exactzero(right.ordp):
            return right
        ans = self._new_c()
        ans.relprec = min(self.relprec, right.relprec)
        if ans.relprec == 0:
            ans._set_inexact_zero(self.ordp + right.ordp)
        else:
            ans.ordp = self.ordp + right.ordp
            cmul(ans.unit, self.unit, right.unit, ans.relprec, ans.prime_pow)
            creduce(ans.unit, ans.unit, ans.relprec, ans.prime_pow)
        check_ordp(ans.ordp)
        return ans

    cpdef RingElement _div_(self, RingElement _right):
        """
        Return the quotient of this element and ``right``.

        .. NOTE::

            The result of division always lives in the fraction field,
            even if the element to be inverted is a unit.

        EXAMPLES::

            sage: R = Zp(5,6)
            sage: R(17) / R(21) #indirect doctest
            2 + 4*5^2 + 3*5^3 + 4*5^4 + O(5^6)
            sage: a = R(50) / R(5); a
            2*5 + O(5^7)
            sage: R(5) / R(50)
            3*5^-1 + 2 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + O(5^5)
            sage: ~a
            3*5^-1 + 2 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + O(5^5)
            sage: 1 / a
            3*5^-1 + 2 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + O(5^5)
        """
        cdef CRElement ans
        cdef CRElement right = _right
        assert_nonzero(right)
        ans = self._new_c()
        if ans.prime_pow.in_field == 0:
            ans._parent = self._parent.fraction_field()
            ans.prime_pow = ans._parent.prime_pow
        if exactzero(self.ordp):
            ans._set_exact_zero()
            return ans
        ans.relprec = min(self.relprec, right.relprec)
        if ans.relprec == 0:
            ans._set_inexact_zero(self.ordp - right.ordp)
        else:
            ans.ordp = self.ordp - right.ordp
            cdivunit(ans.unit, self.unit, right.unit, ans.relprec, ans.prime_pow)
            creduce(ans.unit, ans.unit, ans.relprec, ans.prime_pow)
        check_ordp(ans.ordp)
        return ans

    def __pow__(CRElement self, _right, dummy):
        r"""
        Exponentiation.

        When ``right`` is divisible by `p` then one can get more
        precision than expected.

        Lemma 2.1 [SP]_:

        Let `\alpha` be in `\mathcal{O}_K`.  Let

        ..math ::

            p = -\pi_K^{e_K} \epsilon

        be the factorization of `p` where `\epsilon` is a unit.  Then
        the `p`-th power of `1 + \alpha \pi_K^{\lambda}` satisfies

        ..math ::

            (1 + \alpha \pi^{\lambda})^p \equiv \left{ \begin{array}{lll}
            1 + \alpha^p \pi_K^{p \lambda} &
             \mod \mathfrak{p}_K^{p \lambda + 1} &
             \mbox{if $1 \le \lambda < \frac{e_K}{p-1}$} \\
            1 + (\alpha^p - \epsilon \alpha) \pi_K^{p \lambda} &
             \mod \mathfrak{p}_K^{p \lambda + 1} &
             \mbox{if $\lambda = \frac{e_K}{p-1}$} \\
            1 - \epsilon \alpha \pi_K^{\lambda + e} &
             \mod \mathfrak{p}_K^{\lambda + e + 1} &
             \mbox{if $\lambda > \frac{e_K}{p-1}$}
            \end{array} \right.


        So if ``right`` is divisible by `p^k` we can multiply the
        relative precision by `p` until we exceed `e/(p-1)`, then add
        `e` until we have done a total of `k` things: the precision of
        the result can therefore be greater than the precision of
        ``self``.

        For `\alpha` in `\ZZ_p` we can simplify the result a bit.  In
        this case, the `p`-th power of `1 + \alpha p^{\lambda}`
        satisfies

        .. MATH::

            (1 + \alpha p^{\lambda})^p \equiv 1 + \alpha p^{\lambda + 1} mod p^{\lambda + 2}

        unless `\lambda = 1` and `p = 2`, in which case

        .. MATH::

            (1 + 2 \alpha)^2 \equiv 1 + 4(\alpha^2 + \alpha) mod 8

        So for `p \ne 2`, if right is divisible by `p^k` then we add
        `k` to the relative precision of the answer.

        For `p = 2`, if we start with something of relative precision
        1 (ie `2^m + O(2^{m+1})`), `\alpha^2 + \alpha \equiv 0 \mod
        2`, so the precision of the result is `k + 2`:

        .. MATH::

            (2^m + O(2^{m+1}))^{2^k} = 2^{m 2^k} + O(2^{m 2^k + k + 2})

        For `p`-adic exponents, we define `\alpha^\beta` as
        `\exp(\beta \log(\alpha))`.  The precision of the result is
        determined using the power series expansions for the
        exponential and logarithm maps, together with the notes above.

        .. NOTE::

            For `p`-adic exponents we always need that `a` is a unit.
            For unramified extensions `a^b` will converge as long as
            `b` is integral (though it may converge for non-integral
            `b` as well depending on the value of `a`).  However, in
            highly ramified extensions some bases may be sufficiently
            close to `1` that `exp(b log(a))` does not converge even
            though `b` is integral.

        .. WARNING::

            If `\alpha` is a unit, but not congruent to `1` modulo
            `\pi_K`, the result will not be the limit over integers
            `b` converging to `\beta` since this limit does not exist.
            Rather, the logarithm kills torsion in `\ZZ_p^\times`, and
            `\alpha^\beta` will equal `(\alpha')^\beta`, where
            `\alpha'` is the quotient of `\alpha` by the Teichmuller
            representative congruent to `\alpha` modulo `\pi_K`.  Thus
            the result will always be congruent to `1` modulo `\pi_K`.

        .. REFERENCES::

        .. [SP] Constructing Class Fields over Local Fields.
                Sebastian Pauli.

        INPUT:

        - ``_right`` -- currently integers and `p`-adic exponents are
          supported.

        - ``dummy`` -- not used (Python's ``__pow__`` signature
          includes it)

        EXAMPLES::

            sage: R = Zp(19, 5, 'capped-rel','series')
            sage: a = R(-1); a
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: a^2    # indirect doctest
            1 + O(19^5)
            sage: a^3
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: R(5)^30
            11 + 14*19 + 19^2 + 7*19^3 + O(19^5)
            sage: K = Qp(19, 5, 'capped-rel','series')
            sage: a = K(-1); a
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: a^2
            1 + O(19^5)
            sage: a^3
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: K(5)^30
            11 + 14*19 + 19^2 + 7*19^3 + O(19^5)
            sage: K(5, 3)^19 #indirect doctest
            5 + 3*19 + 11*19^3 + O(19^4)

        `p`-adic exponents are also supported::

            sage: a = K(8/5,4); a
            13 + 7*19 + 11*19^2 + 7*19^3 + O(19^4)
            sage: a^(K(19/7))
            1 + 14*19^2 + 11*19^3 + 13*19^4 + O(19^5)
            sage: (a // K.teichmuller(13))^(K(19/7))
            1 + 14*19^2 + 11*19^3 + 13*19^4 + O(19^5)
            sage: (a.log() * 19/7).exp()
            1 + 14*19^2 + 11*19^3 + 13*19^4 + O(19^5)
        """
        cdef long base_level, exp_prec
        cdef mpz_t tmp
        cdef Integer right
        cdef CRElement base, pright, ans
        cdef bint exact_exp
        if (isinstance(_right, Integer) or isinstance(_right, (int, long)) or isinstance(_right, Rational)):
            if _right < 0:
                base = ~self
                return base.__pow__(-_right, dummy)
            exact_exp = True
        elif self.parent() is _right.parent():
            ## For extension elements, we need to switch to the
            ## fraction field sometimes in highly ramified extensions.
            exact_exp = False
            pright = _right
        else:
            self, _right = canonical_coercion(self, _right)
            return self.__pow__(_right, dummy)
        if exact_exp and _right == 0:
            # return 1 to maximum precision
            ans = self._new_c()
            ans.ordp = 0
            ans.relprec = self.prime_pow.ram_prec_cap
            csetone(ans.unit, ans.prime_pow)
            return ans
        if exactzero(self.ordp):
            if exact_exp:
                # We may assume from above that right > 0
                return self
            else:
                # log(0) is not defined
                raise ValueError("0^x is not defined for p-adic x: log(0) does not converge")
        ans = self._new_c()
        if self.relprec == 0:
            # If a positive integer exponent, return an inexact zero of valuation right * self.ordp.  Otherwise raise an error.
            if isinstance(_right, (int, long)):
                _right = Integer(_right)
            if isinstance(_right, Integer):
                right = <Integer>_right
                mpz_init(tmp)
                mpz_mul_si(tmp, (<Integer>_right).value, self.ordp)
                check_ordp_mpz(tmp)
                ans._set_inexact_zero(mpz_get_si(tmp))
                mpz_clear(tmp)
            else:
                raise PrecisionError
        elif exact_exp:
            # exact_pow_helper is defined in padic_template_element.pxi
            right = exact_pow_helper(&ans.relprec, self.relprec, _right, self.prime_pow)
            if ans.relprec > self.prime_pow.ram_prec_cap:
                ans.relprec = self.prime_pow.ram_prec_cap
            mpz_init(tmp)
            mpz_mul_si(tmp, right.value, self.ordp)
            check_ordp_mpz(tmp)
            ans.ordp = mpz_get_si(tmp)
            mpz_clear(tmp)
            cpow(ans.unit, self.unit, right.value, ans.relprec, ans.prime_pow)
        else:
            # padic_pow_helper is defined in padic_template_element.pxi
            ans.relprec = padic_pow_helper(ans.unit, self.unit, self.ordp, self.relprec,
                                           pright.unit, pright.ordp, pright.relprec, self.prime_pow)
            ans.ordp = 0
        return ans

    cdef pAdicTemplateElement _lshift_c(self, long shift):
        """
        Multiplies by `\pi^{\mbox{shift}}`.

        Negative shifts may truncate the result if the parent is not a
        field.

        TESTS::

            sage: a = Zp(5)(17); a
            2 + 3*5 + O(5^20)
            sage: a << 2 #indirect doctest
            2*5^2 + 3*5^3 + O(5^22)
            sage: a << -2
            O(5^18)
            sage: a << 0 == a
            True
            sage: Zp(5)(0) << -4000
            0
        """
        if exactzero(self.ordp):
            return self
        if self.prime_pow.in_field == 0 and shift < 0 and -shift > self.ordp:
            return self._rshift_c(-shift)
        cdef CRElement ans = self._new_c()
        ans.relprec = self.relprec
        ans.ordp = self.ordp + shift
        check_ordp(ans.ordp)
        ccopy(ans.unit, self.unit, ans.prime_pow)
        return ans

    cdef pAdicTemplateElement _rshift_c(self, long shift):
        """
        Divides by ``\pi^{\mbox{shift}}``.

        Positive shifts may truncate the result if the parent is not a
        field.

        TESTS::

            sage: R = Zp(5); K = Qp(5)
            sage: R(17) >> 1
            3 + O(5^19)
            sage: K(17) >> 1
            2*5^-1 + 3 + O(5^19)
            sage: R(17) >> 40
            O(5^0)
            sage: K(17) >> -5
            2*5^5 + 3*5^6 + O(5^25)
        """
        if exactzero(self.ordp):
            return self
        cdef CRElement ans = self._new_c()
        cdef long diff
        if self.prime_pow.in_field == 1 or shift <= self.ordp:
            ans.relprec = self.relprec
            ans.ordp = self.ordp - shift
            check_ordp(ans.ordp)
            ccopy(ans.unit, self.unit, ans.prime_pow)
        else:
            diff = shift - self.ordp
            if diff >= self.relprec:
                ans._set_inexact_zero(0)
            else:
                ans.relprec = self.relprec - diff
                cshift(ans.unit, self.unit, -diff, ans.relprec, ans.prime_pow, False)
                ans.ordp = 0
                ans._normalize()
        return ans

    cpdef RingElement _floordiv_(self, RingElement _right):
        """
        Floor division.

        TESTS::

            sage: r = Zp(19)
            sage: a = r(1+19+17*19^3+5*19^4); b = r(19^3); a/b
            19^-3 + 19^-2 + 17 + 5*19 + O(19^17)
            sage: a//b     # indirect doctest
            17 + 5*19 + O(19^17)

            sage: R = Zp(19, 5, 'capped-rel','series')
            sage: a = R(-1); a
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: b=R(-2*19^3); b
            17*19^3 + 18*19^4 + 18*19^5 + 18*19^6 + 18*19^7 + O(19^8)
            sage: a//b # indirect doctest
            9 + 9*19 + O(19^2)

            sage: R = Zp(5,5)
            sage: R(28937) // R(75) # indirect doctest
            4 + 3*5 + 3*5^2 + O(5^3)

            sage: R(0,12) // R(175,3)
            O(5^10)
        """
        if exactzero(self.ordp):
            return self
        cdef CRElement right = _right
        assert_nonzero(right)
        cdef CRElement ans = self._new_c()
        cdef long diff = self.ordp - right.ordp
        if self.relprec == 0:
            ans.ordp = diff
            ans.relprec = 0
            csetzero(ans.unit, ans.prime_pow)
        elif diff >= 0 or self.prime_pow.in_field:
            ans.ordp = diff
            ans.relprec = min(self.relprec, right.relprec)
            cdivunit(ans.unit, self.unit, right.unit, ans.relprec, ans.prime_pow)
            creduce(ans.unit, ans.unit, ans.relprec, ans.prime_pow)
        else:
            ans.ordp = 0
            ans.relprec = min(self.relprec, right.relprec) + diff
            if ans.relprec < 0:
                ans.relprec = 0
                csetzero(ans.unit, ans.prime_pow)
            else:
                cdivunit(ans.unit, self.unit, right.unit, ans.relprec - diff, ans.prime_pow)
                cshift(ans.unit, ans.unit, diff, ans.relprec, ans.prime_pow, False)
                ans._normalize()
        return ans

    def add_bigoh(self, absprec):
        """
        Returns a new element with absolute precision decreased to
        ``absprec``.

        INPUT:

        - ``absprec`` -- an integer or infinity

        OUTPUT:

        an equal element with precision set to the minimum of ``self's``
        precision and ``absprec``

        EXAMPLE::

            sage: R = Zp(7,4,'capped-rel','series'); a = R(8); a.add_bigoh(1)
            1 + O(7)
            sage: b = R(0); b.add_bigoh(3)
            O(7^3)
            sage: R = Qp(7,4); a = R(8); a.add_bigoh(1)
            1 + O(7)
            sage: b = R(0); b.add_bigoh(3)
            O(7^3)

            The precision never increases::

            sage: R(4).add_bigoh(2).add_bigoh(4)
            4 + O(7^2)

            Another example that illustrates that the precision does
            not increase::

            sage: k = Qp(3,5)
            sage: a = k(1234123412/3^70); a
            2*3^-70 + 3^-69 + 3^-68 + 3^-67 + O(3^-65)
            sage: a.add_bigoh(2)
            2*3^-70 + 3^-69 + 3^-68 + 3^-67 + O(3^-65)

            sage: k = Qp(5,10)
            sage: a = k(1/5^3 + 5^2); a
            5^-3 + 5^2 + O(5^7)
            sage: a.add_bigoh(2)
            5^-3 + O(5^2)
            sage: a.add_bigoh(-1)
            5^-3 + O(5^-1)
        """
        cdef CRElement ans
        cdef long aprec, newprec
        if absprec is infinity:
            return self
        elif isinstance(absprec, int):
            aprec = absprec
        else:
            if not isinstance(absprec, Integer):
                absprec = Integer(absprec)
            aprec = mpz_get_si((<Integer>absprec).value)
        if aprec < self.ordp:
            ans = self._new_c()
            ans._set_inexact_zero(aprec)
        elif aprec >= self.ordp + self.relprec:
            ans = self
        else:
            ans = self._new_c()
            ans.ordp = self.ordp
            ans.relprec = aprec - self.ordp
            creduce(ans.unit, self.unit, ans.relprec, ans.prime_pow)
        return ans

    cpdef bint _is_exact_zero(self) except -1:
        """
        Returns true if this element is exactly zero.

        EXAMPLES::

            sage: R = Zp(5)
            sage: R(0)._is_exact_zero()
            True
            sage: R(0,5)._is_exact_zero()
            False
            sage: R(17)._is_exact_zero()
            False
        """
        return exactzero(self.ordp)

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Returns True if this element is indistinguishable from zero
        but has finite precision.

        EXAMPLES::

            sage: R = Zp(5)
            sage: R(0)._is_inexact_zero()
            False
            sage: R(0,5)._is_inexact_zero()
            True
            sage: R(17)._is_inexact_zero()
            False
        """
        return self.relprec == 0 and not exactzero(self.ordp)

    def is_zero(self, absprec = None):
        r"""
        Determines whether this element is zero modulo
        `\pi^{\mbox{absprec}}`.

        If ``absprec is None``, returns ``True`` if this element is
        indistinguishable from zero.

        INPUT:

        - ``absprec`` -- an integer, infinity, or ``None``

        EXAMPLES::

            sage: R = Zp(5); a = R(0); b = R(0,5); c = R(75)
            sage: a.is_zero(), a.is_zero(6)
            (True, True)
            sage: b.is_zero(), b.is_zero(5)
            (True, True)
            sage: c.is_zero(), c.is_zero(2), c.is_zero(3)
            (False, True, False)
            sage: b.is_zero(6)
            Traceback (most recent call last):
            ...
            PrecisionError: Not enough precision to determine if element is zero
        """
        if absprec is None:
            return self.relprec == 0
        if exactzero(self.ordp):
            return True
        if absprec is infinity:
            return False
        if isinstance(absprec, int):
            if self.relprec == 0 and absprec > self.ordp:
                raise PrecisionError("Not enough precision to determine if element is zero")
            return self.ordp >= absprec
        if not isinstance(absprec, Integer):
            absprec = Integer(absprec)
        if self.relprec == 0:
            if mpz_cmp_si((<Integer>absprec).value, self.ordp) > 0:
                raise PrecisionError("Not enough precision to determine if element is zero")
            else:
                return True
        return mpz_cmp_si((<Integer>absprec).value, self.ordp) <= 0

    def __nonzero__(self):
        """
        Returns True if self is distinguishable from zero.

        For most applications, explicitly specifying the power of p
        modulo which the element is supposed to be nonzero is
        preferable.

        EXAMPLES::

            sage: R = Zp(5); a = R(0); b = R(0,5); c = R(75)
            sage: bool(a), bool(b), bool(c)
            (False, False, True)
        """
        return self.relprec != 0

    def is_equal_to(self, _right, absprec=None):
        r"""
        Returns whether self is equal to right modulo
        `\pi^{\mbox{absprec}}`.

        If ``absprec is None``, returns True if self and right are
        equal to the minimum of their precisions.

        INPUT:

        - ``right`` -- a `p`-adic element
        - ``absprec`` -- an integer, infinity, or ``None``

        EXAMPLES::

            sage: R = Zp(5, 10); a = R(0); b = R(0, 3); c = R(75, 5)
            sage: aa = a + 625; bb = b + 625; cc = c + 625
            sage: a.is_equal_to(aa), a.is_equal_to(aa, 4), a.is_equal_to(aa, 5)
            (False, True, False)
            sage: a.is_equal_to(aa, 15)
            Traceback (most recent call last):
            ...
            PrecisionError: Elements not known to enough precision

            sage: a.is_equal_to(a, 50000)
            True

            sage: a.is_equal_to(b), a.is_equal_to(b, 2)
            (True, True)
            sage: a.is_equal_to(b, 5)
            Traceback (most recent call last):
            ...
            PrecisionError: Elements not known to enough precision

            sage: b.is_equal_to(b, 5)
            Traceback (most recent call last):
            ...
            PrecisionError: Elements not known to enough precision

            sage: b.is_equal_to(bb, 3)
            True
            sage: b.is_equal_to(bb, 4)
            Traceback (most recent call last):
            ...
            PrecisionError: Elements not known to enough precision

            sage: c.is_equal_to(b, 2), c.is_equal_to(b, 3)
            (True, False)
            sage: c.is_equal_to(b, 4)
            Traceback (most recent call last):
            ...
            PrecisionError: Elements not known to enough precision

            sage: c.is_equal_to(cc, 2), c.is_equal_to(cc, 4), c.is_equal_to(cc, 5)
            (True, True, False)

        TESTS::

            sage: aa.is_equal_to(a), aa.is_equal_to(a, 4), aa.is_equal_to(a, 5)
            (False, True, False)
            sage: aa.is_equal_to(a, 15)
            Traceback (most recent call last):
            ...
            PrecisionError: Elements not known to enough precision

            sage: b.is_equal_to(a), b.is_equal_to(a, 2)
            (True, True)
            sage: b.is_equal_to(a, 5)
            Traceback (most recent call last):
            ...
            PrecisionError: Elements not known to enough precision

            sage: bb.is_equal_to(b, 3)
            True
            sage: bb.is_equal_to(b, 4)
            Traceback (most recent call last):
            ...
            PrecisionError: Elements not known to enough precision

            sage: b.is_equal_to(c, 2), b.is_equal_to(c, 3)
            (True, False)
            sage: b.is_equal_to(c, 4)
            Traceback (most recent call last):
            ...
            PrecisionError: Elements not known to enough precision

            sage: cc.is_equal_to(c, 2), cc.is_equal_to(c, 4), cc.is_equal_to(c, 5)
            (True, True, False)
        """
        cdef CRElement right
        cdef long aprec, rprec
        if self.parent() is _right.parent():
            right = _right
        else:
            right = self.parent().coerce(_right)
        if exactzero(self.ordp) and exactzero(right.ordp):
            return True
        elif absprec is infinity:
            raise PrecisionError("Elements not known to enough precision")
        if absprec is None:
            aprec = min(self.ordp + self.relprec, right.ordp + right.relprec)
        else:
            if not isinstance(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                if mpz_sgn((<Integer>absprec).value) < 0 or \
                   exactzero(self.ordp) and exactzero(right.ordp):
                    return True
                else:
                    raise PrecisionError("Elements not known to enough precision")
            aprec = mpz_get_si((<Integer>absprec).value)
            if aprec > self.ordp + self.relprec or aprec > right.ordp + right.relprec:
                raise PrecisionError("Elements not known to enough precision")
        if self.ordp >= aprec and right.ordp >= aprec:
            return True
        elif self.ordp != right.ordp:
            return False
        rprec = aprec - self.ordp
        return ccmp(self.unit, right.unit, rprec, rprec < self.relprec, rprec < right.relprec, self.prime_pow) == 0

    cdef int _cmp_units(self, pAdicGenericElement _right) except -2:
        """
        Comparison of units, used in equality testing.

        EXAMPLES::

            sage: R = Zp(5)
            sage: a = R(17); b = R(0,3); c = R(85,7); d = R(2, 1)
            sage: any([a == b, a == c, b == c, b == d, c == d])
            False
            sage: all([a == a, b == b, c == c, d == d, a == d])
            True

            sage: sorted([a, b, c, d])
            [2 + 3*5 + O(5^20), 2 + O(5), 2*5 + 3*5^2 + O(5^7), O(5^3)]
        """
        cdef CRElement right = _right
        cdef long rprec = min(self.relprec, right.relprec)
        if rprec == 0:
            return 0
        return ccmp(self.unit, right.unit, rprec, rprec < self.relprec, rprec < right.relprec, self.prime_pow)

    cdef pAdicTemplateElement lift_to_precision_c(self, long absprec):
        """
        Lifts this element to another with precision at least ``absprec``.

        TESTS::

            sage: R = Zp(5); a = R(0); b = R(0,5); c = R(17,3)
            sage: a.lift_to_precision(5)
            0
            sage: b.lift_to_precision(4)
            O(5^5)
            sage: b.lift_to_precision(8)
            O(5^8)
            sage: b.lift_to_precision(40)
            O(5^40)
            sage: c.lift_to_precision(1)
            2 + 3*5 + O(5^3)
            sage: c.lift_to_precision(8)
            2 + 3*5 + O(5^8)
            sage: c.lift_to_precision(40)
            Traceback (most recent call last):
            ...
            PrecisionError: Precision higher than allowed by the precision cap.
        """
        cpdef CRElement ans
        if absprec == maxordp:
            if self.relprec == 0:
                ans = self._new_c()
                ans._set_exact_zero()
                return ans
            else:
                absprec = self.ordp + self.prime_pow.prec_cap
        cdef long relprec = absprec - self.ordp
        if relprec <= self.relprec:
            return self
        ans = self._new_c()
        if self.relprec == 0:
            ans._set_inexact_zero(absprec)
        else:
            ans.ordp = self.ordp
            ans.relprec = relprec
            ccopy(ans.unit, self.unit, ans.prime_pow)
        return ans

    def _cache_key(self):
        r"""
        Return a hashable key which identifies this element for caching.

        TESTS::

            sage: K.<a> = Qq(9)
            sage: (9*a)._cache_key()
            (..., ((0, 1),), 2, 20)

        .. SEEALSO::

            :meth:`sage.misc.cachefunc._cache_key`
        """
        tuple_recursive = lambda l: tuple(tuple_recursive(x) for x in l) if isinstance(l, list) else l
        return (self.parent(), tuple_recursive(self.list()), self.valuation(), self.precision_relative())

    def list(self, lift_mode = 'simple', start_val = None):
        """
        Returns a list of coefficients in a power series expansion of
        self in terms of `\pi`.  If self is a field element, they start at
        `\pi^{\mbox{valuation}}`, if a ring element at `\pi^0`.

        For each lift mode, this function returns a list of `a_i` so
        that this element can be expressed as

        .. MATH::

            \pi^v \cdot \sum_{i=0}^\infty a_i \pi^i

        where `v` is the valuation of this element when the parent is
        a field, and `v = 0` otherwise.

        Different lift modes affect the choice of `a_i`.  When
        ``lift_mode`` is ``'simple'``, the resulting `a_i` will be
        non-negative: if the residue field is `\mathbb{F}_p` then they
        will be integers with `0 \le a_i < p`; otherwise they will be
        a list of integers in the same range giving the coefficients
        of a polynomial in the indeterminant representing the maximal
        unramified subextension.

        Choosing ``lift_mode`` as ``'smallest'`` is similar to
        ``'simple'``, but uses a balanced representation `-p/2 < a_i
        \le p/2`.

        Finally, setting ``lift_mode = 'teichmuller'`` will yield
        Teichmuller representatives for the `a_i`: `a_i^q = a_i`.  In
        this case the `a_i` will also be `p`-adic elements.

        INPUT:

        - ``lift_mode`` -- ``'simple'``, ``'smallest'`` or
          ``'teichmuller'`` (default: ``'simple'``)

        - ``start_val`` -- start at this valuation rather than the
          default (`0` or the valuation of this element).  If
          ``start_val`` is larger than the valuation of this element
          a ``ValueError`` is raised.

        OUTPUT:

        - the list of coefficients of this element.  For base elements
          these will be integers if ``lift_mode`` is ``'simple'`` or
          ``'smallest'``, and elements of ``self.parent()`` if
          ``lift_mode`` is ``'teichmuller'``.

        .. NOTE::

            Use slice operators to get a particular range.

        EXAMPLES::

            sage: R = Zp(7,6); a = R(12837162817); a
            3 + 4*7 + 4*7^2 + 4*7^4 + O(7^6)
            sage: L = a.list(); L
            [3, 4, 4, 0, 4]
            sage: sum([L[i] * 7^i for i in range(len(L))]) == a
            True
            sage: L = a.list('smallest'); L
            [3, -3, -2, 1, -3, 1]
            sage: sum([L[i] * 7^i for i in range(len(L))]) == a
            True
            sage: L = a.list('teichmuller'); L
            [3 + 4*7 + 6*7^2 + 3*7^3 + 2*7^5 + O(7^6),
            0,
            5 + 2*7 + 3*7^3 + O(7^4),
            1 + O(7^3),
            3 + 4*7 + O(7^2),
            5 + O(7)]
            sage: sum([L[i] * 7^i for i in range(len(L))])
            3 + 4*7 + 4*7^2 + 4*7^4 + O(7^6)

            sage: R(0, 7).list()
            []

            sage: R = Qp(7,4); a = R(6*7+7**2); a.list()
            [6, 1]
            sage: a.list('smallest')
            [-1, 2]
            sage: a.list('teichmuller')
            [6 + 6*7 + 6*7^2 + 6*7^3 + O(7^4),
            2 + 4*7 + 6*7^2 + O(7^3),
            3 + 4*7 + O(7^2),
            3 + O(7)]

        TESTS:

        Check to see that :trac:`10292` is resolved::

            sage: E = EllipticCurve('37a')
            sage: R = E.padic_regulator(7)
            sage: len(R.list())
            19
        """
        if start_val is not None and start_val > self.ordp:
            raise ValueError("starting valuation must be smaller than the element's valuation.  See slice()")
        if self.relprec == 0:
            return []
        if lift_mode == 'teichmuller':
            ulist = self.teichmuller_list()
        elif lift_mode == 'simple':
            ulist = clist(self.unit, self.relprec, True, self.prime_pow)
        elif lift_mode == 'smallest':
            ulist = clist(self.unit, self.relprec, False, self.prime_pow)
        else:
            raise ValueError("unknown lift_mode")
        if (self.prime_pow.in_field == 0 and self.ordp > 0) or start_val is not None:
            if lift_mode == 'teichmuller':
                zero = self.parent()(0)
            else:
                # needs to be defined in the linkage file.
                zero = _list_zero
            if start_val is None:
                v = self.ordp
            else:
                v = self.ordp - start_val
            ulist = [zero] * v + ulist
        return ulist

    def teichmuller_list(self):
        r"""
        Returns a list [`a_0`, `a_1`,..., `a_n`] such that

        - `a_i^q = a_i`, where `q` is the cardinality of the residue field,

        - ``self.unit_part() =`` `\sum_{i = 0}^n a_i p^i`, and

        - if `a_i \ne 0`, the absolute precision of `a_i` is
          self.precision_relative() - i

        EXAMPLES::

            sage: R = Qp(5,5); R(70).list('teichmuller') #indirect doctest
            [4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5),
            3 + 3*5 + 2*5^2 + 3*5^3 + O(5^4),
            2 + 5 + 2*5^2 + O(5^3),
            1 + O(5^2),
            4 + O(5)]
        """
        ans = PyList_New(0)
        if self.relprec == 0:
            return ans
        cdef long curpower = self.relprec
        cdef CRElement list_elt
        cdef CRElement tmp = self._new_c()
        ccopy(tmp.unit, self.unit, self.prime_pow)
        while not ciszero(tmp.unit, tmp.prime_pow) and curpower > 0:
            list_elt = self._new_c()
            cteichmuller(list_elt.unit, tmp.unit, curpower, self.prime_pow)
            if ciszero(list_elt.unit, self.prime_pow):
                list_elt._set_exact_zero()
                cshift_notrunc(tmp.unit, tmp.unit, -1, curpower-1, self.prime_pow)
            else:
                list_elt.ordp = 0
                list_elt.relprec = curpower
                csub(tmp.unit, tmp.unit, list_elt.unit, curpower, self.prime_pow)
                cshift_notrunc(tmp.unit, tmp.unit, -1, curpower-1, self.prime_pow)
                creduce(tmp.unit, tmp.unit, curpower-1, self.prime_pow)
            curpower -= 1
            PyList_Append(ans, list_elt)
        return ans

    def _teichmuller_set_unsafe(self):
        """
        Sets this element to the Teichmuller representative with the
        same residue.

        .. WARNING::

            This function modifies the element, which is not safe.
            Elements are supposed to be immutable.

        EXAMPLES::

            sage: R = Zp(17,5); a = R(11)
            sage: a
            11 + O(17^5)
            sage: a._teichmuller_set_unsafe(); a
            11 + 14*17 + 2*17^2 + 12*17^3 + 15*17^4 + O(17^5)
            sage: a.list('teichmuller')
            [11 + 14*17 + 2*17^2 + 12*17^3 + 15*17^4 + O(17^5)]

        Note that if you set an element which is congruent to 0 you
        get an exact 0.

            sage: b = R(17*5); b
            5*17 + O(17^6)
            sage: b._teichmuller_set_unsafe(); b
            0
        """
        if self.ordp > 0:
            self._set_exact_zero()
        elif self.ordp < 0:
            raise ValueError("cannot set negative valuation element to Teichmuller representative.")
        elif self.relprec == 0:
            raise ValueError("not enough precision")
        else:
            cteichmuller(self.unit, self.unit, self.relprec, self.prime_pow)

    def precision_absolute(self):
        """
        Returns the absolute precision of this element.

        This is the power of the maximal ideal modulo which this
        element is defined.

        EXAMPLES::

            sage: R = Zp(7,3,'capped-rel'); a = R(7); a.precision_absolute()
            4
            sage: R = Qp(7,3); a = R(7); a.precision_absolute()
            4
            sage: R(7^-3).precision_absolute()
            0

            sage: R(0).precision_absolute()
            +Infinity
            sage: R(0,7).precision_absolute()
            7
        """
        if exactzero(self.ordp):
            return infinity
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.ordp + self.relprec)
        return ans

    def precision_relative(self):
        """
        Returns the relative precision of this element.

        This is the power of the maximal ideal modulo which the unit
        part of self is defined.

        EXAMPLES::

            sage: R = Zp(7,3,'capped-rel'); a = R(7); a.precision_relative()
            3
            sage: R = Qp(7,3); a = R(7); a.precision_relative()
            3
            sage: a = R(7^-2, -1); a.precision_relative()
            1
            sage: a
            7^-2 + O(7^-1)

            sage: R(0).precision_relative()
            0
            sage: R(0,7).precision_relative()
            0
       """
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.relprec)
        return ans

    cpdef pAdicTemplateElement unit_part(self):
        r"""
        Returns `u`, where this element is `\pi^v u`.

        EXAMPLES::

            sage: R = Zp(17,4,'capped-rel')
            sage: a = R(18*17)
            sage: a.unit_part()
            1 + 17 + O(17^4)
            sage: type(a)
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>
            sage: R = Qp(17,4,'capped-rel')
            sage: a = R(18*17)
            sage: a.unit_part()
            1 + 17 + O(17^4)
            sage: type(a)
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>
            sage: a = R(2*17^2); a
            2*17^2 + O(17^6)
            sage: a.unit_part()
            2 + O(17^4)
            sage: b=1/a; b
            9*17^-2 + 8*17^-1 + 8 + 8*17 + O(17^2)
            sage: b.unit_part()
            9 + 8*17 + 8*17^2 + 8*17^3 + O(17^4)
            sage: Zp(5)(75).unit_part()
            3 + O(5^20)

            sage: R(0).unit_part()
            Traceback (most recent call last):
            ...
            ValueError: unit part of 0 not defined
            sage: R(0,7).unit_part()
            O(17^0)
        """
        if exactzero(self.ordp):
            raise ValueError("unit part of 0 not defined")
        cdef CRElement ans = (<CRElement>self)._new_c()
        ans.ordp = 0
        ans.relprec = (<CRElement>self).relprec
        ccopy(ans.unit, (<CRElement>self).unit, ans.prime_pow)
        return ans

    cdef long valuation_c(self):
        """
        Returns the valuation of this element.

        If self is an exact zero, returns ``maxordp``, which is defined as
        ``(1L << (sizeof(long) * 8 - 2))-1``.

        EXAMPLES::

            sage: R = Qp(5); a = R(1)
            sage: a.valuation() #indirect doctest
            0
            sage: b = (a << 4); b.valuation()
            4
            sage: b = (a << 1073741822); b.valuation()
            1073741822
        """
        return self.ordp

    cpdef val_unit(self, p=None):
        """
        Returns a pair ``(self.valuation(), self.unit_part())``.

        INPUT:

        - ``p`` -- a prime (default: ``None``). If specified, will make sure that p==self.parent().prime()

        .. NOTE::

            The optional argument ``p`` is used for consistency with the
            valuation methods on integer and rational.

        EXAMPLES::

            sage: R = Zp(5); a = R(75, 20); a
            3*5^2 + O(5^20)
            sage: a.val_unit()
            (2, 3 + O(5^18))
            sage: R(0).val_unit()
            Traceback (most recent call last):
            ...
            ValueError: unit part of 0 not defined
            sage: R(0, 10).val_unit()
            (10, O(5^0))
        """
        # Since we keep this element normalized there's not much to do here.
        if p is not None and p != self.parent().prime():
            raise ValueError('Ring (%s) residue field of the wrong characteristic.'%self.parent())
        if exactzero((<CRElement>self).ordp):
            raise ValueError("unit part of 0 not defined")
        cdef Integer val = PY_NEW(Integer)
        mpz_set_si(val.value, (<CRElement>self).ordp)
        cdef CRElement unit = (<CRElement>self)._new_c()
        unit.ordp = 0
        unit.relprec = (<CRElement>self).relprec
        ccopy(unit.unit, (<CRElement>self).unit, unit.prime_pow)
        return val, unit

    def __hash__(self):
        """
        Hashing.

        .. WARNING::

            Hashing of `p`-adic elements will likely be deprecated soon.  See :trac:`11895`.

        EXAMPLES::

            sage: R = Zp(5)
            sage: hash(R(17)) #indirect doctest
            17

            sage: hash(R(-1))
            1977822444 # 32-bit
            95367431640624 # 64-bit
        """
        if exactzero(self.ordp):
            return 0
        return chash(self.unit, self.ordp, self.relprec, self.prime_pow) ^ self.ordp

cdef class pAdicCoercion_ZZ_CR(RingHomomorphism_coercion):
    """
    The canonical inclusion from the integer ring to a capped relative ring.

    EXAMPLES::

        sage: f = Zp(5).coerce_map_from(ZZ); f
        Ring Coercion morphism:
          From: Integer Ring
          To:   5-adic Ring with capped relative precision 20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = Zp(5).coerce_map_from(ZZ); type(f)
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicCoercion_ZZ_CR'>
        """
        RingHomomorphism_coercion.__init__(self, ZZ.Hom(R), check=False)
        self._zero = <CRElement?>R._element_constructor(R, 0)
        self._section = pAdicConvert_CR_ZZ(R)

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Zp(5).coerce_map_from(ZZ)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring Coercion morphism:
              From: Integer Ring
              To:   5-adic Ring with capped relative precision 20
            sage: g == f
            True
            sage: g is f
            False
            sage: g(5)
            5 + O(5^21)
            sage: g(5) == f(5)
            True

        """
        _slots['_zero'] = self._zero
        _slots['_section'] = self._section
        return RingHomomorphism_coercion._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Zp(5).coerce_map_from(ZZ)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring Coercion morphism:
              From: Integer Ring
              To:   5-adic Ring with capped relative precision 20
            sage: g == f
            True
            sage: g is f
            False
            sage: g(5)
            5 + O(5^21)
            sage: g(5) == f(5)
            True

        """
        self._zero = _slots['_zero']
        self._section = _slots['_section']
        RingHomomorphism_coercion._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = Zp(5).coerce_map_from(ZZ)
            sage: f(0).parent()
            5-adic Ring with capped relative precision 20
            sage: f(5)
            5 + O(5^21)
        """
        if mpz_sgn((<Integer>x).value) == 0:
            return self._zero
        cdef CRElement ans = self._zero._new_c()
        ans.relprec = ans.prime_pow.prec_cap
        ans.ordp = cconv_mpz_t(ans.unit, (<Integer>x).value, ans.relprec, False, ans.prime_pow)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both), or an empty element is
        desired.

        See the documentation for
        :meth:`pAdicCappedRelativeElement.__init__` for more details.

        EXAMPLES::

            sage: R = Zp(5,4)
            sage: type(R(10,2))
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>
            sage: R(10,2)
            2*5 + O(5^2)
            sage: R(10,3,1)
            2*5 + O(5^2)
            sage: R(10,absprec=2)
            2*5 + O(5^2)
            sage: R(10,relprec=2)
            2*5 + O(5^3)
            sage: R(10,absprec=1)
            O(5)
            sage: R(10,empty=True)
            O(5^0)
        """
        cdef long val, aprec, rprec
        cdef CRElement ans
        _process_args_and_kwds(&aprec, &rprec, args, kwds, False, self._zero.prime_pow)
        if mpz_sgn((<Integer>x).value) == 0:
            if exactzero(aprec):
                return self._zero
            ans = self._zero._new_c()
            ans._set_inexact_zero(aprec)
        else:
            val = get_ordp(x, self._zero.prime_pow)
            ans = self._zero._new_c()
            if aprec <= val:
                ans._set_inexact_zero(aprec)
            else:
                ans.relprec = min(rprec, aprec - val)
                ans.ordp = cconv_mpz_t(ans.unit, (<Integer>x).value, ans.relprec, False, self._zero.prime_pow)
        return ans

    def section(self):
        """
        Returns a map back to the ring of integers that approximates an element
        by an integer.

        EXAMPLES::

            sage: f = Zp(5).coerce_map_from(ZZ).section()
            sage: f(Zp(5)(-1)) - 5^20
            -1
        """
        return self._section

cdef class pAdicConvert_CR_ZZ(RingMap):
    """
    The map from a capped relative ring back to the ring of integers that
    returns the the smallest non-negative integer approximation to its input
    which is accurate up to the precision.

    Raises a ``ValueError``, if the input is not in the closure of the image of
    the integers.

    EXAMPLES::

        sage: f = Zp(5).coerce_map_from(ZZ).section(); f
        Set-theoretic ring morphism:
          From: 5-adic Ring with capped relative precision 20
          To:   Integer Ring
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(ZZ).section(); type(f)
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicConvert_CR_ZZ'>
            sage: f.category()
            Category of homsets of sets with partial maps
            sage: Zp(5).coerce_map_from(ZZ).section().category()
            Category of homsets of sets
        """
        if R.is_field() or R.degree() > 1 or R.characteristic() != 0 or R.residue_characteristic() == 0:
            RingMap.__init__(self, Hom(R, ZZ, SetsWithPartialMaps()))
        else:
            RingMap.__init__(self, Hom(R, ZZ, Sets()))

    cpdef Element _call_(self, _x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(ZZ).section()
            sage: f(Qp(5)(-1)) - 5^20
            -1
            sage: f(Qp(5)(0))
            0
            sage: f(Qp(5)(1/5))
            Traceback (most recent call last):
            ...
            ValueError: negative valuation
        """
        cdef Integer ans = PY_NEW(Integer)
        cdef CRElement x = _x
        if x.relprec != 0:
            cconv_mpz_t_out(ans.value, x.unit, x.ordp, x.relprec, x.prime_pow)
        return ans

cdef class pAdicCoercion_QQ_CR(RingHomomorphism_coercion):
    """
    The canonical inclusion from the rationals to a capped relative field.

    EXAMPLES::

        sage: f = Qp(5).coerce_map_from(QQ); f
        Ring Coercion morphism:
          From: Rational Field
          To:   5-adic Field with capped relative precision 20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ); type(f)
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicCoercion_QQ_CR'>
        """
        RingHomomorphism_coercion.__init__(self, QQ.Hom(R), check=False)
        self._zero = R._element_constructor(R, 0)
        self._section = pAdicConvert_CR_QQ(R)

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring Coercion morphism:
              From: Rational Field
              To:   5-adic Field with capped relative precision 20
            sage: g == f
            True
            sage: g is f
            False
            sage: g(6)
            1 + 5 + O(5^20)
            sage: g(6) == f(6)
            True

        """
        _slots['_zero'] = self._zero
        _slots['_section'] = self._section
        return RingHomomorphism_coercion._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring Coercion morphism:
              From: Rational Field
              To:   5-adic Field with capped relative precision 20
            sage: g == f
            True
            sage: g is f
            False
            sage: g(6)
            1 + 5 + O(5^20)
            sage: g(6) == f(6)
            True

        """
        self._zero = _slots['_zero']
        self._section = _slots['_section']
        RingHomomorphism_coercion._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ)
            sage: f(0).parent()
            5-adic Field with capped relative precision 20
            sage: f(1/5)
            5^-1 + O(5^19)
            sage: f(1/4)
            4 + 3*5 + 3*5^2 + 3*5^3 + 3*5^4 + 3*5^5 + 3*5^6 + 3*5^7 + 3*5^8 + 3*5^9 + 3*5^10 + 3*5^11 + 3*5^12 + 3*5^13 + 3*5^14 + 3*5^15 + 3*5^16 + 3*5^17 + 3*5^18 + 3*5^19 + O(5^20)
        """
        if mpq_sgn((<Rational>x).value) == 0:
            return self._zero
        cdef CRElement ans = self._zero._new_c()
        ans.relprec = ans.prime_pow.prec_cap
        ans.ordp = cconv_mpq_t(ans.unit, (<Rational>x).value, ans.relprec, False, self._zero.prime_pow)
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both), or an empty element is
        desired.

        See the documentation for
        :meth:`pAdicCappedRelativeElement.__init__` for more details.

        EXAMPLES::

            sage: R = Qp(5,4)
            sage: type(R(10/3,2))
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>
            sage: R(10/3,2)
            4*5 + O(5^2)
            sage: R(10/3,3,1)
            4*5 + O(5^2)
            sage: R(10/3,absprec=2)
            4*5 + O(5^2)
            sage: R(10/3,relprec=2)
            4*5 + 5^2 + O(5^3)
            sage: R(10/3,absprec=1)
            O(5)
            sage: R(10/3,empty=True)
            O(5^0)
            sage: R(3/100,absprec=-1)
            2*5^-2 + O(5^-1)
        """
        cdef long val, aprec, rprec
        cdef CRElement ans
        _process_args_and_kwds(&aprec, &rprec, args, kwds, False, self._zero.prime_pow)
        if mpq_sgn((<Rational>x).value) == 0:
            if exactzero(aprec):
                return self._zero
            ans = self._zero._new_c()
            ans._set_inexact_zero(aprec)
        else:
            val = get_ordp(x, self._zero.prime_pow)
            ans = self._zero._new_c()
            if aprec <= val:
                ans._set_inexact_zero(aprec)
            else:
                ans.relprec = min(rprec, aprec - val)
                ans.ordp = cconv_mpq_t(ans.unit, (<Rational>x).value, ans.relprec, False, self._zero.prime_pow)
        return ans

    def section(self):
        """
        Returns a map back to the rationals that approximates an element by
        a rational number.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ).section()
            sage: f(Qp(5)(1/4))
            1/4
            sage: f(Qp(5)(1/5))
            1/5
        """
        return self._section

cdef class pAdicConvert_CR_QQ(RingMap):
    """
    The map from the capped relative ring back to the rationals that returns a
    rational approximation of its input.

    EXAMPLES::

        sage: f = Qp(5).coerce_map_from(QQ).section(); f
        Set-theoretic ring morphism:
          From: 5-adic Field with capped relative precision 20
          To:   Rational Field
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ).section(); type(f)
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicConvert_CR_QQ'>
            sage: f.category()
            Category of homsets of sets
        """
        if R.degree() > 1 or R.characteristic() != 0 or R.residue_characteristic() == 0:
            RingMap.__init__(self, Hom(R, QQ, SetsWithPartialMaps()))
        else:
            RingMap.__init__(self, Hom(R, QQ, Sets()))

    cpdef Element _call_(self, _x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = Qp(5).coerce_map_from(QQ).section()
            sage: f(Qp(5)(-1))
            -1
            sage: f(Qp(5)(0))
            0
            sage: f(Qp(5)(1/5))
            1/5
        """
        cdef Rational ans = Rational.__new__(Rational)
        cdef CRElement x =  _x
        if x.relprec == 0:
            mpq_set_ui(ans.value, 0, 1)
        else:
            cconv_mpq_t_out(ans.value, x.unit, x.ordp, x.relprec, x.prime_pow)
        return ans

cdef class pAdicConvert_QQ_CR(Morphism):
    """
    The inclusion map from the rationals to a capped relative ring that is
    defined on all elements with non-negative `p`-adic valuation.

    EXAMPLES::

        sage: f = Zp(5).convert_map_from(QQ); f
        Generic morphism:
          From: Rational Field
          To:   5-adic Ring with capped relative precision 20
    """
    def __init__(self, R):
        """
        Initialization.

        EXAMPLES::

            sage: f = Zp(5).convert_map_from(QQ); type(f)
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicConvert_QQ_CR'>
        """
        Morphism.__init__(self, Hom(QQ, R, SetsWithPartialMaps()))
        self._zero = R._element_constructor(R, 0)
        self._section = pAdicConvert_CR_QQ(R)

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Zp(5).convert_map_from(QQ)
            sage: g = copy(f)   # indirect doctest
            sage: g == f
            True
            sage: g(1/6)
            1 + 4*5 + 4*5^3 + 4*5^5 + 4*5^7 + 4*5^9 + 4*5^11 + 4*5^13 + 4*5^15 + 4*5^17 + 4*5^19 + O(5^20)
            sage: g(1/6) == f(1/6)
            True
        """
        _slots['_zero'] = self._zero
        _slots['_section'] = self._section
        return Morphism._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: f = Zp(5).convert_map_from(QQ)
            sage: g = copy(f)   # indirect doctest
            sage: g == f
            True
            sage: g(1/6)
            1 + 4*5 + 4*5^3 + 4*5^5 + 4*5^7 + 4*5^9 + 4*5^11 + 4*5^13 + 4*5^15 + 4*5^17 + 4*5^19 + O(5^20)
            sage: g(1/6) == f(1/6)
            True
        """
        self._zero = _slots['_zero']
        self._section = _slots['_section']
        Morphism._update_slots(self, _slots)

    cpdef Element _call_(self, x):
        """
        Evaluation.

        EXAMPLES::

            sage: f = Zp(5,4).convert_map_from(QQ)
            sage: f(1/7)
            3 + 3*5 + 2*5^3 + O(5^4)
            sage: f(0)
            0
        """
        if mpq_sgn((<Rational>x).value) == 0:
            return self._zero
        cdef CRElement ans = self._zero._new_c()
        ans.relprec = ans.prime_pow.prec_cap
        ans.ordp = cconv_mpq_t(ans.unit, (<Rational>x).value, ans.relprec, False, self._zero.prime_pow)
        if ans.ordp < 0:
            raise ValueError("p divides the denominator")
        return ans

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both), or an empty element is
        desired.

        See the documentation for
        :meth:`pAdicCappedRelativeElement.__init__` for more details.

        EXAMPLES::

            sage: R = Zp(5,4)
            sage: type(R(10/3,2))
            <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>
            sage: R(10/3,2)
            4*5 + O(5^2)
            sage: R(10/3,3,1)
            4*5 + O(5^2)
            sage: R(10/3,absprec=2)
            4*5 + O(5^2)
            sage: R(10/3,relprec=2)
            4*5 + 5^2 + O(5^3)
            sage: R(10/3,absprec=1)
            O(5)
            sage: R(10/3,empty=True)
            O(5^0)
            sage: R(3/100,relprec=3)
            Traceback (most recent call last):
            ...
            ValueError: p divides the denominator
        """
        cdef long val, aprec, rprec
        cdef CRElement ans
        _process_args_and_kwds(&aprec, &rprec, args, kwds, False, self._zero.prime_pow)
        if mpq_sgn((<Rational>x).value) == 0:
            if exactzero(aprec):
                return self._zero
            ans = self._zero._new_c()
            ans._set_inexact_zero(aprec)
        else:
            val = get_ordp(x, self._zero.prime_pow)
            ans = self._zero._new_c()
            if aprec <= val:
                ans._set_inexact_zero(aprec)
            else:
                ans.relprec = min(rprec, aprec - val)
                ans.ordp = cconv_mpq_t(ans.unit, (<Rational>x).value, ans.relprec, False, self._zero.prime_pow)
        if ans.ordp < 0:
            raise ValueError("p divides the denominator")
        return ans

    def section(self):
        """
        Returns the map back to the rationals that returns the smallest
        non-negative integer approximation to its input which is accurate up to
        the precision.

        EXAMPLES::

            sage: f = Zp(5,4).convert_map_from(QQ).section()
            sage: f(Zp(5,4)(-1))
            -1
        """
        return self._section

cdef class pAdicCoercion_CR_frac_field(RingHomomorphism_coercion):
    """
    The canonical inclusion of Zq into its fraction field.

    EXAMPLES::

        sage: R.<a> = ZqCR(27, implementation='FLINT')
        sage: K = R.fraction_field()
        sage: K.coerce_map_from(R)
        Ring Coercion morphism:
          From: Unramified Extension of 3-adic Ring with capped relative precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
          To:   Unramified Extension of 3-adic Field with capped relative precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
    """
    def __init__(self, R, K):
        """
        Initialization.

        EXAMPLES::

            sage: R.<a> = ZqCR(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R); type(f)
            <type 'sage.rings.padics.qadic_flint_CR.pAdicCoercion_CR_frac_field'>
        """
        RingHomomorphism_coercion.__init__(self, R.Hom(K), check=False)
        self._zero = K(0)
        self._section = pAdicConvert_CR_frac_field(K, R)

    cpdef Element _call_(self, _x):
        """
        Evaluation.

        EXAMPLES::

            sage: R.<a> = ZqCR(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(a)
            a + O(3^20)
            sage: f(R(0))
            0
        """
        cdef CRElement x = _x
        cdef CRElement ans = self._zero._new_c()
        ans.ordp = x.ordp
        ans.relprec = x.relprec
        cshift(ans.unit, x.unit, 0, ans.relprec, x.prime_pow, False)
        return ans

    cpdef Element _call_with_args(self, _x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both).

        See the documentation for
        :meth:`pAdicCappedAbsoluteElement.__init__` for more details.

        EXAMPLES::

            sage: R.<a> = ZqCR(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(a, 3)
            a + O(3^3)
            sage: b = 9*a
            sage: f(b, 3)
            a*3^2 + O(3^3)
            sage: f(b, 4, 1)
            a*3^2 + O(3^3)
            sage: f(b, 4, 3)
            a*3^2 + O(3^4)
            sage: f(b, absprec=4)
            a*3^2 + O(3^4)
            sage: f(b, relprec=3)
            a*3^2 + O(3^5)
            sage: f(b, absprec=1)
            O(3)
            sage: f(R(0))
            0
        """
        cdef long aprec, rprec
        cdef CRElement x = _x
        cdef CRElement ans = self._zero._new_c()
        cdef bint reduce = False
        _process_args_and_kwds(&aprec, &rprec, args, kwds, False, ans.prime_pow)
        if aprec <= x.ordp:
            csetzero(ans.unit, x.prime_pow)
            ans.relprec = 0
            ans.ordp = aprec
        else:
            if rprec < x.relprec:
                reduce = True
            else:
                rprec = x.relprec
            if aprec < rprec + x.ordp:
                rprec = aprec - x.ordp
                reduce = True
            ans.ordp = x.ordp
            ans.relprec = rprec
            cshift(ans.unit, x.unit, 0, rprec, x.prime_pow, reduce)
        return ans

    def section(self):
        """
        Returns a map back to the ring that converts elements of
        non-negative valuation.

        EXAMPLES::

            sage: R.<a> = ZqCR(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: f(K.gen())
            a + O(3^20)
        """
        return self._section

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        TESTS::

            sage: R.<a> = ZqCR(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring Coercion morphism:
              From: Unramified Extension of 3-adic Ring with capped relative precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
              To:   Unramified Extension of 3-adic Field with capped relative precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
            sage: g == f
            True
            sage: g is f
            False
            sage: g(a)
            a + O(3^20)
            sage: g(a) == f(a)
            True

        """
        _slots['_zero'] = self._zero
        _slots['_section'] = self._section
        return RingHomomorphism_coercion._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        TESTS::

            sage: R.<a> = ZqCR(9, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = K.coerce_map_from(R)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Ring Coercion morphism:
              From: Unramified Extension of 3-adic Ring with capped relative precision 20 in a defined by (1 + O(3^20))*x^2 + (2 + O(3^20))*x + (2 + O(3^20))
              To:   Unramified Extension of 3-adic Field with capped relative precision 20 in a defined by (1 + O(3^20))*x^2 + (2 + O(3^20))*x + (2 + O(3^20))
            sage: g == f
            True
            sage: g is f
            False
            sage: g(a)
            a + O(3^20)
            sage: g(a) == f(a)
            True

        """
        self._zero = _slots['_zero']
        self._section = _slots['_section']
        RingHomomorphism_coercion._update_slots(self, _slots)

cdef class pAdicConvert_CR_frac_field(Morphism):
    """
    The section of the inclusion from `\ZZ_q`` to its fraction field.

    EXAMPLES::

        sage: R.<a> = ZqCR(27, implementation='FLINT')
        sage: K = R.fraction_field()
        sage: f = R.convert_map_from(K); f
        Generic morphism:
          From: Unramified Extension of 3-adic Field with capped relative precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
          To:   Unramified Extension of 3-adic Ring with capped relative precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
    """
    def __init__(self, K, R):
        """
        Initialization.

        EXAMPLES::

            sage: R.<a> = ZqCR(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K); type(f)
            <type 'sage.rings.padics.qadic_flint_CR.pAdicConvert_CR_frac_field'>
        """
        Morphism.__init__(self, Hom(K, R, SetsWithPartialMaps()))
        self._zero = R(0)

    cpdef Element _call_(self, _x):
        """
        Evaluation.

        EXAMPLES::

            sage: R.<a> = ZqCR(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K)
            sage: f(K.gen())
            a + O(3^20)
        """
        cdef CRElement x = _x
        if x.ordp < 0: raise ValueError("negative valuation")
        cdef CRElement ans = self._zero._new_c()
        ans.relprec = x.relprec
        ans.ordp = x.ordp
        cshift(ans.unit, x.unit, 0, ans.relprec, ans.prime_pow, False)
        return ans

    cpdef Element _call_with_args(self, _x, args=(), kwds={}):
        """
        This function is used when some precision cap is passed in
        (relative or absolute or both).

        See the documentation for
        :meth:`pAdicCappedAbsoluteElement.__init__` for more details.

        EXAMPLES::

            sage: R.<a> = ZqCR(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K); a = K(a)
            sage: f(a, 3)
            a + O(3^3)
            sage: b = 9*a
            sage: f(b, 3)
            a*3^2 + O(3^3)
            sage: f(b, 4, 1)
            a*3^2 + O(3^3)
            sage: f(b, 4, 3)
            a*3^2 + O(3^4)
            sage: f(b, absprec=4)
            a*3^2 + O(3^4)
            sage: f(b, relprec=3)
            a*3^2 + O(3^5)
            sage: f(b, absprec=1)
            O(3)
            sage: f(K(0))
            0
        """
        cdef long aprec, rprec
        cdef CRElement x = _x
        if x.ordp < 0: raise ValueError("negative valuation")
        cdef CRElement ans = self._zero._new_c()
        cdef bint reduce = False
        _process_args_and_kwds(&aprec, &rprec, args, kwds, False, ans.prime_pow)
        if aprec <= x.ordp:
            csetzero(ans.unit, x.prime_pow)
            ans.relprec = 0
            ans.ordp = aprec
        else:
            if rprec < x.relprec:
                reduce = True
            else:
                rprec = x.relprec
            if aprec < rprec + x.ordp:
                rprec = aprec - x.ordp
                reduce = True
            ans.ordp = x.ordp
            ans.relprec = rprec
            cshift(ans.unit, x.unit, 0, rprec, x.prime_pow, reduce)
        return ans

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        TESTS::

            sage: R.<a> = ZqCR(27, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K)
            sage: a = K(a)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Generic morphism:
              From: Unramified Extension of 3-adic Field with capped relative precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
              To:   Unramified Extension of 3-adic Ring with capped relative precision 20 in a defined by (1 + O(3^20))*x^3 + (O(3^20))*x^2 + (2 + O(3^20))*x + (1 + O(3^20))
            sage: g == f
            True
            sage: g is f
            False
            sage: g(a)
            a + O(3^20)
            sage: g(a) == f(a)
            True

        """
        _slots['_zero'] = self._zero
        return Morphism._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        TESTS::

            sage: R.<a> = ZqCR(9, implementation='FLINT')
            sage: K = R.fraction_field()
            sage: f = R.convert_map_from(K)
            sage: a = K(a)
            sage: g = copy(f)   # indirect doctest
            sage: g
            Generic morphism:
              From: Unramified Extension of 3-adic Field with capped relative precision 20 in a defined by (1 + O(3^20))*x^2 + (2 + O(3^20))*x + (2 + O(3^20))
              To:   Unramified Extension of 3-adic Ring with capped relative precision 20 in a defined by (1 + O(3^20))*x^2 + (2 + O(3^20))*x + (2 + O(3^20))
            sage: g == f
            True
            sage: g is f
            False
            sage: g(a)
            a + O(3^20)
            sage: g(a) == f(a)
            True

        """
        self._zero = _slots['_zero']
        Morphism._update_slots(self, _slots)

def unpickle_cre_v2(cls, parent, unit, ordp, relprec):
    """
    Unpickles a capped relative element.

    EXAMPLES::

        sage: from sage.rings.padics.padic_capped_relative_element import unpickle_cre_v2
        sage: R = Zp(5); a = R(85,6)
        sage: b = unpickle_cre_v2(a.__class__, R, 17, 1, 5)
        sage: a == b
        True
        sage: a.precision_relative() == b.precision_relative()
        True
    """
    cdef CRElement ans = cls.__new__(cls)
    ans._parent = parent
    ans.prime_pow = <PowComputer_?>parent.prime_pow
    cconstruct(ans.unit, ans.prime_pow)
    cunpickle(ans.unit, unit, ans.prime_pow)
    ans.ordp = ordp
    ans.relprec = relprec
    return ans
