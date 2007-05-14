"""
Elements of p-Adic Rings with Capped Relative Precision

AUTHOR:
    -- David Roe
    -- Genya Zaytman: documentation
    -- David Harvey: doctests
"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../../ext/gmp.pxi"
include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"

cimport sage.rings.padics.padic_generic_element
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement

cimport sage.rings.rational
cimport sage.rings.padics.pow_computer
from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.rational cimport Rational

import sage.rings.padics.padic_generic_element
import sage.rings.padics.padic_lazy_element
import sage.rings.integer_mod
import sage.rings.integer
import sage.rings.rational

from sage.rings.infinity import infinity
from sage.rings.integer_mod import Mod
from sage.rings.padics.precision_error import PrecisionError

from sage.rings.padics.padic_lazy_element import pAdicLazyElement

pari = sage.libs.pari.gen.pari
pari_gen = sage.libs.pari.gen.gen
PariError = sage.libs.pari.gen.PariError

cdef class pAdicCappedRelativeElement(pAdicBaseGenericElement):
    def __init__(pAdicCappedRelativeElement self, parent, x, absprec=infinity, relprec=infinity, empty = False):
        """
        Constructs new element with given parent and value.

        INPUT:
            x -- value to coerce into a capped relative ring or field
            absprec -- maximum number of digits of absolute precision
            relprec -- maximum number of digits of relative precision
            construct -- boolean, default False. True is for internal use,
                in which case x is a triple to be assigned directly.

        EXAMPLES:
            sage: R = Zp(5, 10, 'capped-rel')

        Construct from integers:
            sage: R(3)
            3 + O(5^10)
            sage: R(75)
            3*5^2 + O(5^12)
            sage: R(0)
            0
            sage: R(-1)
            4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)
            sage: R(-5)
            4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + O(5^11)
            sage: R(-7*25)
            3*5^2 + 3*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + O(5^12)

        Construct from rationals:
            sage: R(1/2)
            3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 2*5^8 + 2*5^9 + O(5^10)
            sage: R(-7875/874)
            3*5^3 + 2*5^4 + 2*5^5 + 5^6 + 3*5^7 + 2*5^8 + 3*5^10 + 3*5^11 + 3*5^12 + O(5^13)
            sage: R(15/425)
            Traceback (most recent call last):
            ...
            ValueError: p divides the denominator

        Construct from IntegerMod:
            sage: R(Integers(125)(3))
            3 + O(5^3)
            sage: R(Integers(5)(3))
            3 + O(5)
            sage: R(Integers(5^30)(3))
            3 + O(5^10)
            sage: R(Integers(5^30)(1+5^23))
            1 + O(5^10)
            sage: R(Integers(49)(3))
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce from the given integer mod ring (not a power of the same prime)

        # todo: should the above TypeError be another type of error?

            sage: R(Integers(48)(3))
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce from the given integer mod ring (not a power of the same prime)

        # todo: the error message for the above TypeError is not quite accurate

        Some other conversions:
            sage: R(R(5))
            5 + O(5^11)

        # todo: doctests for converting from other types of p-adic rings

        """
        cdef RingElement ordp
        mpz_init(self.unit)
        pAdicGenericElement.__init__(self, parent)
        if empty:
            self._normalized = 0
            return
        self._normalized = 1
        if relprec > parent.precision_cap():
            relprec = parent.precision_cap()
        if not absprec is infinity and not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)

        if isinstance(x, pAdicGenericElement):
            if self.prime_pow.in_field == 0 and x.valuation() < 0:
                raise ValueError, "element has negative valuation."
            if parent.prime() != x.parent().prime():
                raise TypeError, "Cannot coerce between p-adic parents with different primes"
        if isinstance(x, pAdicLazyElement):
            # One can do this in a better way to minimize the amount of increasing precision on x.
            if absprec is infinity:
                try:
                    x.set_precision_relative(relprec)
                except PrecisionError:
                    pass
            else:
                try:
                    x.set_precision_absolute(absprec)
                except PrecisionError:
                    pass
            if x.precision_relative() < relprec:
                try:
                    x.set_precision_relative(relprec)
                except PrecisionError:
                    pass
        if isinstance(x, pAdicBaseGenericElement):
            ordp = x.valuation()
            if ordp >= absprec:
                if absprec is infinity:
                    self.set_exact_zero()
                else:
                    self.set_inexact_zero(mpz_get_si((<Integer>absprec).value))
            else:
                unit = x.unit_part().lift()
                relprec = min(relprec, absprec - ordp, x.precision_relative())
                self.set_from_Integers(ordp, unit, relprec)
            return

        if isinstance(x, pari_gen):
            if x.type() == "t_PADIC":
                absprec = min(x.padicprec(parent.prime()), absprec)
                x = x.lift()
            if x.type() == "t_INT":
                x = Integer(x)
            elif x.type() == "t_FRAC":
                x = Rational(x)
            else:
                raise TypeError, "unsupported coercion from pari: only p-adics, integers and rationals allowed"

        #if sage.rings.finite_field_element.is_FiniteFieldElement(x):
        #    if x.parent().order() != parent.prime():
        #        raise TypeError, "can only create p-adic element out of finite field when order of field is p"
        #    #prec = min(prec, 1)
        #    x = x.lift()

        cdef mpz_t modulus
        cdef Integer tmp
        cdef int k
        if sage.rings.integer_mod.is_IntegerMod(x):
            mpz_init_set(modulus, (<Integer>x.modulus()).value)
            k = mpz_remove(modulus, modulus, self.prime_pow.prime.value)
            if mpz_cmp_ui(modulus, 1) == 0:
                absprec = min(k, absprec)
                x = x.lift()
                mpz_clear(modulus)
            else:
                mpz_clear(modulus)
                raise TypeError, "cannot coerce from the given integer mod ring (not a power of the same prime)"

            # We now use the code, below, so don't make the next line elif
        if isinstance(x, (int, long)):
            x = Integer(x)
        #print "before integer"
        if isinstance(x, Integer):
            ordp, unit = x.val_unit(parent.prime())
            #print "after pair"
            if ordp >= absprec:
                if absprec is infinity:
                    self.set_exact_zero()
                else:
                    self.set_inexact_zero(mpz_get_si((<Integer>absprec).value))
                return
            #print "now 1"
            relprec = min(relprec, absprec - ordp)
            #print "now 2"
            self.set_from_Integers(ordp, unit, relprec)
            #print "now 3"
        elif isinstance(x, Rational):
            ordp, unit = x.val_unit(parent.prime())
            if ordp >= absprec:
                if absprec is infinity:
                    self.set_exact_zero()
                else:
                    self.set_inexact_zero(mpz_get_si((<Integer>absprec).value))
                return
            if self.prime_pow.in_field == 0 and ordp < 0:
                raise ValueError, "p divides the denominator"
            # todo: make the following faster
            relprec = min(relprec, absprec - ordp)
            unit = Mod(unit, parent.prime_pow(relprec)).lift()
            self.set_from_Integers(ordp, unit, relprec)
        else:
            raise TypeError, "cannot create a p-adic out of %s"%(type(x))
        if self.ordp < 0 and self.prime_pow.in_field == 0:
            raise ValueError, "element has negative valuation."

    cdef void set_exact_zero(pAdicCappedRelativeElement self):
        mpz_set_si(self.unit, -1)
        self._normalized = 1

    cdef void set_inexact_zero(pAdicCappedRelativeElement self, long absprec):
        self.relprec = 0
        mpz_set_ui(self.unit, 0)
        self.ordp = absprec
        self._normalized = 1

    cdef void set_precs(pAdicCappedRelativeElement self, unsigned long relprec):
        self.relprec = relprec

    cdef void set_from_Integers(pAdicCappedRelativeElement self, Integer ordp, Integer unit, Integer relprec):
        self.set_precs(mpz_get_ui(relprec.value))
        self.ordp = mpz_get_si(ordp.value)
        mpz_set(self.unit, unit.value)
        if mpz_sgn(self.unit) == -1 or mpz_cmp(self.unit, self.prime_pow.dense_list[self.relprec]) >= 0:
            mpz_mod(self.unit, self.unit, self.prime_pow.dense_list[self.relprec])

    cdef pAdicCappedRelativeElement _new_c(pAdicCappedRelativeElement self):
        cdef pAdicCappedRelativeElement ans
        ans = PY_NEW(pAdicCappedRelativeElement)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        mpz_init(ans.unit)
        return ans

    cdef void _normalize(pAdicCappedRelativeElement self):
        cdef long diff
        if self._normalized == 0:
            if mpz_sgn(self.unit) > 0:
                if mpz_divisible_p(self.unit, self.prime_pow.prime.value) != 0:
                    diff = mpz_remove(self.unit, self.unit, self.prime_pow.prime.value)
                    if self.relprec > diff:
                        self.ordp += diff
                        self.relprec -= diff
                    else:
                        self.set_inexact_zero(self.ordp + self.relprec)
                if mpz_cmp(self.unit, self.prime_pow.dense_list[self.relprec]) >= 0:
                    mpz_mod(self.unit, self.unit, self.prime_pow.dense_list[self.relprec])
            elif mpz_sgn(self.unit) == 0:
                self.ordp = self.ordp + self.relprec
                self.relprec = 0
        self._normalized = 1

    def __dealloc__(pAdicCappedRelativeElement self):
        mpz_clear(self.unit)

    def __reduce__(self):
        """
        sage: a = ZpCR(5)(-3)
        sage: type(a)
        <type 'sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement'>
        sage: loads(dumps(a)) == a
        True
        """
        # Necessary for pickling.  See integer.pyx for more info.
        cdef Integer relprec, unit, ordp
        relprec = PY_NEW(Integer)
        mpz_set_ui(relprec.value, self.relprec)
        unit = PY_NEW(Integer)
        mpz_set(unit.value, self.unit)
        ordp = PY_NEW(Integer)
        mpz_set_si(ordp.value, self.ordp)
        return unpickle_pcre_v1, (self.parent(), unit, ordp, relprec)

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cdef ModuleElement _neg_c_impl(self):
        """
        EXAMPLES:
            sage: R = Zp(5, 20, 'capped-rel', 'val-unit')
            sage: -R(1)
            95367431640624 + O(5^20)
            sage: -R(5)
            5 * 95367431640624 + O(5^21)
            sage: -R(0)
            0
        """
        cdef pAdicCappedRelativeElement ans
        self._normalize()
        if mpz_sgn(self.unit) <= 0: # an exact or inexact zero
            return self
        ans = self._new_c()
        ans._normalized = self._normalized
        ans.relprec = self.relprec
        ans.ordp = self.ordp
        mpz_sub(ans.unit, ans.prime_pow.dense_list[self.relprec], self.unit)
        return ans

    def __pow__(pAdicCappedRelativeElement self, right, dummy):
        """
        EXAMPLES:
            sage: R = Zp(19, 5, 'capped-rel','series')
            sage: a = R(-1); a
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: a^2
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
        """
        cdef pAdicCappedRelativeElement ans
        cdef mpz_t tmp
        if not PY_TYPE_CHECK(right, Integer):
            right = Integer(right)
        # if right < 0, we return (~self)^(-right)
        if mpz_sgn((<Integer>right).value) == -1:
            return self._invert_c_impl().__pow__(-right)
        # We now check to see if right is the Integer 0.  In that case, we always return 1 to maximum precision.
        elif mpz_sgn((<Integer>right).value) == 0:
            ans = self._new_c()
            # note that right = Integer(0)
            ans.set_from_Integers(right, Integer(1), self.parent().precision_cap())
            ans._normalized = 1
            return ans
        # We next check if self is an exact zero
        elif mpz_sgn((<pAdicCappedRelativeElement>self).unit) == -1:
            return self
        # We next check if self is an inexact zero
        elif mpz_sgn((<pAdicCappedRelativeElement>self).unit) == 0:
            ans = self._new_c()
            mpz_init_set_si(tmp, (<pAdicCappedRelativeElement>self).ordp)
            mpz_mul(tmp, tmp, (<Integer>right).value)
            if mpz_fits_slong_p(tmp):
                ans.set_inexact_zero(mpz_get_ui(tmp))
            else:
                raise ValueError, "Valuation too large"
            mpz_clear(tmp)
            return ans
        else:
            self._normalize()
            ans = self._new_c()
            # You actually get a bit more relative precision when the exponent is divisible by p...
            ans.set_precs(self.relprec)
            mpz_init_set_si(tmp, (<pAdicCappedRelativeElement>self).ordp)
            mpz_mul(tmp, (<Integer>right).value, tmp)
            if mpz_fits_slong_p(tmp) == 0:
                raise ValueError, "Valuation too large"
            ans.ordp = mpz_get_si(tmp)
            _sig_on
            mpz_powm(ans.unit, (<pAdicCappedRelativeElement>self).unit, (<Integer>right).value, ans.prime_pow.dense_list[ans.relprec])
            _sig_off
            return ans

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        # Once the code for _add_c_impl has stabilized, it may be worth getting rid of the extra function call.
        return self._add_c(right._neg_c())

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        """
        EXAMPLES:
            sage: R = Zp(19, 5, 'capped-rel','series')
            sage: a = R(-1); a
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: b=R(-5/2); b
            7 + 9*19 + 9*19^2 + 9*19^3 + 9*19^4 + O(19^5)
            sage: a+b
            6 + 9*19 + 9*19^2 + 9*19^3 + 9*19^4 + O(19^5)
        """
        cdef pAdicCappedRelativeElement ans
        cdef long tmpL
        # If either self or right is an exact zero, we do the appropriate thing
        if mpz_sgn(self.unit) == -1:
            return right
        if mpz_sgn((<pAdicCappedRelativeElement>right).unit) == -1:
            return self
        if self.ordp == (<pAdicCappedRelativeElement>right).ordp:
            ans = self._new_c()
            # The relative precision of the sum is the minimum of the relative precisions in this case, possibly decreasing if we got cancellation
            if self.relprec < (<pAdicCappedRelativeElement>right).relprec:
                ans.relprec = self.relprec
            else:
                ans.relprec = (<pAdicCappedRelativeElement>right).relprec
            # Since the valuations are the same, we can just add the units
            mpz_add(ans.unit, self.unit, (<pAdicCappedRelativeElement>right).unit)
            ans.ordp = self.ordp
            ans._normalized = 0
        elif self.ordp < (<pAdicCappedRelativeElement>right).ordp:
            tmpL = (<pAdicCappedRelativeElement>right).ordp - self.ordp
            if tmpL >= self.relprec:
                return self
            ans = self._new_c()
            ans.ordp = self.ordp
            if mpz_size(self.prime_pow.modulus) > 10000:
                # We only enable the signal handler if the product will take a while.
                _sig_on
                mpz_mul(ans.unit, (<pAdicCappedRelativeElement>right).unit, self.prime_pow.dense_list[tmpL])
                _sig_off
            else:
                mpz_mul(ans.unit, (<pAdicCappedRelativeElement>right).unit, self.prime_pow.dense_list[tmpL])
            mpz_add(ans.unit, ans.unit, self.unit)
            if self.relprec <= tmpL + (<pAdicCappedRelativeElement>right).relprec:
                ans.set_precs(self.relprec)
            else:
                ans.set_precs(tmpL + (<pAdicCappedRelativeElement>right).relprec)
            ans._normalized = 0
        else:
            tmpL = self.ordp - (<pAdicCappedRelativeElement>right).ordp
            if tmpL >= (<pAdicCappedRelativeElement>right).relprec:
                return right
            ans = self._new_c()
            ans.ordp = (<pAdicCappedRelativeElement>right).ordp
            if mpz_size(self.prime_pow.modulus) > 10000:
                # We only enable the signal handler if the product will take a while.
                _sig_on
                mpz_mul(ans.unit, self.unit, self.prime_pow.dense_list[tmpL])
                _sig_off
            else:
                mpz_mul(ans.unit, self.unit, self.prime_pow.dense_list[tmpL])
            mpz_add(ans.unit, ans.unit, (<pAdicCappedRelativeElement>right).unit)
            if (<pAdicCappedRelativeElement>right).relprec <= tmpL + self.relprec:
                ans.set_precs((<pAdicCappedRelativeElement>right).relprec)
            else:
                ans.set_precs(tmpL + self.relprec)
            ans._normalized = 0
        return ans

    def __invert__(self):
        r"""
        Returns the multiplicative inverse of self.

        EXAMPLES:
            sage: R = Qp(7,4,'capped-rel','series'); a = R(3); a
            3 + O(7^4)
            sage: ~a
            5 + 4*7 + 4*7^2 + 4*7^3 + O(7^4)
        """
        return self._invert_c_impl()

    cdef RingElement _invert_c_impl(self):
        self._normalize()
        if mpz_sgn(self.unit) == -1:
            raise ZeroDivisionError, "cannot divide by zero"
        if self.relprec == 0:
            raise PrecisionError, "cannot divide by something indistinguishable from zero."
        cdef pAdicCappedRelativeElement ans
        ans = self._new_c()
        if ans.prime_pow.in_field == 0:
            ans._parent = self.parent().fraction_field()
            ans.prime_pow = ans._parent.prime_pow
        ans.ordp = -self.ordp
        ans.set_precs(self.relprec)
        _sig_on
        mpz_invert(ans.unit, self.unit, ans.prime_pow.dense_list[ans.relprec])
        _sig_off
        return ans

    def __floordiv__(pAdicCappedRelativeElement self, right):
        """
        If a,b are p-adic integers, then floordiv (//) satisfies the following equation:
        a%b + b*(a//b) == a

        EXAMPLES:
            sage: r = Zp(19)
            sage: a = r(1+19+17*19^3+5*19^4); b = r(19^3); a/b
            19^-3 + 19^-2 + 17 + 5*19 + O(19^17)
            sage: a/b
            19^-3 + 19^-2 + 17 + 5*19 + O(19^17)
            sage: a//b
            17 + 5*19 + O(19^17)

            sage: R = Zp(19, 5, 'capped-rel','series')
            sage: a = R(-1); a
            18 + 18*19 + 18*19^2 + 18*19^3 + 18*19^4 + O(19^5)
            sage: b=R(-2*19^3); b
            17*19^3 + 18*19^4 + 18*19^5 + 18*19^6 + 18*19^7 + O(19^8)
            sage: a//b
            9 + 9*19 + O(19^2)
        """
        if self.parent() is not right.parent():
            right = pAdicCappedRelativeElement(self.parent(), right)
        return self._floordiv_c_impl(right)

    cdef RingElement _floordiv_c_impl(self, RingElement right):
        cdef pAdicCappedRelativeElement ans
        cdef unsigned long relprec, diff
        (<pAdicCappedRelativeElement>right)._normalize()
        self._normalize()
        # For fields, we define floor division as normal division and % as always 0
        if self.prime_pow.in_field == 1:
            return self._div_c(right)
        # We check to see if right is an exact or inexact zero.
        if mpz_sgn((<pAdicCappedRelativeElement>right).unit) == -1:
            raise ZeroDivisionError, "cannot divide by zero"
        elif mpz_sgn((<pAdicCappedRelativeElement>right).unit) == 0:
            raise PrecisionError, "cannot divide by something indistinguishable from zero"
        elif mpz_sgn(self.unit) == -1:
            # If the numerator is an exact zero, we can now just return it.
            return self
        ans = self._new_c()
        # We compute the quotient of the unit parts.
        # One might want to compute this inverse to the precision it's actually needed to for speed's sake...
        _sig_on
        mpz_invert(ans.unit, (<pAdicCappedRelativeElement>right).unit, self.prime_pow.dense_list[(<pAdicCappedRelativeElement>right).relprec])
        mpz_mul(ans.unit, ans.unit, self.unit)
        _sig_off
        # The relative precision is now the minimum of the relative precisions of self and right (though this may decrease)
        if self.relprec <= (<pAdicCappedRelativeElement>right).relprec:
            relprec = self.relprec
        else:
            relprec = (<pAdicCappedRelativeElement>right).relprec
        # If right's valuation is less than self's, then the normal quotient is integral and we return that.
        if self.ordp >= (<pAdicCappedRelativeElement>right).ordp:
            ans.ordp = self.ordp -(<pAdicCappedRelativeElement>right).ordp
            ans.set_precs(relprec)
            if mpz_cmp(ans.unit, self.prime_pow.dense_list[relprec]) >= 0:
                ans._normalized = 0
            else:
                ans._normalized = 1
        # Otherwise, we have to do floordivision on the unit part of the answer
        else:
            diff = (<pAdicCappedRelativeElement>right).ordp - self.ordp
            # If we're dividing by something where the difference in valuations is bigger than the relative precision, we can only get zero.
            if diff >= relprec:
                ans.set_inexact_zero(0)
            else:
                # Otherwise, our relative precision goes down by the difference in valuations, and we set ans.unit to ans.unit // ppow.
                relprec = relprec - diff
                _sig_on
                mpz_fdiv_q(ans.unit, ans.unit, self.prime_pow.dense_list[diff])
                _sig_off
                ans.ordp = 0
                ans._normalized = 0
                ans.set_precs(relprec)
        return ans

    def __lshift__(self, shift):
        r"""
        Multiplies self by p^shift.

        EXAMPLES:
        sage: a = Zp(5)(17); a
        2 + 3*5 + O(5^20)
        sage: a << 2
        2*5^2 + 3*5^3 + O(5^22)
        """
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        return (<pAdicCappedRelativeElement>self)._lshift_c(mpz_get_si((<Integer>shift).value))

    cdef pAdicCappedRelativeElement _lshift_c(pAdicCappedRelativeElement self, long shift):
        cdef pAdicCappedRelativeElement ans
        if mpz_sgn(self.unit) == -1:
            return self
        if self.prime_pow.in_field == 0 and shift < 0 and -shift > self.ordp:
            ans = self._rshift_c(-shift)
        else:
            ans = self._new_c()
            ans.relprec = self.relprec
            mpz_set(ans.unit, self.unit)
            ans.ordp = self.ordp + shift
            ans._normalized = self._normalized
        return ans

    def __rshift__(pAdicCappedRelativeElement self, shift):
        r"""
        Divides self by p^shift.  If self is a ring element, throws away the non-integral part.

        EXAMPLES:
        sage: a = Zp(5)(17); a
        2 + 3*5 + O(5^20)
        sage: a >> 1
        3 + O(5^19)
        sage: a = Qp(5)(17); a >> 1
        2*5^-1 + 3 + O(5^19)
        """
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        return self._rshift_c(mpz_get_si((<Integer>shift).value))

    cdef pAdicCappedRelativeElement _rshift_c(pAdicCappedRelativeElement self, long shift):
        cdef pAdicCappedRelativeElement ans
        cdef long relprec, diff
        if mpz_sgn(self.unit) == -1:
            return self
        ans = self._new_c()
        if self.prime_pow.in_field == 1 or shift <= self.ordp:
            ans.relprec = self.relprec
            mpz_set(ans.unit, self.unit)
            ans.ordp = self.ordp - shift
            ans._normalized = self._normalized
        else:
            diff = shift - self.ordp
            if diff >= self.relprec:
                ans.set_inexact_zero(0)
            else:
                relprec = self.relprec - diff
                _sig_on
                mpz_fdiv_q(ans.unit, self.unit, self.prime_pow.dense_list[diff])
                _sig_off
                ans.ordp = 0
                ans.set_precs(relprec)
                ans._normalized = 0
        return ans

    def _integer_(pAdicCappedRelativeElement self):
        r"""
        Converts self to an integer

        sage: R = Zp(5, prec = 4); a = R(642); a
        2 + 3*5 + O(5^4)
        sage: Integer(a)
        17
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        if mpz_sgn(self.unit) == -1:
            mpz_set_ui(ans.value, 0)
            return ans
        if self.ordp < 0:
            raise ValueError, "Cannot form an integer out of a p-adic field element with negative valuation"
        if self.ordp <= self.prime_pow._cache_limit:
            if mpz_size(self.prime_pow.modulus) > 10000:
                _sig_on
                mpz_mul(ans.value, self.prime_pow.dense_list[self.ordp], self.unit)
                _sig_off
            else:
                mpz_mul(ans.value, self.prime_pow.dense_list[self.ordp], self.unit)
        else:
            _sig_on
            mpz_pow_ui(ans.value, self.prime_pow.prime.value, self.ordp)
            mpz_mul(ans.value, ans.value, self.unit)
            _sig_off
        return ans

    cdef RingElement _mul_c_impl(self, RingElement right):
        r"""
        sage: R = Zp(5)
        sage: a = R(2385,11); a
        2*5 + 4*5^3 + 3*5^4 + O(5^11)
        sage: b = R(2387625, 16); b
        5^3 + 4*5^5 + 2*5^6 + 5^8 + 5^9 + O(5^16)
        sage: a * b
        2*5^4 + 2*5^6 + 4*5^7 + 2*5^8 + 3*5^10 + 5^11 + 3*5^12 + 4*5^13 + O(5^14)
        """
        cdef pAdicCappedRelativeElement ans
        if mpz_sgn(self.unit) == -1:
            return self
        if mpz_sgn((<pAdicCappedRelativeElement>right).unit) == -1:
            return right
        self._normalize()
        (<pAdicCappedRelativeElement>right)._normalize()
        ans = self._new_c()
        # Do we need to do overflow checking here?
        ans.ordp = self.ordp + (<pAdicCappedRelativeElement>right).ordp
        if self.relprec <= (<pAdicCappedRelativeElement>right).relprec:
            ans.set_precs(self.relprec)
        else:
            ans.set_precs((<pAdicCappedRelativeElement>right).relprec)
        mpz_mul(ans.unit, self.unit, (<pAdicCappedRelativeElement>right).unit)
        if mpz_cmp(ans.unit, self.prime_pow.dense_list[ans.relprec]) >= 0:
            ans._normalized = 0
        else:
            ans._normalized = 1
        return ans

    cdef RingElement _div_c_impl(self, RingElement right):
        cdef pAdicCappedRelativeElement ans
        if mpz_sgn((<pAdicCappedRelativeElement>right).unit) == -1:
            raise ZeroDivisionError, "cannot divide by zero"
        elif mpz_sgn((<pAdicCappedRelativeElement>right).unit) == 0:
            raise PrecisionError, "cannot divide by something indistinguishable from zero"
        ans = self._new_c()
        self._normalize()
        (<pAdicCappedRelativeElement>right)._normalize()
        if ans.prime_pow.in_field == 0:
            ans._parent = self._parent.fraction_field()
            ans.prime_pow = ans._parent.prime_pow
        if mpz_sgn(self.unit) == -1:
            mpz_set(ans.unit, self.unit)
        else:
            ans.ordp = self.ordp - (<pAdicCappedRelativeElement>right).ordp
            if self.relprec <= (<pAdicCappedRelativeElement>right).relprec:
                ans.set_precs(self.relprec)
            else:
                ans.set_precs((<pAdicCappedRelativeElement>right).relprec)
            _sig_on
            mpz_invert(ans.unit, (<pAdicCappedRelativeElement>right).unit, self.prime_pow.dense_list[ans.relprec])
            mpz_mul(ans.unit, ans.unit, self.unit)
            _sig_off
            if mpz_cmp(ans.unit, self.prime_pow.dense_list[ans.relprec]) >= 0:
                ans._normalized = 0
            else:
                ans._normalized = 1
        return ans

    def add_bigoh(self, absprec):
        """
        Returns a new element with absolute precision decreased to absprec
        INPUT:
            self -- a p-adic element
            absprec -- an integer
        OUTPUT:
            element -- self with precision set to the minimum of self's precision and absprec

        EXAMPLE:
            sage: R = Zp(7,4,'capped-rel','series'); a = R(8); a.add_bigoh(1)
            1 + O(7)
            sage: b = R(0); b.add_bigoh(3)
            O(7^3)
            sage: R = Qp(7,4); a = R(8); a.add_bigoh(1)
            1 + O(7)
            sage: b = R(0); b.add_bigoh(3)
            O(7^3)
        """
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        cdef pAdicCappedRelativeElement ans
        cdef long aprec
        if mpz_sgn(self.unit) == -1 or mpz_cmp_ui((<Integer>absprec).value, self.ordp) <= 0:
            ans = self._new_c()
            ans.set_inexact_zero(mpz_get_si((<Integer>absprec).value))
            return ans
        # Do we need to worry about overflow?
        aprec = self.ordp + self.relprec
        if mpz_cmp_ui((<Integer>absprec).value, aprec) >= 0:
            return self
        ans = self._new_c()
        ans.ordp = self.ordp
        ans.set_precs(mpz_get_ui((<Integer>absprec).value) - self.ordp)
        mpz_set(ans.unit, self.unit)
        if mpz_cmp(self.unit, self.prime_pow.dense_list[ans.relprec]) >= 0:
            ans._normalized = 0
        else:
            ans._normalized = self._normalized
        return ans

    def copy(self):
        cdef pAdicCappedRelativeElement ans
        ans = self._new_c()
        ans.relprec = self.relprec
        ans.ordp = self.ordp
        ans._normalized = self._normalized
        mpz_set(ans.unit, self.unit)
        return ans

    def exp_artin_hasse(self):
        raise NotImplementedError
        #E_p(x) = exp(x + x^p/p + x^(p^2)/p^2 + ...)

    def gamma(self):
        raise NotImplementedError

    def _is_exact_zero(self):
        return bool(mpz_sgn(self.unit) == -1)

    def is_zero(self, absprec = None):
        r"""
        Returns whether self is zero modulo $p^{\mbox{absprec}}$.

        INPUT:
            self -- a p-adic element
            absprec -- an integer
        OUTPUT:
            boolean -- whether self is zero

        """
        if absprec is None:
            return bool(mpz_sgn(self.unit) <= 0)
        if mpz_sgn(self.unit) == -1:
            return True
        self._normalize()
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        elif mpz_sgn(self.unit) == 0:
            if mpz_fits_slong_p((<Integer>absprec).value) == 0 or self.ordp < mpz_get_si((<Integer>absprec).value):
                raise PrecisionError, "Not enough precision to determine if element is zero"
            else:
                return True
        return bool(self.ordp >= mpz_get_si((<Integer>absprec).value))

    def is_equal_to(self, right, absprec=None):
        r"""
        Returns whether self is equal to right modulo $p^{\mbox{absprec}}$.

        INPUT:
            self -- a p-adic element
            right -- a p-addic element
            absprec -- an integer
        OUTPUT:
            boolean -- whether self is equal to right

        """
        cdef mpz_t tmp, tmp2
        if not self.parent() is right.parent():
            right = self.parent()(right)
        (<pAdicCappedRelativeElement>self)._normalize()
        (<pAdicCappedRelativeElement>right)._normalize()
        if absprec is None:
            if mpz_sgn(self.unit) <= 0:
                if mpz_sgn((<pAdicCappedRelativeElement>right).unit) <= 0:
                    return True
                else:
                    return False
            else:
                if mpz_sgn((<pAdicCappedRelativeElement>right).unit) <= 0:
                    return False
                elif self.ordp != (<pAdicCappedRelativeElement>right).ordp:
                    return False
                elif self.relprec <= (<pAdicCappedRelativeElement>right).relprec:
                    if mpz_cmp((<pAdicCappedRelativeElement>right).unit, self.prime_pow.dense_list[self.relprec]) >= 0:
                        mpz_init(tmp)
                        mpz_mod(tmp, (<pAdicCappedRelativeElement>right).unit, self.prime_pow.dense_list[self.relprec])
                        if mpz_cmp(tmp, self.unit) == 0:
                            mpz_clear(tmp)
                            return True
                        else:
                            mpz_clear(tmp)
                            return False
                    else:
                        return bool(mpz_cmp(self.unit, (<pAdicCappedRelativeElement>right).unit) == 0)
                else:
                    if mpz_cmp(self.prime_pow.dense_list[(<pAdicCappedRelativeElement>right).relprec], self.unit) <= 0:
                        mpz_init(tmp)
                        mpz_mod(tmp, self.unit, self.prime_pow.dense_list[(<pAdicCappedRelativeElement>right).relprec])
                        if mpz_cmp(tmp, (<pAdicCappedRelativeElement>right).unit) == 0:
                            mpz_clear(tmp)
                            return True
                        else:
                            mpz_clear(tmp)
                            return False
                    else:
                        return bool(mpz_cmp(self.unit, (<pAdicCappedRelativeElement>right).unit) == 0)
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        cdef unsigned long aprec
        aprec = mpz_get_ui((<Integer>absprec).value)
        if mpz_sgn(self.unit) == -1:
            if mpz_sgn((<pAdicCappedRelativeElement>right).unit) == -1:
                return True
            if aprec > (<pAdicCappedRelativeElement>right).ordp + (<pAdicCappedRelativeElement>right).relprec:
                raise PrecisionError, "Elements not known to enough precision"
            elif mpz_sgn((<pAdicCappedRelativeElement>right).unit) == 0:
                return True
            else:
                return False
        if aprec > (<pAdicCappedRelativeElement>self).ordp + (<pAdicCappedRelativeElement>self).relprec:
            raise PrecisionError, "Elements not known to enough precision"
        if mpz_sgn((<pAdicCappedRelativeElement>right).unit) == -1:
            if mpz_sgn((<pAdicCappedRelativeElement>self).unit) == 0:
                return True
            else:
                return False
        if aprec > (<pAdicCappedRelativeElement>right).ordp + (<pAdicCappedRelativeElement>right).relprec:
            raise PrecisionError, "Elements not known to enough precision"
        # We now know that both self and right have enough precision to determine modulo p^absprec
        if self.ordp != (<pAdicCappedRelativeElement>right).ordp:
            return False
        mpz_init(tmp)
        mpz_init(tmp2)
        aprec = aprec - self.ordp
        mpz_mod(tmp, self.unit, self.prime_pow.dense_list[aprec])
        mpz_mod(tmp2, (<pAdicCappedRelativeElement>right).unit, self.prime_pow.dense_list[aprec])
        if mpz_cmp(tmp, tmp2) == 0:
            mpz_clear(tmp)
            mpz_clear(tmp2)
            return True
        else:
            mpz_clear(tmp)
            mpz_clear(tmp2)
            return False

    def lift(self):
        """
        Return an integer or rational congruent to self modulo self's precision.  If a rational is returned, its denominator will eqaul p^ordp(self).

        INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- a integer congruent to self mod $p^{\mbox{prec}}$
        EXAMPLE:
            sage: R = Zp(7,4,'capped-rel'); a = R(8); a.lift()
            8
            sage: R = Qp(7,4); a = R(8); a.lift()
            8
            sage: R = Qp(7,4); a = R(8/7); a.lift()
            8/7
        """
        cdef Integer ans
        cdef Rational ansr
        self._normalize()
        if mpz_sgn(self.unit) <= 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 0)
            return ans
        elif self.ordp < 0:
            ansr = PY_NEW(Rational)
            mpz_set(mpq_numref(ansr.value), self.unit)
            if self.ordp < -self.prime_pow._cache_limit:
                _sig_on
                mpz_pow_ui(mpq_denref(ansr.value), self.prime_pow.prime.value, -self.ordp)
                _sig_off
            else:
                mpz_set(mpq_denref(ansr.value), self.prime_pow.dense_list[-self.ordp])
            return ansr
        elif self.ordp == 0:
            ans = PY_NEW(Integer)
            mpz_set(ans.value, self.unit)
            return ans
        else:
            ans = PY_NEW(Integer)
            if self.ordp > self.prime_pow._cache_limit:
                _sig_on
                mpz_pow_ui(ans.value, self.prime_pow.prime.value, self.ordp)
                mpz_mul(ans.value, ans.value, self.unit)
                _sig_off
            else:
                mpz_mul(ans.value, self.unit, self.prime_pow.dense_list[self.ordp])
            return ans

    def lift_to_precision(self, absprec):
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        return self.lift_to_precision_c(mpz_get_si((<Integer>absprec).value))

    cdef pAdicCappedRelativeElement lift_to_precision_c(pAdicCappedRelativeElement self, long absprec):
        cdef long relprec
        cdef pAdicCappedRelativeElement ans
        if mpz_sgn(self.unit) == -1:
            return self
        elif mpz_sgn(self.unit) == 0:
            if self.ordp + self.relprec >= absprec: # have to add these since self might not be normalized
                return self
            else:
                ans = self._new_c()
                ans.set_inexact_zero(absprec)
                return ans
        relprec = absprec - self.ordp
        if relprec <= self.relprec:
            return self
        else:
            ans = self._new_c()
            mpz_set(ans.unit, self.unit)
            ans.set_precs(relprec)
            ans.ordp = self.ordp
            ans._normalized = self._normalized
            return ans

    def list(self, list_mode = 'simple'):
        """
        Returns a list of coefficients in a power series expansion of self in terms of p.  If self is a field element, they start at p^valuation, if a ring element at p^0.
        INPUT:
            self -- a p-adic element
        OUTPUT:
            list -- the list of coeficients of self
        EXAMPLES:
            sage: R = Zp(7,3,'capped-rel'); a = R(2*7+7**2); a.list()
            [0, 2, 1]
            sage: R = Qp(7,4); a = R(2*7+7**2); a.list()
            [2, 1]

        NOTE:
            use slice operators to get a particular range

        """
        cdef Integer ordp
        cdef pAdicCappedRelativeElement zero
        if mpz_sgn(self.unit) <= 0:
            return []
        if list_mode == 'teichmuller':
            ulist = self.teichmuller_list()
        else:
            ulist = self.base_p_list(self.unit, list_mode)
        if self.prime_pow.in_field == 0 and self.ordp > 0:
            zero = self._new_c()
            zero.set_exact_zero()
            ordp = PY_NEW(Integer)
            mpz_set_si(ordp.value, self.ordp)
            return [zero]*ordp + ulist
        else:
            return ulist

# Working on this
    cdef teichmuller_list(pAdicCappedRelativeElement self):
        # May eventually want to add a dict to store teichmuller lifts already seen, if p small enough
        cdef unsigned long curpower, preccap, i
        cdef mpz_t tmp
        cdef pAdicCappedRelativeElement list_elt
        cdef PowComputer_class powerer
        ans = PyList_New(0)
        preccap = self.relprec
        curpower = preccap
        self._normalize()
        if mpz_sgn(self.unit) <= 0:
            return ans
        mpz_init_set(tmp, self.unit)
        if self.prime_pow.in_field == 0 and self.ordp > 0:
            list_elt = self._new_c()
            list_elt.set_exact_zero()
            for i from 0 <= i < self.ordp:
                PyList_Append(ans, list_elt)
        while mpz_sgn(tmp) != 0:
            curpower -= 1
            list_elt = self._new_c()
            mpz_mod(list_elt.unit, tmp, self.prime_pow.prime.value)
            if mpz_sgn(list_elt.unit) == 0:
                list_elt.set_exact_zero()
                mpz_divexact(tmp, tmp, self.prime_pow.prime.value)
            else:
                list_elt.ordp = 0
                list_elt.set_precs(preccap)
                sage.rings.padics.padic_generic_element.teichmuller_set_c(list_elt.unit, self.prime_pow.prime.value, self.prime_pow.dense_list[preccap])
                list_elt._normalized = 1
                mpz_sub(tmp, tmp, list_elt.unit)
                mpz_divexact(tmp, tmp, self.prime_pow.prime.value)
                mpz_mod(tmp, tmp, self.prime_pow.dense_list[curpower])
            PyList_Append(ans, list_elt)
        mpz_clear(tmp)
        return ans

    def log_artin_hasse(self):
        raise NotImplementedError

    def padded_list(self, n, list_mode = 'simple'):
        """
        Returns a list of coeficiants of p starting with $p^0$ up to $p^n$ exclusive (padded with zeros if needed).  If a field element, starts at p^val instead.
        INPUT:
            self -- a p-adic element
            n - an integer
        OUTPUT:
            list -- the list of coeficients of self
        EXAMPLES:
            sage: R = Zp(7,3,'capped-rel'); a = R(2*7+7**2); a.padded_list(5)
            [0, 2, 1, 0, 0]
            sage: R = Qp(7,3); a = R(2*7+7**2); a.padded_list(5)
            [2, 1, 0, 0]
            sage: a.padded_list(3)
            [2, 1]

        NOTE:
            this differs from the padded_list method of padic_field_element
            the slice operators throw an error if asked for a slice above the precision
        """
        if list_mode == 'simple' or list_mode == 'smallest':
            zero = Integer(0)
        else:
            zero = self.parent()(0)
        self._normalize()
        L = self.list()
        if self.prime_pow.in_field == 1:
            if mpz_sgn(self.unit) > 0:
                n -= self.valuation()
            else:
                n = 0
        return L[:n] + [zero] * (n - len(L))

    def precision_absolute(pAdicCappedRelativeElement self):
        """
        Returns the absolute precision of self
         INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the absolute precision of self
        EXAMPLES:
            sage: R = Zp(7,3,'capped-rel'); a = R(7); a.precision_absolute()
            4
            sage: R = Qp(7,3); a = R(7); a.precision_absolute()
            4
        """
        cdef Integer ans
        if mpz_sgn(self.unit) == -1:
            return infinity
        ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.ordp + self.relprec)
        return ans

    def precision_relative(self):
        """
        Returns the relative precision of self
         INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the relative precision of self
        EXAMPLES:
            sage: R = Zp(7,3,'capped-rel'); a = R(7); a.precision_relative()
            3
            sage: R = Qp(7,3); a = R(7); a.precision_relative()
            3
       """
        self._normalize()
        cdef Integer ans
        ans = PY_NEW(Integer)
        if mpz_sgn(self.unit) < 0:
            mpz_set_ui(ans.value, 0)
        else:
            mpz_set_si(ans.value, self.relprec)
        return ans

    def residue(self, absprec):
        """
        Reduces this mod $p^absprec$
        INPUT:
            self -- a p-adic element
            absprec - an integer
        OUTPUT:
            element of Z/(p^absprec Z) -- self reduced mod p^absprec
        EXAMPLES:
            sage: R = Zp(7,4,'capped-rel'); a = R(8); a.residue(1)
            1
            sage: R = Qp(7,4,'capped-rel'); a = R(8); a.residue(1)
            1
        """
        cdef Integer selfvalue, modulus
        cdef PowComputer_class powerer
        cdef unsigned long aprec
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if absprec > self.precision_absolute():
            raise PrecisionError, "Not enough precision known in order to compute residue."
        elif absprec < 0:
            raise ValueError, "cannot reduce modulo a negative power of p"
        aprec = mpz_get_ui((<Integer>absprec).value)
        if self.ordp < 0:
            self._normalize()
        if self.ordp < 0:
            raise ValueError, "Element must have non-negative valuation in order to compute residue."
        if mpz_sgn(self.unit) <= 0:
            modulus = PY_NEW(Integer)
            mpz_set(modulus.value, self.prime_pow.dense_list[aprec])
            selfvalue = PY_NEW(Integer)
            mpz_set_ui(selfvalue.value, 0)
            return Mod(selfvalue, modulus)
        else:
            # Need to do this better.
            selfvalue = PY_NEW(Integer)
            modulus = PY_NEW(Integer)
            if self.ordp < self.prime_pow._cache_limit:
                mpz_mul(selfvalue.value, self.prime_pow.dense_list[self.ordp], self.unit)
            else:
                _sig_on
                mpz_pow_ui(selfvalue.value, self.prime_pow.prime.value, self.ordp)
                mpz_mul(selfvalue.value, selfvalue.value, self.unit)
                _sig_off
            mpz_set(modulus.value, self.prime_pow.dense_list[aprec])
            return Mod(selfvalue, modulus)

    def unit_part(self):
        r"""
        Returns the unit part of self.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            p-adic element -- the unit part of self
        EXAMPLES:
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
        """
        return self.unit_part_c()

    cdef pAdicCappedRelativeElement unit_part_c(self):
        self._normalize()
        if mpz_sgn(self.unit) < 0:
            raise ValueError, "unit part of 0 not defined"
        cdef pAdicCappedRelativeElement ans
        ans = self._new_c()
        mpz_set(ans.unit, self.unit)
        ans.set_precs(self.relprec)
        ans.ordp = 0
        ans._normalized = 1
        return ans

    def valuation(self):
        """
        Returns the valuation of self.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the valuation of self

        EXAMPLES:
            sage: R = Zp(17, 4,'capped-rel')
            sage: a = R(2*17^2)
            sage: a.valuation()
            2
            sage: R = Zp(5, 4,'capped-rel')
            sage: R(0).valuation()
            +Infinity

        TESTS:
            sage: R(1).valuation()
            0
            sage: R(2).valuation()
            0
            sage: R(5).valuation()
            1
            sage: R(10).valuation()
            1
            sage: R(25).valuation()
            2
            sage: R(50).valuation()
            2
            sage: R = Qp(17, 4)
            sage: a = R(2*17^2)
            sage: a.valuation()
            2
            sage: R = Qp(5, 4)
            sage: R(0).valuation()
            +Infinity
            sage: R(1).valuation()
            0
            sage: R(2).valuation()
            0
            sage: R(5).valuation()
            1
            sage: R(10).valuation()
            1
            sage: R(25).valuation()
            2
            sage: R(50).valuation()
            2
            sage: R(1/2).valuation()
            0
            sage: R(1/5).valuation()
            -1
            sage: R(1/10).valuation()
            -1
            sage: R(1/25).valuation()
            -2
            sage: R(1/50).valuation()
            -2
        """
        return self.valuation_c()

    cdef valuation_c(self):
        cdef Integer ans
        if mpz_sgn(self.unit) == -1:
            return infinity
        else:
            self._normalize()
            ans = PY_NEW(Integer)
            mpz_set_si(ans.value, self.ordp)
            return ans

    cdef val_unit_c(self):
        cdef Integer val
        cdef pAdicCappedRelativeElement unit
        if mpz_sgn(self.unit) == -1:
            unit = self._new_c()
            unit.set_inexact_zero(0)
            return (infinity, unit)
        self._normalize()
        unit = self._new_c()
        mpz_set(unit.unit, self.unit)
        unit.set_precs(self.relprec)
        unit.ordp = 0
        unit._normalized = 1
        val = PY_NEW(Integer)
        mpz_set_si(val.value, self.ordp)
        return (val, unit)

    def val_unit(self):
        """
        Returns an pair (self.valuation(), self.unit_part()).

        If self is zero (either inexact or exact), the second return value is 0 + O(p^0).

        EXAMPLES:
        sage: R = Zp(5); a = R(75, 20); a
        3*5^2 + O(5^20)
        sage: a.val_unit()
        (2, 3 + O(5^18))
        sage: R(0).val_unit()
        (+Infinity, O(5^0))
        sage: R(0, 10).val_unit()
        (10, O(5^0))
        """
        return self.val_unit_c()

    def __hash__(self):
        return self._hash()

    cdef long _hash(self) except -1:
        cdef Integer ans
        self._normalize()
        if mpz_sgn(self.unit) <= 0:
            return 0
        else:
            ans = PY_NEW(Integer)
            if self.ordp <= self.prime_pow._cache_limit:
                mpz_mul(ans.value, self.unit, self.prime_pow.dense_list[self.ordp])
            else:
                _sig_on
                mpz_pow_ui(ans.value, self.prime_pow.prime.value, self.ordp)
                mpz_mul(ans.value, ans.value, self.unit)
                _sig_off
            return hash(ans)

    def _teichmuller_set(self, Integer n, Integer absprec):
        cdef unsigned long aprec
        if mpz_divisible_p(n.value, self.prime_pow.prime.value):
            self.set_exact_zero()
        else:
            mpz_set(self.unit, n.value)
            self.ordp = 0
            if mpz_fits_ulong_p(absprec.value) == 0:
                aprec = self.prime_pow._cache_limit
            if mpz_sgn(absprec.value) != 1:
                raise ValueError, "can only compute to positive precision"
            aprec = mpz_get_ui(absprec.value)
            if aprec > self.prime_pow._cache_limit:
                aprec = self.prime_pow._cache_limit
            self.set_precs(aprec)
            sage.rings.padics.padic_generic_element.teichmuller_set_c(self.unit, self.prime_pow.prime.value, self.prime_pow.dense_list[aprec])

def unpickle_pcre_v1(R, unit, ordp, relprec):
    cdef pAdicCappedRelativeElement ans
    ans = PY_NEW(pAdicCappedRelativeElement)
    ans._parent = R
    ans.prime_pow = <PowComputer_class>R.prime_pow
    mpz_init_set(ans.unit, (<Integer>unit).value)
    ans.ordp = mpz_get_si((<Integer>ordp).value)
    ans.relprec = mpz_get_ui((<Integer>relprec).value)
    ans._normalized = 0
    return ans
