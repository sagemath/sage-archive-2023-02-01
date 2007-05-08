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
        #print "beginning"
        cdef RingElement ordp
        mpz_init(self.unit)
        mpz_init(self.ordp)
        mpz_init_set(self.p, (<Integer>parent.prime()).value)
        mpz_init(self.modulus)
        if parent.is_field():
            self.in_field = 1
        else:
            self.in_field = 0
        CommutativeRingElement.__init__(self, parent)
        if empty:
            return
        cdef PowComputer_class powerer
        powerer = <PowComputer_class> parent.prime_pow
        if not absprec is infinity and not relprec is infinity:
            raise ValueError, "can only specify one of absprec and relprec"
        if absprec is infinity:
            if relprec > parent.precision_cap():
                relprec = parent.precision_cap()

        if isinstance(x, pAdicGenericElement):
            if self.in_field == 0 and x.valuation() < 0:
                raise ValueError, "element has negative valuation."
            if parent.prime() != x.parent().prime():
                raise TypeError, "Cannot coerce between p-adic parents with different primes"
        if isinstance(x, pAdicLazyElement):
            if absprec is infinity:
                relprec = min(relprec, parent.precision_cap())
                try:
                    x.set_precision_relative(relprec)
                except PrecisionError:
                    pass
            else:
                try:
                    x.set_precision_absolute(absprec)
                except PrecisionError:
                    pass
        if isinstance(x, pAdicBaseGenericElement):
            ordp = x.valuation()
            if ordp >= absprec:
                if absprec is infinity:
                    self.set_exact_zero()
                else:
                    if not PY_TYPE_CHECK(absprec, Integer):
                        absprec = Integer(absprec)
                    self.set_inexact_zero((<Integer>absprec).value)
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
            k = mpz_remove(modulus, modulus, self.p)
            if mpz_cmp_ui(modulus, 1) == 0:
                if absprec is infinity and relprec == parent.precision_cap(): #this allows you to lift integer_mod elements to higher precision than they are defined.  Subtle bug maybe: user tries to specify relprec = parent.precision_cap?
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
            if ordp is infinity:
                if absprec is infinity:
                    self.set_exact_zero()
                else:
                    self.set_inexact_zero((<Integer>absprec).value)
                return
            #print "now 1"
            relprec = min(relprec, absprec - ordp, parent.precision_cap())
            #print "now 2"
            self.set_from_Integers(ordp, unit, relprec)
            #print "now 3"
        elif isinstance(x, Rational):
            ordp, unit = x.val_unit(parent.prime())
            if ordp is infinity:
                if absprec is infinity:
                    self.set_exact_zero()
                else:
                    self.set_inexact_zero((<Integer>absprec).value)
                return
            if self.in_field == 0 and ordp < 0:
                raise ValueError, "p divides the denominator"
            # todo: make the following faster
            relprec = min(relprec, absprec - ordp, parent.precision_cap())
            unit = Mod(unit, parent.prime_pow(relprec)).lift()
            self.set_from_Integers(ordp, unit, relprec)
        else:
            raise TypeError, "cannot create a p-adic out of %s"%(type(x))
        if mpz_sgn(self.ordp) == -1 and self.in_field == 0:
            raise ValueError, "element has negative valuation."

    cdef void set_exact_zero(pAdicCappedRelativeElement self):
        mpz_set_si(self.unit, -1)

    cdef void set_inexact_zero(pAdicCappedRelativeElement self, mpz_t absprec):
        self.relprec = 0
        mpz_set_ui(self.unit, 0)
        mpz_set(self.ordp, absprec)
        mpz_set_ui(self.modulus, 1)

    cdef void set_precs(pAdicCappedRelativeElement self, unsigned int relprec):
        cdef PowComputer_class powerer
        # Need to check that relprec fits in an int
        self.relprec = relprec
        #print "inner"
        powerer = self._parent.prime_pow
        powerer.pow_mpz_ui(self.modulus, self.relprec)

    cdef void set_from_Integers(pAdicCappedRelativeElement self, Integer ordp, Integer unit, Integer relprec):
        #print "here 1"
        #print type(relprec)
        self.set_precs(mpz_get_ui(relprec.value))
        #print "here 2"
        mpz_set(self.ordp, ordp.value)
        #print "here 3"
        mpz_set(self.unit, unit.value)
        #print "here 4"
        if mpz_sgn(self.unit) == -1 or mpz_cmp(self.unit, self.modulus) >= 0:
            #print "here 5"
            mpz_mod(self.unit, self.unit, self.modulus)

    cdef pAdicCappedRelativeElement _new_c(pAdicCappedRelativeElement self):
        cdef pAdicCappedRelativeElement ans
        ans = PY_NEW(pAdicCappedRelativeElement)
        ans._parent = self._parent
        ans.in_field = self.in_field
        mpz_init_set(ans.p, self.p)
        mpz_init(ans.unit)
        mpz_init(ans.modulus)
        mpz_init(ans.ordp)
        return ans

    def __dealloc__(pAdicCappedRelativeElement self):
        mpz_clear(self.p)
        mpz_clear(self.unit)
        mpz_clear(self.modulus)
        mpz_clear(self.ordp)

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    def __reduce__(self):
        # Necessary for pickling.  See integer.pyx for more info.
        cdef Integer relprec, unit, ordp
        relprec = PY_NEW(Integer)
        mpz_set_ui(relprec.value, self.relprec)
        unit = PY_NEW(Integer)
        mpz_set(unit.value, self.unit)
        ordp = PY_NEW(Integer)
        mpz_set(ordp.value, self.ordp)
        return unpickle_pcre_v1, (self.parent(), unit, ordp, relprec)

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
        if mpz_sgn(self.unit) <= 0: # an exact or inexact zero
            return self
        ans = self._new_c()
        ans.relprec = self.relprec
        mpz_set(ans.modulus, self.modulus)
        mpz_set(ans.ordp, self.ordp)
        mpz_sub(ans.unit, ans.modulus, self.unit)
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
        cdef mpz_t absprec
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
            return ans
        # We next check if self is an exact zero
        elif mpz_sgn((<pAdicCappedRelativeElement>self).unit) == -1:
            return self
        # We next check if self is an inexact zero
        elif mpz_sgn((<pAdicCappedRelativeElement>self).unit) == 0:
            ans = self._new_c()
            mpz_init(absprec)
            mpz_mul(absprec, (<pAdicCappedRelativeElement>self).ordp, (<Integer>right).value)
            ans.set_inexact_zero(absprec)
            mpz_clear(absprec)
            return ans
        else:
            ans = self._new_c()
            # You actually get a bit more relative precision when the exponent is divisible by p...
            ans.set_precs(self.relprec)
            mpz_mul(ans.ordp, (<Integer>right).value, (<pAdicCappedRelativeElement>self).ordp)
            mpz_powm(ans.unit, (<pAdicCappedRelativeElement>self).unit, (<Integer>right).value, ans.modulus)
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
        cdef PowComputer_class powerer
        cdef mpz_t ppow, tmp, modified
        cdef int val_cmp, diff, tmpi
        # If either self or right is an exact zero, we do the appropriate thing
        if mpz_sgn(self.unit) == -1:
            return right
        if mpz_sgn((<pAdicCappedRelativeElement>right).unit) == -1:
            return self
        val_cmp = mpz_cmp(self.ordp, (<pAdicCappedRelativeElement>right).ordp)
        if val_cmp == 0:
            ans = self._new_c()
            # The relative precision of the sum is the minimum of the relative precisions in this case, possibly decreasing if we got cancellation
            if self.relprec < (<pAdicCappedRelativeElement>right).relprec:
                ans.relprec = self.relprec
            else:
                ans.relprec = (<pAdicCappedRelativeElement>right).relprec
            # Since the valuations are the same, we can just add the units
            mpz_add(ans.unit, self.unit, (<pAdicCappedRelativeElement>right).unit)
            if mpz_sgn(ans.unit) == 0:
                mpz_add_ui(ans.ordp, self.ordp, ans.relprec)
                ans.set_precs(0)
            else:
                if mpz_divisible_p(ans.unit, ans.p) != 0: # the valuation of the sum increases
                    diff = mpz_remove(ans.unit, ans.unit, ans.p)
                    mpz_add_ui(ans.ordp, self.ordp, diff)
                    ans.set_precs(ans.relprec - diff)
                else:
                    mpz_set(ans.ordp, self.ordp)
                    ans.set_precs(ans.relprec)
                if mpz_cmp(ans.unit, ans.modulus) >= 0:
                    mpz_mod(ans.unit, ans.unit, ans.modulus)
        elif val_cmp < 0: # self has lower valuation
            ans = self._new_c()
            mpz_init(tmp)
            mpz_sub(tmp, (<pAdicCappedRelativeElement>right).ordp, self.ordp)
            # if the difference in valuations is bigger than self's relprec, return self
            if mpz_cmp_si(tmp, self.relprec) >= 0:
                mpz_clear(tmp)
                return self
            # The valuation of the answer is self.ordp
            mpz_set(ans.ordp, self.ordp)
            # The relative precision of the answer is the minimum of self.relprec and tmp + right.relprec
            tmpi = mpz_get_si(tmp)
            if self.relprec <= tmpi + (<pAdicCappedRelativeElement>right).relprec:
                ans.set_precs(self.relprec)
            else:
                ans.set_precs(tmpi + (<pAdicCappedRelativeElement>right).relprec)
            # Now we compute p^tmp
            mpz_init(ppow)
            powerer = <PowComputer_class> self.parent().prime_pow
            powerer.pow_mpz_ui(ppow, tmpi)
            # ans.unit = self.unit + ppow * right.unit
            mpz_mul(tmp, ppow, (<pAdicCappedRelativeElement>right).unit)
            mpz_add(ans.unit, self.unit, tmp)
            if mpz_cmp(ans.unit, ans.modulus) >= 0:
                mpz_mod(ans.unit, ans.unit, ans.modulus)
            mpz_clear(ppow)
            mpz_clear(tmp)
        else: # right has lower valuation
            ans = self._new_c()
            mpz_init(tmp)
            mpz_sub(tmp, self.ordp, (<pAdicCappedRelativeElement>right).ordp)
            # if the difference in valuations is bigger than right's relprec, return right
            if mpz_cmp_si(tmp, (<pAdicCappedRelativeElement>right).relprec) >= 0:
                mpz_clear(tmp)
                return right
            # The valuation of the answer is right.ordp
            mpz_set(ans.ordp, (<pAdicCappedRelativeElement>right).ordp)
            # The relative precision of the answer is the minimum of right.relprec and tmp + self.relprec
            tmpi = mpz_get_si(tmp)
            if (<pAdicCappedRelativeElement>right).relprec <= self.relprec + tmpi:
                ans.set_precs((<pAdicCappedRelativeElement>right).relprec)
            else:
                ans.set_precs(tmpi + self.relprec)
            # Now we compute p^tmp
            mpz_init(ppow)
            powerer = <PowComputer_class> self.parent().prime_pow
            powerer.pow_mpz_mpz(ppow, tmp)
            # ans.unit = right.unit + ppow * self.unit
            mpz_mul(tmp, ppow, self.unit)
            mpz_add(ans.unit, (<pAdicCappedRelativeElement>right).unit, tmp)
            if mpz_cmp(ans.unit, ans.modulus) >= 0:
                mpz_mod(ans.unit, ans.unit, ans.modulus)
            mpz_clear(ppow)
            mpz_clear(tmp)
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
        if mpz_sgn(self.unit) == -1:
            raise ZeroDivisionError, "cannot divide by zero"
        if self.relprec == 0:
            raise PrecisionError, "cannot divide by something indistinguishable from zero."
        cdef pAdicCappedRelativeElement ans
        ans = self._new_c()
        if ans.in_field == 0:
            ans._parent = self.parent().fraction_field()
            ans.in_field = 1
        mpz_neg(ans.ordp, (<pAdicCappedRelativeElement>self).ordp)
        ans.set_precs((<pAdicCappedRelativeElement>self).relprec)
        mpz_invert(ans.unit, (<pAdicCappedRelativeElement>self).unit, ans.modulus)
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
        cdef int relprec, diff
        cdef mpz_t tmp
        cdef PowComputer_class powerer
        cdef mpz_t absprec
        # For fields, we define floor division as normal division and % as always 0
        if self.in_field == 1:
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
        mpz_invert(ans.unit, (<pAdicCappedRelativeElement>right).unit, (<pAdicCappedRelativeElement>right).modulus)
        mpz_mul(ans.unit, ans.unit, self.unit)
        # The relative precision is now the minimum of the relative precisions of self and right (though this may decrease)
        if self.relprec < (<pAdicCappedRelativeElement>right).relprec:
            relprec = self.relprec
        else:
            relprec = (<pAdicCappedRelativeElement>right).relprec
        # If right's valuation is less than self's, then the normal quotient is integral and we return that.
        if mpz_cmp(self.ordp, (<pAdicCappedRelativeElement>right).ordp) >= 0:
            mpz_sub(ans.ordp, self.ordp, (<pAdicCappedRelativeElement>right).ordp)
            ans.set_precs(relprec)
            mpz_mod(ans.unit, ans.unit, ans.modulus)
        # Otherwise, we have to do floordivision on the unit part of the answer
        else:
            mpz_init(tmp)
            mpz_sub(tmp, (<pAdicCappedRelativeElement>right).ordp, self.ordp)
            # If we're dividing by something where the difference in valuations is bigger than the relative precision, we can only get zero.
            if mpz_cmp_ui(tmp, relprec) >= 0:
                mpz_init(absprec)
                mpz_set_ui(absprec, 0)
                ans.set_inexact_zero(absprec)
                mpz_clear(absprec)
            else:
                # Otherwise, our relative precision goes down by the difference in valuations, and we set ans.unit to ans.unit // ppow.
                relprec = relprec - mpz_get_si(tmp)
                powerer = <PowComputer_class> self.parent().prime_pow
                powerer.pow_mpz_mpz(tmp, tmp)
                mpz_fdiv_q(ans.unit, ans.unit, tmp)
                # If this floor division has given us zero, we set ans to be an inexact zero.
                if mpz_sgn(ans.unit) == 0:
                    mpz_init(absprec)
                    mpz_set_ui(absprec, relprec)
                    ans.set_inexact_zero(absprec)
                    mpz_clear(absprec)
                else:
                    # Otherwise, we pull off powers of p into ans.ordp.
                    if mpz_divisible_p(ans.unit, ans.p) != 0:
                        diff = mpz_remove(ans.unit, ans.unit, ans.p)
                        relprec = relprec - diff
                        mpz_set_ui(ans.ordp, diff)
                    else:
                        mpz_set_ui(ans.ordp, 0)
                    ans.set_precs(relprec)
                    mpz_mod(ans.unit, ans.unit, ans.modulus)
            mpz_clear(tmp)
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
        return (<pAdicCappedRelativeElement>self)._lshift_c((<Integer>shift).value)

    cdef pAdicCappedRelativeElement _lshift_c(pAdicCappedRelativeElement self, mpz_t shift):
        cdef pAdicCappedRelativeElement ans
        cdef mpz_t tmp
        if mpz_sgn(self.unit) == -1:
            return self
        if self.in_field == 0 and mpz_sgn(shift) == -1 and mpz_cmpabs(shift, self.ordp) > 0:
            mpz_init(tmp)
            mpz_neg(tmp, shift)
            ans = self._rshift_c(tmp)
            mpz_clear(tmp)
        else:
            ans = self._new_c()
            ans.relprec = self.relprec
            mpz_set(ans.modulus, self.modulus)
            mpz_set(ans.unit, self.unit)
            mpz_add(ans.ordp, self.ordp, shift)
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
        return self._rshift_c((<Integer>shift).value)

    cdef pAdicCappedRelativeElement _rshift_c(pAdicCappedRelativeElement self, mpz_t shift):
        cdef pAdicCappedRelativeElement ans
        cdef mpz_t tmp, absprec
        cdef int relprec
        cdef PowComputer_class powerer
        if mpz_sgn(self.unit) == -1:
            return self
        ans = self._new_c()
        if self.in_field == 1 or mpz_cmp(shift, self.ordp) <= 0:
            ans.relprec = self.relprec
            mpz_set(ans.modulus, self.modulus)
            mpz_set(ans.unit, self.unit)
            mpz_sub(ans.ordp, self.ordp, shift)
        else:
            mpz_init(tmp)
            mpz_sub(tmp, shift, self.ordp)
            if mpz_cmp_ui(tmp, self.relprec) >= 0:
                mpz_init(absprec)
                mpz_set_ui(absprec, 0)
                ans.set_inexact_zero(absprec)
                mpz_clear(absprec)
            else:
                relprec = self.relprec - mpz_get_si(tmp)
                powerer = <PowComputer_class> self.parent().prime_pow
                powerer.pow_mpz_mpz(tmp, tmp)
                mpz_fdiv_q(ans.unit, self.unit, tmp)
                if mpz_sgn(ans.unit) == 0:
                    mpz_init(absprec)
                    mpz_set_ui(absprec, relprec)
                    ans.set_inexact_zero(absprec)
                    mpz_clear(absprec)
                else:
                    if mpz_divisible_p(ans.unit, ans.p) != 0:
                        diff = mpz_remove(ans.unit, ans.unit, ans.p)
                        relprec = relprec - diff
                        mpz_set_ui(ans.ordp, diff)
                    else:
                        mpz_set_ui(ans.ordp, 0)
                    ans.set_precs(relprec)
            mpz_clear(tmp)
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
        cdef PowComputer_class powerer
        ans = PY_NEW(Integer)
        if mpz_sgn(self.unit) == -1:
            mpz_set_ui(ans.value, 0)
            return ans
        if mpz_sgn(self.ordp) == -1:
            raise ValueError, "Cannot form an integer out of a p-adic field element with negative valuation"
        powerer = self.parent().prime_pow
        powerer.pow_mpz_mpz(ans.value, self.ordp)
        _sig_on
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
        ans = self._new_c()
        mpz_add(ans.ordp, self.ordp, (<pAdicCappedRelativeElement>right).ordp)
        if self.relprec <= (<pAdicCappedRelativeElement>right).relprec:
            ans.set_precs(self.relprec)
        else:
            ans.set_precs((<pAdicCappedRelativeElement>right).relprec)
        mpz_mul(ans.unit, self.unit, (<pAdicCappedRelativeElement>right).unit)
        if mpz_cmp(ans.unit, ans.modulus) >= 0:
            mpz_mod(ans.unit, ans.unit, ans.modulus)
        return ans

    cdef RingElement _div_c_impl(self, RingElement right):
        cdef pAdicCappedRelativeElement ans
        if mpz_sgn((<pAdicCappedRelativeElement>right).unit) == -1:
            raise ZeroDivisionError, "cannot divide by zero"
        elif mpz_sgn((<pAdicCappedRelativeElement>right).unit) == 0:
            raise PrecisionError, "cannot divide by something indistinguishable from zero"
        ans = self._new_c()
        if ans.in_field == 0:
            ans._parent = self._parent.fraction_field()
            ans.in_field = 1
        if mpz_sgn(self.unit) == -1:
            mpz_set(ans.unit, self.unit)
        else:
            mpz_sub(ans.ordp, self.ordp, (<pAdicCappedRelativeElement>right).ordp)
            if self.relprec <= (<pAdicCappedRelativeElement>right).relprec:
                ans.set_precs(self.relprec)
            else:
                ans.set_precs((<pAdicCappedRelativeElement>right).relprec)
            mpz_invert(ans.unit, (<pAdicCappedRelativeElement>right).unit, ans.modulus)
            mpz_mul(ans.unit, ans.unit, self.unit)
            if mpz_cmp(ans.unit, ans.modulus) >= 0:
                mpz_mod(ans.unit, ans.unit, ans.modulus)
        return ans

    def add_bigoh(self, prec):
        """
        Returns a new element with absolute precision decreased to prec
        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            element -- self with precision set to the minimum of  self's precision and prec

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
        if not PY_TYPE_CHECK(prec, Integer):
            prec = Integer(prec)
        cdef pAdicCappedRelativeElement ans
        cdef mpz_t aprec
        if mpz_sgn(self.unit) == -1 or mpz_cmp(self.ordp, (<Integer>prec).value) >= 0:
            ans = self._new_c()
            ans.set_inexact_zero((<Integer>prec).value)
            return ans
        mpz_init(aprec)
        mpz_add_ui(aprec, self.ordp, self.relprec)
        if mpz_cmp(aprec, (<Integer>prec).value) <= 0:
            mpz_clear(aprec)
            return self
        ans = self._new_c()
        mpz_set(ans.ordp, self.ordp)
        mpz_sub(aprec, (<Integer>prec).value, self.ordp)
        ans.set_precs(mpz_get_ui(aprec))
        if mpz_cmp(self.unit, ans.modulus) >= 0:
            mpz_mod(ans.unit, self.unit, ans.modulus)
        else:
            mpz_set(ans.unit, self.unit)
        return ans

    def copy(self):
        cdef pAdicCappedRelativeElement ans
        ans = self._new_c()
        ans.relprec = self.relprec
        mpz_set(ans.modulus, self.modulus)
        mpz_set(ans.ordp, self.ordp)
        mpz_set(ans.unit, self.unit)
        return ans

    def exp_artin_hasse(self):
        raise NotImplementedError
        #E_p(x) = exp(x + x^p/p + x^(p^2)/p^2 + ...)

    def gamma(self):
        raise NotImplementedError

    def _is_exact_zero(self):
        return bool(mpz_sgn(self.unit) == -1)

    def is_zero(self, prec = None):
        r"""
        Returns whether self is zero modulo $p^{\mbox{prec}}$.

        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            boolean -- whether self is zero

        """
        if prec is None:
            return bool(mpz_sgn(self.unit) <= 0)
        if mpz_sgn(self.unit) == -1:
            return True
        if not PY_TYPE_CHECK(prec, Integer):
            prec = Integer(prec)
        elif mpz_sgn(self.unit) == 0:
            if mpz_cmp(self.ordp, (<Integer>prec).value) >= 0:
                return True
            else:
                raise PrecisionError, "Not enough precision to determine if element is zero"
        return bool(mpz_cmp(self.ordp, (<Integer>prec).value) >= 0)

    def is_equal_to(self, right, prec=None):
        r"""
        Returns whether self is equal to right modulo $p^{\mbox{prec}}$.

        INPUT:
            self -- a p-adic element
            right -- a p-addic element
            prec -- an integer
        OUTPUT:
            boolean -- whether self is equal to right

        """
        cdef mpz_t tmp, tmp2, ppow
        cdef PowComputer_class powerer
        if not self.parent() is right.parent():
            right = self.parent()(right)
        if prec is None:
            if mpz_sgn(self.unit) <= 0:
                if mpz_sgn((<pAdicCappedRelativeElement>right).unit) <= 0:
                    return True
                else:
                    return False
            else:
                if mpz_sgn((<pAdicCappedRelativeElement>right).unit) <= 0:
                    return False
                elif mpz_cmp(self.ordp, (<pAdicCappedRelativeElement>right).ordp) != 0:
                    return False
                elif self.relprec <= (<pAdicCappedRelativeElement>right).relprec:
                    if mpz_cmp((<pAdicCappedRelativeElement>right).unit, self.modulus) >= 0:
                        mpz_init(tmp)
                        mpz_mod(tmp, (<pAdicCappedRelativeElement>right).unit, self.modulus)
                        if mpz_cmp(tmp, self.unit) == 0:
                            mpz_clear(tmp)
                            return True
                        else:
                            mpz_clear(tmp)
                            return False
                    else:
                        return bool(mpz_cmp(self.unit, (<pAdicCappedRelativeElement>right).unit) == 0)
                else:
                    if mpz_cmp((<pAdicCappedRelativeElement>right).modulus, self.unit) <= 0:
                        mpz_init(tmp)
                        mpz_mod(tmp, self.unit, (<pAdicCappedRelativeElement>right).modulus)
                        if mpz_cmp(tmp, (<pAdicCappedRelativeElement>right).unit) == 0:
                            mpz_clear(tmp)
                            return True
                        else:
                            mpz_clear(tmp)
                            return False
                    else:
                        return bool(mpz_cmp(self.unit, (<pAdicCappedRelativeElement>right).unit) == 0)
        if not PY_TYPE_CHECK(prec, Integer):
            prec = Integer(prec)
        if mpz_sgn(self.unit) == -1:
            if mpz_sgn((<pAdicCappedRelativeElement>right).unit) == -1:
                return True
            if mpz_cmp((<Integer>right.precision_absolute()).value, (<Integer>prec).value) < 0:
                raise PrecisionError, "Elements not known to enough precision"
            elif mpz_sgn((<pAdicCappedRelativeElement>right).unit) == 0:
                return True
            else:
                return False
        if mpz_cmp((<Integer>self.precision_absolute()).value, (<Integer>prec).value) < 0:
            raise PrecisionError, "Elements not known to enough precision"
        if mpz_sgn((<pAdicCappedRelativeElement>right).unit) == -1:
            if mpz_sgn((<pAdicCappedRelativeElement>self).unit) == 0:
                return True
            else:
                return False
        if mpz_cmp((<Integer>right.precision_absolute()).value, (<Integer>prec).value) < 0:
            raise PrecisionError, "Elements not known to enough precision"
        # We now know that both self and right have enough precision to determine modulo p^prec
        if mpz_cmp((<pAdicCappedRelativeElement>right).ordp, self.ordp) != 0:
            return False
        mpz_init(ppow)
        powerer = <PowComputer_class> self.parent().prime_pow
        mpz_init(tmp)
        mpz_init(tmp2)
        mpz_sub(tmp, (<Integer>prec).value, self.ordp)
        powerer.pow_mpz_mpz(ppow, tmp)
        mpz_mod(tmp, self.unit, ppow)
        mpz_mod(tmp2, (<pAdicCappedRelativeElement>right).unit, ppow)
        if mpz_cmp(tmp, tmp2) == 0:
            return True
        else:
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
        cdef PowComputer_class powerer
        cdef mpz_t tmp
        if mpz_sgn(self.unit) <= 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 0)
            return ans
        elif mpz_sgn(self.ordp) == -1:
            mpz_init(tmp)
            mpz_neg(tmp, self.ordp)
            ansr = PY_NEW(Rational)
            powerer = self.parent().prime_pow
            powerer.pow_mpz_mpz(mpq_denref(ansr.value), tmp)
            mpz_set(mpq_numref(ansr.value), self.unit)
            return ansr
        elif mpz_sgn(self.ordp) == 0:
            ans = PY_NEW(Integer)
            mpz_set(ans.value, self.unit)
            return ans
        else:
            ans = PY_NEW(Integer)
            powerer = self.parent().prime_pow
            powerer.pow_mpz_mpz(ans.value, self.ordp)
            mpz_mul(ans.value, ans.value, self.unit)
            return ans

    def lift_to_precision(self, absprec):
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        return self.lift_to_precision_c((<Integer>absprec).value)

    cdef pAdicCappedRelativeElement lift_to_precision_c(pAdicCappedRelativeElement self, mpz_t absprec):
        cdef mpz_t relprec
        cdef int irelprec, usign
        cdef pAdicCappedRelativeElement ans
        usign = mpz_sgn(self.unit)
        if usign == -1:
            return self
        elif usign == 0:
            if mpz_cmp(self.ordp, absprec) >= 0:
                return self
            else:
                ans = self._new_c()
                ans.set_inexact_zero(absprec)
                return ans
        mpz_init(relprec)
        mpz_sub(relprec, self.ordp, absprec)
        if mpz_sgn(relprec) <= 0:
            mpz_clear(relprec)
            return self
        if mpz_cmp(relprec, (<Integer>self.parent().precision_cap()).value) > 0:
            irelprec = mpz_get_si((<Integer>self.parent().precision_cap()).value)
        else:
            irelprec = mpz_get_si(relprec)
        mpz_clear(relprec)
        if irelprec <= self.relprec:
            return self
        else:
            ans = self._new_c()
            mpz_set(ans.unit, self.unit)
            ans.set_precs(irelprec)
            mpz_set(ans.ordp, self.ordp)
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
            ulist = self.base_p_list(self.unit, self.p, list_mode, self.parent().prime_pow, mpz_get_si((<Integer>self.parent().precision_cap()).value))
        if self.in_field == 0 and mpz_sgn(self.ordp) == 1:
            zero = self._new_c()
            zero.set_exact_zero()
            ordp = PY_NEW(Integer)
            mpz_set(ordp.value, self.ordp)
            return [zero]*ordp + ulist
        else:
            return ulist

    cdef teichmuller_list(pAdicCappedRelativeElement self):
        # May eventually want to add a dict to store teichmuller lifts already seen, if p small enough
        cdef int curpower, preccap
        cdef mpz_t ppow, tmp
        cdef pAdicCappedRelativeElement list_elt
        cdef PowComputer_class powerer
        powerer = <PowComputer_class>self.parent().prime_pow
        ans = PyList_New(0)
        preccap = mpz_get_si((<Integer>self.parent().precision_cap()).value)
        curpower = preccap
        mpz_init(ppow)
        mpz_init_set(tmp, self.unit)
        while mpz_sgn(tmp) != 0:
            curpower -= 1
            list_elt = self._new_c()
            mpz_mod(list_elt.unit, tmp, self.p)
            if mpz_sgn(list_elt.unit) == 0:
                list_elt.set_exact_zero()
                mpz_divexact(tmp, tmp, self.p)
                powerer.pow_mpz_ui(ppow, curpower)
            else:
                mpz_set_ui(list_elt.ordp, 0)
                # ??? mpz_set_precs(preccap)
                sage.rings.padics.padic_generic_element.teichmuller_set_c(list_elt.unit, self.p, list_elt.modulus)
                mpz_sub(tmp, tmp, list_elt.unit)
                mpz_divexact(tmp, tmp, self.p)
                powerer.pow_mpz_ui(ppow, curpower)
                mpz_mod(tmp, tmp, ppow)
            PyList_Append(ans, list_elt)
        mpz_clear(ppow)
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
        L = self.list()
        if self.in_field == 1:
            n -= self.valuation()
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
        mpz_add_ui(ans.value, self.ordp, self.relprec)
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
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.relprec)
        return ans

    def residue(self, prec):
        """
        Reduces this mod $p^prec$
        INPUT:
            self -- a p-adic element
            prec - an integer
        OUTPUT:
            element of Z/(p^prec Z) -- self reduced mod p^prec
        EXAMPLES:
            sage: R = Zp(7,4,'capped-rel'); a = R(8); a.residue(1)
            1
            sage: R = Qp(7,4,'capped-rel'); a = R(8); a.residue(1)
            1
        """
        cdef Integer selfvalue, modulus
        cdef mpz_t ppow
        cdef PowComputer_class powerer
        if prec > self.precision_absolute():
            raise PrecisionError, "Not enough precision known in order to compute residue."
        elif prec < 0:
            raise ValueError, "cannot reduce modulo a negative power of p"
        elif self.valuation() >= 0:
            # Need to do this better.
            selfvalue = PY_NEW(Integer)
            modulus = PY_NEW(Integer)
            mpz_init(ppow)
            powerer = self.parent().prime_pow
            powerer.pow_mpz_mpz(selfvalue.value, self.ordp)
            mpz_mul(selfvalue.value, selfvalue.value, self.unit)
            if not PY_TYPE_CHECK(prec, Integer):
                prec = Integer(prec)
            powerer.pow_mpz_mpz(modulus.value, (<Integer>prec).value)
            return Mod(selfvalue, modulus)
        else:
            raise ValueError, "Element must have non-negative valuation in order to compute residue."

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
        if mpz_sgn(self.unit) == -1:
            raise ValueError, "unit part of 0 not defined"
        cdef pAdicCappedRelativeElement ans
        ans = self._new_c()
        mpz_set(ans.unit, self.unit)
        ans.set_precs(self.relprec)
        mpz_set_ui(ans.ordp, 0)
        return ans

    #def _unit_part(self):
    #    """
    #    Returns the unit part of self.
    #
    #    INPUT:
    #       self -- a p-adic element
    #    OUTPUT:
    #       the unit part of self -- a p-adic element
    #
    #    EXAMPLES:
    #    sage: R = Zp(17, 4,'capped-rel')
    #    """
    #    return self._unit

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
            ans = PY_NEW(Integer)
            mpz_set(ans.value, self.ordp)
            return ans

    cdef val_unit_c(self):
        cdef Integer val
        cdef pAdicCappedRelativeElement unit
        cdef mpz_t zero
        if mpz_sgn(self.unit) == -1:
            unit = self._new_c()
            mpz_init_set_ui(zero, 0)
            unit.set_inexact_zero(zero)
            mpz_clear(zero)
            return (infinity, unit)
        unit = self._new_c()
        mpz_set(unit.unit, self.unit)
        unit.set_precs(self.relprec)
        mpz_set_ui(unit.ordp, 0)
        val = PY_NEW(Integer)
        mpz_set(val.value, self.ordp)
        return (val, unit)

    def val_unit(self):
        return self.val_unit_c()

    def __hash__(self):
        return self._hash()

    cdef long _hash(self) except -1:
        cdef Integer ans
        if mpz_sgn(self.unit) <= 0:
            return 0
        else:
            ans = PY_NEW(Integer)
            mpz_xor(ans.value, self.unit, self.ordp)
            mpz_xor(ans.value, ans.value, self.modulus)
            return hash(ans)

    def _teichmuller_set(self, Integer n, Integer absprec):
        cdef mpz_t ppow
        if mpz_divisible_p(n.value, self.p):
            self.set_exact_zero()
        else:
            mpz_set(self.unit, n.value)
            mpz_set_ui(self.ordp, 0)
            if mpz_fits_sint_p(absprec.value) == 0:
                raise ValueError, "cannot compute teichmuller lift to that high precision"
            if mpz_sgn(absprec.value) != 1:
                raise ValueError, "can only compute to positive precision"
            self.set_precs(mpz_get_si(absprec.value))
            sage.rings.padics.padic_generic_element.teichmuller_set_c(self.unit, self.p, self.modulus)

def unpickle_pcre_v1(R, unit, ordp, relprec):
    cdef pAdicCappedRelativeElement ans
    ans = PY_NEW(pAdicCappedRelativeElement)
    ans._parent = R
    if R.is_field():
        ans.in_field = 1
    else:
        ans.in_field = 0
    mpz_init_set(ans.p, (<Integer>R.prime()).value)
    mpz_init_set(ans.unit, (<Integer>unit).value)
    mpz_init_set(ans.ordp, (<Integer>ordp).value)
    cdef PowComputer_class powerer
    powerer = R.prime_pow
    powerer.pow_mpz_mpz(ans.modulus, (<Integer>relprec).value)
    return ans
