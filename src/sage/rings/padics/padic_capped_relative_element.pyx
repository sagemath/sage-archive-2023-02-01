"""
Elements of p-Adic Rings with Capped Relative Precision

AUTHOR:
    -- David Roe
    -- Genya Zaytman: documentation
    -- David Harvey: doctests

TESTS:
    sage: M = MatrixSpace(pAdicField(3,100),2)
    sage: (M([1,0,0,90]) - (1+O(3^100)) * M(1)).left_kernel()
    Vector space of degree 2 and dimension 1 over 3-adic Field with capped relative precision 100
    Basis matrix:
    [1 + O(3^100)            0]
"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../../libs/ntl/decl.pxi"
include "../../ext/gmp.pxi"
include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"
include "../../libs/pari/decl.pxi"

cdef long maxint
from sys import maxint

cimport sage.rings.padics.padic_generic_element
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement

cimport sage.rings.rational
cimport sage.rings.padics.pow_computer
from sage.rings.padics.pow_computer cimport PowComputer_base
from sage.rings.padics.padic_printing cimport pAdicPrinter_class
from sage.rings.rational cimport Rational

import sage.rings.padics.padic_generic_element
import sage.rings.padics.padic_lazy_element
import sage.rings.integer_mod
import sage.rings.integer
import sage.rings.rational

cdef object infinity
from sage.rings.infinity import infinity
from sage.rings.integer_mod import Mod
from sage.rings.padics.precision_error import PrecisionError

from sage.rings.padics.padic_lazy_element import pAdicLazyElement

cdef PariInstance P = sage.libs.pari.all.pari

cdef long maxordp = (1L << (sizeof(long) * 8 - 2)) -1

cdef extern from "convert.h":
    cdef void t_INT_to_ZZ( mpz_t value, GEN g )

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

        Construct from Pari objects:
            sage: R = Zp(5)
            sage: x = pari(123123) ; R(x)
            3 + 4*5 + 4*5^2 + 4*5^3 + 5^4 + 4*5^5 + 2*5^6 + 5^7 + O(5^20)
            sage: R(pari(R(5252)))
            2 + 2*5^3 + 3*5^4 + 5^5 + O(5^20)
            sage: R = Zp(5,prec=5)
            sage: R(pari(-1))
            4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5)
            sage: pari(R(-1))
            4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5)
            sage: pari(R(0))
            0
            sage: R(pari(R(0) + O(5^5)))
            O(5^5)

        # todo: doctests for converting from other types of p-adic rings

        """
        #print "x = %s, type = %s, absprec = %s, relprec = %s"%(x, type(x),absprec, relprec)
        cpdef RingElement ordp
        cdef mpz_t modulus, tmp2
        cdef GEN pari_tmp
        cdef Integer tmp
        cdef unsigned long k
        mpz_init(self.unit)
        pAdicBaseGenericElement.__init__(self, parent)
        if empty:
            self._normalized = 0
            return
        self._normalized = 1
        if relprec is not infinity and not PY_TYPE_CHECK(relprec, Integer):
            relprec = Integer(relprec)
        if (relprec is infinity) or (relprec > parent.precision_cap()):
            relprec = parent.precision_cap()
        if not absprec is infinity and not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if isinstance(x, pAdicGenericElement):
            if self.prime_pow.in_field == 0 and x.valuation() < 0:
                raise ValueError, "element has negative valuation."
            if parent.prime() != x.parent().prime():
                raise TypeError, "Cannot coerce between p-adic parents with different primes"
        if isinstance(x, pAdicLazyElement):
            ## One can do this in a better way to minimize the amount of
            ## increasing precision on x.
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

            if (relprec is infinity) or (x.precision_relative() < relprec):
                try:
                    x.set_precision_relative(relprec)
                except PrecisionError:
                    pass

        if PY_TYPE_CHECK(x, Integer):
            self._set_from_Integer(x, absprec, relprec)
            return

        elif PY_TYPE_CHECK(x, Rational):
            self._set_from_Rational(x, absprec, relprec)
            return

        elif isinstance(x, (int, long)):
            x = Integer(x)
            self._set_from_Integer(x, absprec, relprec)
            return

        elif PY_TYPE_CHECK(x, pAdicBaseGenericElement):
            ## this case could be rethought to remove the
            ## potential arithmetic with infinity.
            ordp = x.valuation()
            if (ordp is infinity) and (absprec is infinity):
                self._set_exact_zero()
            elif ordp is infinity or ordp >= absprec:
                self._set_inexact_zero(mpz_get_si((<Integer>absprec).value))
            else:
                unit = x.unit_part().lift()
                if absprec is infinity:
                    relprec = min(relprec, x.precision_relative(), self.parent().precision_cap())
                else:
                    relprec = min(relprec, absprec - ordp, x.precision_relative(), self.parent().precision_cap())
                if ordp < 0 and self.prime_pow.in_field == 0:
                    raise ValueError, "negative valuation"
                self._set(mpz_get_si((<Integer>ordp).value), (<Integer>unit).value, mpz_get_si((<Integer>relprec).value))
                self._normalized = 1
            return

        elif sage.rings.integer_mod.is_IntegerMod(x):
            mpz_init_set(modulus, (<Integer>x.modulus()).value)
            k = mpz_remove(modulus, modulus, self.prime_pow.prime.value)
            if mpz_cmp_ui(modulus, 1) == 0:
                tmp = PY_NEW(Integer)
                mpz_set_ui(tmp.value, k)
                absprec = Integer(min(tmp, absprec))
                x = x.lift()
                mpz_clear(modulus)
                self._set_from_Integer(x, absprec, relprec)
                return
            else:
                mpz_clear(modulus)
                raise TypeError, "cannot coerce from the given integer mod ring (not a power of the same prime)"

        elif isinstance(x, pari_gen):

            pari_tmp = (<pari_gen>x).g

            #if x.type() == "t_PADIC":
            if typ(pari_tmp) == t_PADIC:
                self.relprec = precp(pari_tmp)
                if self.relprec > relprec:
                    self.relprec = relprec
                self.ordp = valp(pari_tmp)
                t_INT_to_ZZ(self.unit, <GEN>(pari_tmp[4]))
                self._normalized = 1
                if mpz_sgn(self.unit) == -1 or mpz_cmp(self.unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0]) >= 0:
                    mpz_mod(self.unit, self.unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0])
                if mpz_sgn(self.unit) == 0:
                    self._set_inexact_zero(self.ordp + self.relprec)

            elif typ(pari_tmp) == t_INT:
                mpz_init(tmp2)
                t_INT_to_ZZ(tmp2, pari_tmp)
                ### This code is duplicated from _set_from_Integer
                if absprec is infinity:
                    if relprec is infinity or mpz_fits_slong_p((<Integer>relprec).value) == 0:
                        self._set_from_mpz_rel(tmp2, self.prime_pow.prec_cap)
                    else:
                        self._set_from_mpz_rel(tmp2, mpz_get_si((<Integer>relprec).value))
                else:
                    if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                        raise ValueError, "absprec must fit in a long"
                    if relprec is infinity or mpz_fits_slong_p((<Integer>relprec).value) == 0:
                        self._set_from_mpz_both(tmp2, mpz_get_si((<Integer>absprec).value), self.prime_pow.prec_cap)
                    else:
                        self._set_from_mpz_both(tmp2, mpz_get_si((<Integer>absprec).value), mpz_get_si((<Integer>relprec).value))

            elif x.type() == "t_FRAC":
                self._set_from_Rational(Rational(x), absprec, relprec)

            else:
                raise TypeError, "unsupported coercion from pari: only p-adics, integers and rationals allowed"

        ##
        ## If this case gets included, it should be rewritten to
        ## have the following new ending:
        ##
        #if sage.rings.finite_field_element.is_FiniteFieldElement(x):
        #    if x.parent().order() != parent.prime():
        #        raise TypeError, "can only create p-adic element out of finite field when order of field is p"
        #    #prec = min(prec, 1)
        #    x = x.lift()
        ##  self._set_from_Integer(parent,x,absprec,relprec)

        else:
            self._set_from_Rational(Rational(x), absprec, relprec)

    cdef int _set_from_Integer(pAdicCappedRelativeElement self,
                               Integer x, absprec, relprec) except -1:
        """
        self.prime_pow should already be set.
        absprec should be infinity or an Integer.
        relprec should be infinity or a nonnegative Integer.
        """
        if absprec is infinity:
            if relprec is infinity or mpz_fits_slong_p((<Integer>relprec).value) == 0:
                return self._set_from_mpz_rel(x.value, self.prime_pow.prec_cap)
            else:
                return self._set_from_mpz_rel(x.value, mpz_get_si((<Integer>relprec).value))
        else:
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                raise ValueError, "absprec must fit in a long"
            if relprec is infinity or mpz_fits_slong_p((<Integer>relprec).value) == 0:
                return self._set_from_mpz_both(x.value, mpz_get_si((<Integer>absprec).value), self.prime_pow.prec_cap)
            else:
                return self._set_from_mpz_both(x.value, mpz_get_si((<Integer>absprec).value), mpz_get_si((<Integer>relprec).value))

    cdef int _set_from_Rational(pAdicCappedRelativeElement self,
                                Rational x, absprec, relprec) except -1:
        """
        self.prime_pow should already be set.
        absprec should be infinity or an Integer.
        relprec should be infinity or a nonnegative Integer.
        """
        if absprec is infinity:
            if relprec is infinity or mpz_fits_slong_p((<Integer>relprec).value) == 0:
                return self._set_from_mpq_rel(x.value, self.prime_pow.prec_cap)
            else:
                return self._set_from_mpq_rel(x.value, mpz_get_si((<Integer>relprec).value))
        else:
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                raise ValueError, "absprec must fit in a long"
            if relprec is infinity or mpz_fits_slong_p((<Integer>relprec).value) == 0:
                return self._set_from_mpq_both(x.value, mpz_get_si((<Integer>absprec).value), self.prime_pow.prec_cap)
            else:
                return self._set_from_mpq_both(x.value, mpz_get_si((<Integer>absprec).value), mpz_get_si((<Integer>relprec).value))

    cdef int _set_from_mpz_rel(pAdicCappedRelativeElement self, mpz_t value, long relprec) except -1:
        """
        self.prime_pow should already be set.
        relprec should be in the range [0,self.prime_pow.prec_cap]
        """
        if mpz_sgn(value) == 0:
            self._set_exact_zero()
            return 0
        _sig_on
        self.ordp = mpz_remove(self.unit, value, self.prime_pow.prime.value)
        self.relprec = relprec
        _sig_off
        if mpz_sgn(self.unit) == -1 or \
               (mpz_cmp(self.unit, self.prime_pow.pow_mpz_t_tmp(relprec)) >= 0):
            _sig_on
            mpz_mod(self.unit, self.unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0])
            _sig_off
        self._normalized = 1
        return 0

    cdef int _set_from_mpz_both(pAdicCappedRelativeElement self, mpz_t value, long absprec, long relprec) except -1:
        """
        self.prime_pow should already be set.
        relprec should be in the range [0,self.prime_pow.prec_cap]
        """
        if mpz_sgn(value) == 0:
            self._set_inexact_zero(absprec)
            return 0
        _sig_on
        self.ordp = mpz_remove(self.unit, value, self.prime_pow.prime.value)
        _sig_off
        if self.ordp >= absprec:
            self._set_inexact_zero(absprec)
            return 0
        self.relprec = absprec - self.ordp
        if self.relprec > relprec:
            self.relprec = relprec
        if mpz_sgn(self.unit) == -1 or \
               (mpz_cmp(self.unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)) >= 0):
            _sig_on
            mpz_mod(self.unit, self.unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0])
            _sig_off
        self._normalized = 1
        return 0

    cdef int _set_from_mpq_rel(pAdicCappedRelativeElement self, mpq_t value, long relprec) except -1:
        """
        self.prime_pow should already be set.
        relprec should be in the range [0,self.prime_pow.prec_cap]
        """
        cdef mpz_t tmp
        if mpq_sgn(value) == 0:
            self._set_exact_zero()
            return 0
        _sig_on
        self.ordp = mpz_remove(self.unit, mpq_numref(value), self.prime_pow.prime.value)
        _sig_off
        if self.ordp == 0:
            self.ordp = mpz_remove(self.unit, mpq_denref(value), self.prime_pow.prime.value)
            if self.ordp > 0 and self.prime_pow.in_field == 0:
                raise ValueError, "p divides the denominator"
            self.ordp = -self.ordp
            mpz_invert(self.unit, self.unit, self.prime_pow.pow_mpz_t_tmp(relprec)[0])
            mpz_mul(self.unit, self.unit, mpq_numref(value))
            mpz_mod(self.unit, self.unit, self.prime_pow.pow_mpz_t_tmp(relprec)[0])
        else:
            mpz_init(tmp)
            mpz_invert(tmp, mpq_denref(value), self.prime_pow.pow_mpz_t_tmp(relprec)[0])
            mpz_mul(self.unit, self.unit, tmp)
            mpz_mod(self.unit, self.unit, self.prime_pow.pow_mpz_t_tmp(relprec)[0])
            mpz_clear(tmp)
        self.relprec = relprec
        self._normalized = 1
        return 0

    cdef int _set_from_mpq_both(pAdicCappedRelativeElement self, mpq_t value, long absprec, long relprec) except -1:
        """
        self.prime_pow should already be set.
        relprec should be in the range [0,self.prime_pow.prec_cap]
        """
        cdef mpz_t num_unit, den_unit
        cdef long num_ordp, den_ordp
        if mpq_sgn(value) == 0:
            self._set_inexact_zero(absprec)
            return 0
        _sig_on
        mpz_init(num_unit)
        mpz_init(den_unit)
        num_ordp = mpz_remove(num_unit, mpq_numref(value), self.prime_pow.prime.value)
        den_ordp = mpz_remove(den_unit, mpq_denref(value), self.prime_pow.prime.value)
        self.ordp = num_ordp - den_ordp
        _sig_off
        if self.ordp < 0 and self.prime_pow.in_field == 0:
            mpz_clear(num_unit)
            mpz_clear(den_unit)
            raise ValueError, "p divides the denominator"
        if self.ordp >= absprec:
            mpz_clear(num_unit)
            mpz_clear(den_unit)
            self._set_inexact_zero(absprec)
            return 0
        self.relprec = absprec - self.ordp
        if self.relprec > relprec:
            self.relprec = relprec
        if mpz_sgn(num_unit) == -1 or mpz_cmp(num_unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0]) >= 0:
            mpz_mod(num_unit, num_unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0])
        if mpz_cmp_ui(den_unit, 1) == 0:
            mpz_set(self.unit, num_unit)
        else:
            mpz_invert(self.unit, den_unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0])
            mpz_mul(self.unit, self.unit, num_unit)
            mpz_mod(self.unit, self.unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0])
        mpz_clear(num_unit)
        mpz_clear(den_unit)
        self._normalized = 1

    cdef int _set_from_CR(pAdicCappedRelativeElement self, pAdicCappedRelativeElement other) except -1:
        mpz_set(self.unit, other.unit)
        if mpz_sgn(self.unit) >= 0:
            self.ordp = other.ordp
            if self.ordp < 0 and self.prime_pow.in_field == 0:
                raise ValueError, "negative valuation"
            self.relprec = other.relprec
            if self.relprec > self.prime_pow.prec_cap:
                if not other._normalized:
                    other._normalize()
                    self.ordp = other.ordp
                    self.relprec = other.relprec
                    mpz_set(self.unit, other.unit)
                if self.relprec > self.prime_pow.prec_cap and mpz_sgn(self.unit) == 1 and mpz_cmp(self.unit, self.prime_pow.pow_mpz_t_top()[0]) >= 0:
                    _sig_on
                    mpz_mod(self.unit, self.unit, self.prime_pow.pow_mpz_t_top()[0])
                    _sig_off
            self._normalized = other._normalized
        else:
            self._normalized = 1
        return 0

    cdef int _set_exact_zero(pAdicCappedRelativeElement self) except -1:
        mpz_set_si(self.unit, -1)
        self._normalized = 1

    cdef int _set_inexact_zero(pAdicCappedRelativeElement self, long absprec) except -1:
        self.relprec = 0
        mpz_set_ui(self.unit, 0)
        self.ordp = absprec
        self._normalized = 1

    cdef int _set_zero(pAdicCappedRelativeElement self, absprec) except -1:
        if absprec is infinity:
            mpz_set_si(self.unit, -1)
        else:
            mpz_set_ui(self.unit, 0)
            self.relprec = 0
            self.ordp = mpz_get_si((<Integer>absprec).value)
        self._normalized = 1

    cdef int _set_prec(pAdicCappedRelativeElement self, long relprec) except -1:
        self.relprec = relprec

    cdef int _set(pAdicCappedRelativeElement self, long ordp, mpz_t unit, long relprec) except -1:
        self.relprec = relprec
        self.ordp = ordp
        mpz_set(self.unit, unit)
        if mpz_sgn(self.unit) == -1 or mpz_cmp(self.unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0]) >= 0:
            mpz_mod(self.unit, self.unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0])
        if mpz_sgn(self.unit) == 0:
            self._set_inexact_zero(self.ordp + self.relprec)

    cdef int _set_to_mpz(pAdicCappedRelativeElement self, mpz_t dest) except -1:
        if mpz_sgn(self.unit) == -1:
            mpz_set_ui(dest, 0)
        elif self.ordp < 0:
            raise ValueError, "negative valuation"
        elif mpz_sgn(self.unit) == 0:
            mpz_set(dest, self.prime_pow.pow_mpz_t_tmp(self.ordp))
        else:
            mpz_set(dest, self.unit)
            if self.ordp > 0:
                mpz_mul(dest, dest, self.prime_pow.pow_mpz_t_tmp(self.ordp))
        return 0

    cdef int _set_to_mpq(pAdicCappedRelativeElement self, mpq_t dest) except -1:
        if mpz_sgn(self.unit) == -1:
            mpq_set_ui(dest, 0, 1)
        elif self.ordp < 0:
            mpz_set(mpq_denref(dest), self.prime_pow.pow_mpz_t_tmp(-self.ordp))
            mpz_set(mpq_numref(dest), self.unit)
        elif mpz_sgn(self.unit) == 0:
            mpq_set_z(dest, self.prime_pow.pow_mpz_t_tmp(self.ordp))
        else:
            mpq_set_z(dest, self.unit)
            if self.ordp > 0:
                mpz_mul(mpq_numref(dest), mpq_numref(dest), self.prime_pow.pow_mpz_t_tmp(self.ordp))
        return 0

    cdef pAdicCappedRelativeElement _new_c(pAdicCappedRelativeElement self):
        cdef pAdicCappedRelativeElement ans
        ans = PY_NEW(pAdicCappedRelativeElement)
        ans._parent = self._parent
        ans.prime_pow = self.prime_pow
        mpz_init(ans.unit)
        return ans

    cdef int _normalize(pAdicCappedRelativeElement self) except -1:
        cdef long diff
        if self._normalized == 0:
            if mpz_sgn(self.unit) > 0:
                if mpz_divisible_p(self.unit, self.prime_pow.prime.value) != 0:
                    diff = mpz_remove(self.unit, self.unit, self.prime_pow.prime.value)
                    if self.relprec > diff:
                        self.ordp += diff
                        self.relprec -= diff
                    else:
                        self._set_inexact_zero(self.ordp + self.relprec)
                if mpz_cmp(self.unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0]) >= 0:
                    mpz_mod(self.unit, self.unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0])
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
        cdef Integer unit
        unit = <Integer>PY_NEW(Integer)
        mpz_set(unit.value, self.unit)
        return unpickle_pcre_v1, (self.parent(), unit, self.ordp, self.relprec)

    def __richcmp__(left, right, int op):
        return (<Element>left)._richcmp(right, op)

    cpdef ModuleElement _neg_(self):
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
        mpz_sub(ans.unit, ans.prime_pow.pow_mpz_t_tmp(self.relprec)[0], self.unit)
        return ans

    def __pow__(pAdicCappedRelativeElement self, _right, dummy):
        r"""
        Returns self^right.

        Note: when right is divisible by p then one can get more precision than expected.
        Lemma 2.1 (Constructing Class Fields over Local Fields, Sebastian Pauli):
        [modified from original for Qp.  See padic_ZZ_pX_CR_element for original]
        Let $\alpha$ be in $\mathbb{Z}_p$. The $p$-th power of $1 + \alpha p^{\lambda}$ satisifes
            (1 + \alpha p^{\lambda})^p \equiv 1 + \alpha p^{\lambda + 1} mod p^{\lambda + 2}

            unless $\lambda = 1$ and $p = 2$, in which case
            (1 + 2 \alpha)^2 \equiv 1 + 4(\alpha^2 + \alpha) mod 8

        So for $p \ne 2$, if right is divisible by $p^k$ then we add $k$ to the relative precision
        of the answer.

        For $p = 2$, if we start with something of relative precision 1 (ie $2^m + O(2^{m+1})$),
        $\alpha^2 + \alpha \equiv 0 \mod 2$, so the precision of the result is $k + 2$:
        $(2^m + O(2^{m+1}))^{2^k} = 2^{m 2^k} + O(2^{m 2^k + k + 2})

        There is also the issue of $p$-adic exponents, and determining how the precision of the exponent affects the precision of the result.
        In computing $(a + O(p^k))^{b + O(p^m)}$, we can factor out the Teichmuller part and use the above lemma to find the first spot where
        $(1 + \alpha p^{\lambda})^(p^m)$ differs from 1.  This a relative precision of $\lambda + m$ except in the case $p = 2$ and $\lambda = 1$,
        where it gives $m + 2$.  We compare this with the precision bound given by computing $(a + O(p^k))^b$ (ie $k + b.valuation(p)$
        or $2 + b.valuation(2)$ if $p = 2$ and $k = 1$) and take the lesser of the two.

        In order to do this we need to compute the valuation of (self / self.parent().teichmuller(self)) - 1.  This takes a reasonable amount of time: we cache the result
        as __pow_level.

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
            sage: K(5, 3)^19
            5 + 3*19 + 11*19^3 + O(19^4)
        """
        self._normalize()
        cdef pAdicCappedRelativeElement ans
        cdef mpz_t tmp
        cdef Integer right
        cdef Integer exp_val, exp_prec
        cdef Integer base_level
        cdef pAdicCappedRelativeElement base
        cdef bint padic_exp
        cdef long relprec
        if mpz_sgn(self.unit) < 0:
            # Return 0 except for 0^0 error or type error on the exponent.
            if PY_TYPE_CHECK(_right, Integer) or PY_TYPE_CHECK(_right, Rational) or (PY_TYPE_CHECK(_right, pAdicBaseGenericElement) and _right.parent().prime() == self.prime_pow.prime)  or isinstance(_right, (int, long)):
                if _right == 0:
                    raise ArithmeticError, "0^0 is undefined"
                return self
            else:
                raise TypeError, "exponent must be an integer, rational or base p-adic with the same prime"
        elif mpz_sgn(self.unit) == 0:
            # If an integer exponent, return an inexact zero of valuation right * self.ordp.  Otherwise raise an error.
            if isinstance(_right, (int, long)):
                _right = Integer(_right)
            if PY_TYPE_CHECK(_right, Integer):
                ans = self._new_c()
                mpz_init_set_si(tmp, self.ordp)
                mpz_mul(tmp, tmp, (<Integer>_right).value)
                if mpz_cmp_si(tmp, maxordp) >= 0 or mpz_cmp_si(tmp, -maxordp) <= 0:
                    raise ValueError, "valuation overflow"
                ans._set_inexact_zero(mpz_get_si(tmp))
                mpz_clear(tmp)
                return ans
            elif PY_TYPE_CHECK(_right, Rational) or (PY_TYPE_CHECK(_right, pAdicBaseGenericElement) and _right.parent().prime() == self.prime_pow.prime):
                raise ValueError, "Need more precision"
            else:
                raise TypeError, "exponent must be an integer, rational or base p-adic with the same prime"
        if isinstance(_right, (int, long)):
            _right = Integer(_right)
        if PY_TYPE_CHECK(_right, Integer):
            right = <Integer> _right
            if right == 0:
                # return 1 to maximum precision
                ans = self._new_c()
                ans.ordp = 0
                ans.relprec = self.prime_pow.prec_cap
                ans._normalized = True
                mpz_set_ui(ans.unit, 1)
                return ans
            exp_val = right.valuation(self.prime_pow.prime)
            padic_exp = False
        elif PY_TYPE_CHECK(_right, pAdicBaseGenericElement) and _right.parent().prime() == self.prime_pow.prime:
            if self.ordp != 0:
                raise ValueError, "in order to raise to a p-adic exponent, base must be a unit"
            right = Integer(_right)
            padic_exp = True
            exp_prec = _right.precision_absolute()
            exp_val = _right.valuation()
            if exp_val < 0:
                raise NotImplementedError, "negative valuation exponents not yet supported"
            try:
                base_level = self.__pow_level
            except AttributeError:
                # compute the "level"
                teich_part = self.parent().teichmuller(self)
                base_level = (self / teich_part - 1).valuation() ##
                self.__pow_level = base_level
        elif PY_TYPE_CHECK(_right, Rational):
            raise NotImplementedError
        else:
            raise TypeError, "exponent must be an integer, rational or base p-adic with the same prime"
        # if right < 0, we return (~self)^(-right)
        # Now we compute the increased relprec due to the exponent having positive p-adic valuation
        if exp_val > 0:
            mpz_init_set_si(tmp, self.relprec)
            mpz_add(tmp, tmp, exp_val.value)
            if mpz_cmp_ui(self.prime_pow.prime.value, 2) == 0 and self.relprec == 1:
                mpz_add_ui(tmp, tmp, 1)
            if mpz_cmp_si(tmp, self.prime_pow.prec_cap) > 0:
                relprec = self.prime_pow.prec_cap
            else:
                relprec = mpz_get_si(tmp)
            mpz_clear(tmp)
        else:
            relprec = self.relprec
        # Now we compute the limit on relprec due to a non-infinite precision on the exponent.
        if padic_exp:
            if exp_prec > 0:
                # I can change base_level, so I use it in place of tmp above.
                if mpz_cmp_ui(self.prime_pow.prime.value, 2) == 0 and base_level == 1:
                    mpz_add_ui(base_level.value, base_level.value, 1)
                mpz_add(base_level.value, base_level.value, exp_prec.value)
                if mpz_cmp_si(base_level.value, relprec) < 0:
                    relprec = mpz_get_si(base_level.value)
            else:
                ans = self._new_c()
                ans._set_inexact_zero(0)
                return ans
        if right < 0:
            base = ~self
            right = -right
        else:
            base = self
        ans = self._new_c()
        ans._set_prec(relprec)
        mpz_init_set_si(tmp, base.ordp)
        mpz_mul(tmp, right.value, tmp)
        if mpz_cmp_si(tmp, maxordp) >= 0:
            raise ValueError, "Valuation too large"
        ans.ordp = mpz_get_si(tmp)
        mpz_clear(tmp)
        _sig_on
        mpz_powm(ans.unit, base.unit, right.value, ans.prime_pow.pow_mpz_t_tmp(ans.relprec)[0])
        _sig_off
        return ans

    # Once the code for _add_ has stabilized, it may be worth getting rid of the extra function call for _sub_.

    cpdef ModuleElement _add_(self, ModuleElement _right):
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
        # TODO: add more examples/tests that verify all the cases below (w.r.t. precision, valuation)
        cdef pAdicCappedRelativeElement ans
        cdef pAdicCappedRelativeElement right = _right
        cdef long tmpL
        # If either self or right is an exact zero, we do the appropriate thing
        if mpz_sgn(self.unit) == -1:
            return right
        if mpz_sgn(right.unit) == -1:
            return self
        if self.ordp == right.ordp:
            ans = self._new_c()
            # The relative precision of the sum is the minimum of the relative precisions in this case, possibly decreasing if we got cancellation
            if self.relprec < right.relprec:
                ans.relprec = self.relprec
            else:
                ans.relprec = right.relprec
            # Since the valuations are the same, we can just add the units
            mpz_add(ans.unit, self.unit, right.unit)
            ans.ordp = self.ordp
            ans._normalized = 0
        else:
            if self.ordp > right.ordp:
                # Addition is commutative, swap so self.ordp < right.ordp
                ans = right; right = self; self = ans
            tmpL = right.ordp - self.ordp
            if tmpL >= self.relprec:
                return self
            ans = self._new_c()
            ans.ordp = self.ordp
            if mpz_size(self.prime_pow.pow_mpz_t_top()[0]) > 10000:
                # We only enable the signal handler if the product will take a while.
                _sig_on
                mpz_mul(ans.unit, right.unit, self.prime_pow.pow_mpz_t_tmp(tmpL)[0])
                _sig_off
            else:
                mpz_mul(ans.unit, right.unit, self.prime_pow.pow_mpz_t_tmp(tmpL)[0])
            mpz_add(ans.unit, ans.unit, self.unit)
            if self.relprec <= tmpL + right.relprec:
                ans._set_prec(self.relprec)
            else:
                ans._set_prec(tmpL + right.relprec)
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
        ans._set_prec(self.relprec)
        _sig_on
        mpz_invert(ans.unit, self.unit, ans.prime_pow.pow_mpz_t_tmp(ans.relprec)[0])
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

    cpdef RingElement _floordiv_c_impl(self, RingElement right):
        cdef pAdicCappedRelativeElement ans
        cdef long relprec, diff
        (<pAdicCappedRelativeElement>right)._normalize()
        self._normalize()
        # For fields, we define floor division as normal division and % as always 0
        if self.prime_pow.in_field == 1:
            return self._div_(right)
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
        mpz_invert(ans.unit, (<pAdicCappedRelativeElement>right).unit, self.prime_pow.pow_mpz_t_tmp((<pAdicCappedRelativeElement>right).relprec)[0])
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
            ans._set_prec(relprec)
            if mpz_cmp(ans.unit, self.prime_pow.pow_mpz_t_tmp(relprec)[0]) >= 0:
                ans._normalized = 0
            else:
                ans._normalized = 1
        # Otherwise, we have to do floordivision on the unit part of the answer
        else:
            diff = (<pAdicCappedRelativeElement>right).ordp - self.ordp
            # If we're dividing by something where the difference in valuations is bigger than the relative precision, we can only get zero.
            if diff >= relprec:
                ans._set_inexact_zero(0)
            else:
                # Otherwise, our relative precision goes down by the difference in valuations, and we set ans.unit to ans.unit // ppow.
                relprec = relprec - diff
                _sig_on
                mpz_fdiv_q(ans.unit, ans.unit, self.prime_pow.pow_mpz_t_tmp(diff)[0])
                _sig_off
                ans.ordp = 0
                ans._normalized = 0
                ans._set_prec(relprec)
        return ans

    def __lshift__(self, shift):
        # TODO: move this up the hierarchy, perhaps this should go all the way to element?
        # The "verify that shift is an integer" part could be shared
        # should accept int (rather than do int -> Integer -> int every time)
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
        # TODO: move this up the hierarchy
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
                ans._set_inexact_zero(0)
            else:
                relprec = self.relprec - diff
                _sig_on
                mpz_fdiv_q(ans.unit, self.unit, self.prime_pow.pow_mpz_t_tmp(diff)[0])
                _sig_off
                ans.ordp = 0
                ans._set_prec(relprec)
                ans._normalized = 0
        return ans

    def _integer_(self, Z=None):
        r"""
        Converts self to an integer

        sage: R = Zp(5, prec = 4); a = R(642); a
        2 + 3*5 + O(5^4)
        sage: Integer(a)
        17
        """
        if self.ordp < 0:
            raise ValueError, "Cannot form an integer out of a p-adic field element with negative valuation"
        return self.lift_c()

    cpdef RingElement _mul_(self, RingElement _right):
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
        cdef pAdicCappedRelativeElement right = <pAdicCappedRelativeElement>_right
        if mpz_sgn(self.unit) == -1:
            return self
        if mpz_sgn(right.unit) == -1:
            return right
        self._normalize()
        right._normalize()
        ans = self._new_c()
        # Do we need to do overflow checking here?
        ans.ordp = self.ordp + right.ordp
        if self.relprec <= right.relprec:
            ans._set_prec(self.relprec)
        else:
            ans._set_prec(right.relprec)
        mpz_mul(ans.unit, self.unit, right.unit)
        if mpz_cmp(ans.unit, self.prime_pow.pow_mpz_t_tmp(ans.relprec)[0]) >= 0:
            ans._normalized = 0
        else:
            ans._normalized = 1
        return ans

    cpdef RingElement _div_(self, RingElement right):
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
                ans._set_prec(self.relprec)
            else:
                ans._set_prec((<pAdicCappedRelativeElement>right).relprec)
            _sig_on
            mpz_invert(ans.unit, (<pAdicCappedRelativeElement>right).unit, self.prime_pow.pow_mpz_t_tmp(ans.relprec)[0])
            mpz_mul(ans.unit, ans.unit, self.unit)
            _sig_off
            if mpz_cmp(ans.unit, self.prime_pow.pow_mpz_t_tmp(ans.relprec)[0]) >= 0:
                ans._normalized = 0
            else:
                ans._normalized = 1
        return ans

    def add_bigoh(self, absprec):
        """
        Returns a new element with absolute precision decreased to
        absprec.

        INPUT:
            self -- a p-adic element
            absprec -- an integer

        OUTPUT:
            element -- self with precision set to the minimum of
                       self's precision and absprec

        EXAMPLE:
            sage: R = Zp(7,4,'capped-rel','series'); a = R(8); a.add_bigoh(1)
            1 + O(7)
            sage: b = R(0); b.add_bigoh(3)
            O(7^3)
            sage: R = Qp(7,4); a = R(8); a.add_bigoh(1)
            1 + O(7)
            sage: b = R(0); b.add_bigoh(3)
            O(7^3)

        The precision never increases:
            sage: R(4).add_bigoh(2).add_bigoh(4)
            4 + O(7^2)

        Another example that illustrates that the precision does not
        increase:
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
        cdef pAdicCappedRelativeElement ans
        cdef long aprec, newprec
        if PY_TYPE_CHECK(absprec, int):
            aprec = absprec
        else:
            if not PY_TYPE_CHECK(absprec, Integer):
                absprec = Integer(absprec)
            aprec = mpz_get_si((<Integer>absprec).value)
        if mpz_sgn(self.unit) == -1 or aprec < self.ordp:
            ans = self._new_c()
            ans._set_inexact_zero(aprec)
            return ans
        # Do we still need to worry about overflow?
        if aprec > self.ordp + self.relprec:
            return self

        ans = self._new_c()
        ans.ordp = self.ordp
        newprec = aprec - self.ordp
        if newprec >= self.relprec:
            return self
        ans._set_prec(newprec)
        mpz_set(ans.unit, self.unit)
        if mpz_cmp(self.unit, self.prime_pow.pow_mpz_t_tmp(ans.relprec)[0]) >= 0:
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

    cpdef bint _is_exact_zero(self) except -1:
        return mpz_sgn(self.unit) == -1

    cpdef bint _is_inexact_zero(self) except -1:
        self._normalize()
        return self.relprec == 0 and not self._is_exact_zero()

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
            return mpz_sgn(self.unit) <= 0
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
        return self.ordp >= mpz_get_si((<Integer>absprec).value)

    def __nonzero__(self):
        self._normalize()
        return mpz_sgn(self.unit) > 0

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
        # TODO: lots of examples (this is a non-trivial function)
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
                    if mpz_cmp((<pAdicCappedRelativeElement>right).unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0]) >= 0:
                        mpz_init(tmp)
                        mpz_mod(tmp, (<pAdicCappedRelativeElement>right).unit, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0])
                        if mpz_cmp(tmp, self.unit) == 0:
                            mpz_clear(tmp)
                            return True
                        else:
                            mpz_clear(tmp)
                            return False
                    else:
                        return mpz_cmp(self.unit, (<pAdicCappedRelativeElement>right).unit) == 0
                else:
                    if mpz_cmp(self.prime_pow.pow_mpz_t_tmp((<pAdicCappedRelativeElement>right).relprec)[0], self.unit) <= 0:
                        mpz_init(tmp)
                        mpz_mod(tmp, self.unit, self.prime_pow.pow_mpz_t_tmp((<pAdicCappedRelativeElement>right).relprec)[0])
                        if mpz_cmp(tmp, (<pAdicCappedRelativeElement>right).unit) == 0:
                            mpz_clear(tmp)
                            return True
                        else:
                            mpz_clear(tmp)
                            return False
                    else:
                        return mpz_cmp(self.unit, (<pAdicCappedRelativeElement>right).unit) == 0
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        cdef long aprec
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
        mpz_mod(tmp, self.unit, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
        mpz_mod(tmp2, (<pAdicCappedRelativeElement>right).unit, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
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
        return self.lift_c()

    cdef lift_c(self):
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
            mpz_set(mpq_denref(ansr.value), self.prime_pow.pow_mpz_t_tmp(-self.ordp)[0])
            return ansr
        elif self.ordp == 0:
            ans = PY_NEW(Integer)
            mpz_set(ans.value, self.unit)
            return ans
        else:
            ans = PY_NEW(Integer)
            mpz_mul(ans.value, self.unit, self.prime_pow.pow_mpz_t_tmp(self.ordp)[0])
            return ans

    def lift_to_precision(self, absprec):
        # TODO: docs ans examples
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
                ans._set_inexact_zero(absprec)
                return ans
        relprec = absprec - self.ordp
        if relprec <= self.relprec:
            return self
        else:
            ans = self._new_c()
            mpz_set(ans.unit, self.unit)
            ans._set_prec(relprec)
            ans.ordp = self.ordp
            ans._normalized = self._normalized
            return ans

    cdef pari_gen _to_gen(pAdicCappedRelativeElement self):
        return P.new_gen_from_padic(self.ordp, self.relprec,
                                    self.prime_pow.prime.value,
                                    self.prime_pow.pow_mpz_t_tmp(self.relprec)[0],
                                    self.unit)

    def _pari_(self):
        if mpz_sgn(self.unit) == -1:
            return P.new_gen_from_int(0)
        else:
            return P.new_gen_from_padic(self.ordp, self.relprec,
                                        self.prime_pow.prime.value,
                                        self.prime_pow.pow_mpz_t_tmp(self.relprec)[0],
                                        self.unit)

    def list(self, lift_mode = 'simple'):
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
        if lift_mode == 'teichmuller':
            ulist = self.teichmuller_list()
        elif lift_mode == 'simple':
            ulist = (<pAdicPrinter_class>self.parent()._printer).base_p_list(self.unit, True)
        elif lift_mode == 'smallest':
            ulist = (<pAdicPrinter_class>self.parent()._printer).base_p_list(self.unit, False)
        else:
            raise ValueError
        if self.prime_pow.in_field == 0 and self.ordp > 0:
            zero = self._new_c()
            zero._set_exact_zero()
            ordp = PY_NEW(Integer)
            mpz_set_si(ordp.value, self.ordp)
            return [zero]*ordp + ulist
        else:
            return ulist

# Working on this
    cdef teichmuller_list(pAdicCappedRelativeElement self):
        # May eventually want to add a dict to store teichmuller lifts already seen, if p small enough
        cdef unsigned long curpower, preccap, i
        cdef mpz_t tmp, tmp2
        cdef pAdicCappedRelativeElement list_elt
        cdef PowComputer_class powerer
        ans = PyList_New(0)
        preccap = self.relprec
        curpower = preccap
        self._normalize()
        if mpz_sgn(self.unit) <= 0:
            return ans
        mpz_init_set(tmp, self.unit)
        mpz_init(tmp2)
        if self.prime_pow.in_field == 0 and self.ordp > 0:
            list_elt = self._new_c()
            list_elt._set_exact_zero()
            for i from 0 <= i < self.ordp:
                PyList_Append(ans, list_elt)
        while mpz_sgn(tmp) != 0:
            curpower -= 1
            list_elt = self._new_c()
            mpz_mod(list_elt.unit, tmp, self.prime_pow.prime.value)
            if mpz_sgn(list_elt.unit) == 0:
                list_elt._set_exact_zero()
                mpz_divexact(tmp, tmp, self.prime_pow.prime.value)
            else:
                list_elt.ordp = 0
                list_elt._set_prec(preccap)
                mpz_set(tmp2, self.prime_pow.pow_mpz_t_tmp(preccap)[0])
                self.teichmuller_set_c(list_elt.unit, tmp2)
                list_elt._normalized = 1
                mpz_sub(tmp, tmp, list_elt.unit)
                mpz_divexact(tmp, tmp, self.prime_pow.prime.value)
                mpz_mod(tmp, tmp, self.prime_pow.pow_mpz_t_tmp(curpower)[0])
            PyList_Append(ans, list_elt)
        mpz_clear(tmp)
        mpz_clear(tmp2)
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
        cdef long aprec
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
            mpz_set(modulus.value, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
            selfvalue = PY_NEW(Integer)
            mpz_set_ui(selfvalue.value, 0)
            return Mod(selfvalue, modulus)
        else:
            # Need to do this better.
            selfvalue = PY_NEW(Integer)
            modulus = PY_NEW(Integer)
            mpz_mul(selfvalue.value, self.prime_pow.pow_mpz_t_tmp(self.ordp)[0], self.unit)
            mpz_set(modulus.value, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
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
        ans._set_prec(self.relprec)
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
        cdef long v = self.valuation_c()
        cdef long maxlong
        if sizeof(long) == 4:
            maxlong = 2147483647
        elif sizeof(long) == 8:
            maxlong = 9223372036854775807
        else:
            raise RuntimeError
        if v == maxlong:
            return infinity
        cdef Integer ans = PY_NEW(Integer)
        mpz_set_si(ans.value, v)
        return ans

    cdef long valuation_c(self):
        cdef Integer ans
        cdef long maxlong
        self._normalize()
        if mpz_sgn(self.unit) == -1:
            if sizeof(long) == 4:
                maxlong = 2147483647
            elif sizeof(long) == 8:
                maxlong = 9223372036854775807
            else:
                raise RuntimeError
            return maxlong
        else:
            return self.ordp

    cdef val_unit_c(self):
        cdef Integer val
        cdef pAdicCappedRelativeElement unit
        if mpz_sgn(self.unit) == -1:
            unit = self._new_c()
            unit._set_inexact_zero(0)
            return (infinity, unit)
        self._normalize()
        unit = self._new_c()
        mpz_set(unit.unit, self.unit)
        unit._set_prec(self.relprec)
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
        return hash(self.lift_c())

    def _teichmuller_set(self):
        cdef mpz_t tmp
        self._normalize()
        if mpz_sgn(self.unit) < 0 or self.ordp > 0:
            self._set_exact_zero()
        elif self.ordp < 0:
            raise ValueError, "cannot set negative valuation element to Teichmuller representative."
        elif self.relprec == 0:
            raise ValueError, "not enough precision"
        else:
            mpz_init_set(tmp, self.prime_pow.pow_mpz_t_tmp(self.relprec)[0])
            self.teichmuller_set_c(self.unit, tmp)
            mpz_clear(tmp)

def unpickle_pcre_v1(R, unit, ordp, relprec):
    cdef pAdicCappedRelativeElement ans
    ans = PY_NEW(pAdicCappedRelativeElement)
    ans._parent = R
    ans.prime_pow = <PowComputer_class>R.prime_pow
    mpz_init_set(ans.unit, (<Integer>unit).value)
    ans.ordp = ordp
    ans.relprec = relprec
    ans._normalized = 0
    return ans
