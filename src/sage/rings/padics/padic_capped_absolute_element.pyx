"""
Elements of p-Adic Rings with Absolute Precision Cap

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

cimport sage.rings.rational
cimport sage.rings.padics.pow_computer
from sage.rings.rational cimport Rational
from sage.rings.padics.pow_computer cimport PowComputer_class

import sage.rings.padics.padic_generic_element
import sage.rings.padics.padic_lazy_element
import sage.rings.integer_mod
import sage.libs.pari.gen
import sage.rings.integer
import sage.rings.rational

from sage.rings.infinity import infinity
from sage.rings.integer_mod import Mod
from sage.rings.padics.precision_error import PrecisionError

pAdicLazyElement = sage.rings.padics.padic_lazy_element.pAdicLazyElement
pAdicGenericElement = sage.rings.padics.padic_generic_element.pAdicGenericElement
pari = sage.libs.pari.gen.pari
pari_gen = sage.libs.pari.gen.gen
PariError = sage.libs.pari.gen.PariError

cdef class pAdicCappedAbsoluteElement(pAdicBaseGenericElement):
    def __init__(pAdicCappedAbsoluteElement self, parent, x, absprec=infinity, relprec = infinity, empty=False):
        mpz_init(self.value)
        pAdicGenericElement.__init__(self,parent)
        if empty:
            return
        if absprec > parent.precision_cap():
            absprec = parent.precision_cap()
        elif not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if not relprec is infinity and not PY_TYPE_CHECK(relprec, Integer):
            relprec = Integer(relprec)
        if isinstance(x, pAdicGenericElement):
            if x.valuation() < 0:
                raise ValueError, "element valuation cannot be negative."
            if parent.prime() != x.parent().prime():
                raise TypeError, "Cannot coerce between p-adic parents with different primes"
        if isinstance(x, pAdicLazyElement):
            # We may be doing unnecessary precision increases here...
            if relprec is infinity:
                try:
                    x.set_precision_absolute(absprec)
                except PrecisionError:
                    pass
            else:
                try:
                    x.set_precision_relative(relprec)
                except PrecisionError:
                    try:
                        x.set_precision_absolute(absprec)
                    except PrecisionError:
                        pass
        if isinstance(x, pAdicBaseGenericElement):
            if relprec is infinity:
                absprec = min(x.precision_absolute(), absprec)
            else:
                absprec = min(x.precision_absolute(), x._min_valuation() + relprec, parent.precision_cap())
            self.set_from_Integers(x._integer_(), absprec)
            return

        if isinstance(x, pari_gen):
            if x.type() == "t_PADIC":
                absprec = min(Integer(x.padicprec(parent.prime())), absprec)
                x = x.lift()
            if x.type() == "t_INT":
                x = Integer(x)
            elif x.type() == "t_FRAC":
                x = Rational(x)
            else:
                raise TypeError, "unsupported coercion from pari: only p-adics, integers and rationals allowed"

        cdef mpz_t modulus
        cdef Integer tmp
        cdef unsigned long k
        if sage.rings.integer_mod.is_IntegerMod(x):
            mpz_init_set(modulus, (<Integer>x.modulus()).value)
            k = mpz_remove(modulus, modulus, self.prime_pow.prime.value)
            if mpz_cmp_ui(modulus, 1) == 0:
                if mpz_cmp_ui((<Integer>absprec).value, k) > 0:
                    absprec = PY_NEW(Integer)
                    mpz_set_ui((<Integer>absprec).value, k)
                x = x.lift()
                mpz_clear(modulus)
            else:
                mpz_clear(modulus)
                raise TypeError, "cannot coerce from the given integer mod ring (not a power of the same prime)"

        #Now use the code below to convert from integer or rational, so don't make the next line elif

        if isinstance(x, (int, long)):
            x = Integer(x)
        if isinstance(x, Integer):
            val = x.valuation(parent.prime())
        elif isinstance(x, Rational):
            val = x.valuation(parent.prime())
            if val < 0:
                raise ValueError, "p divides the denominator"
            # todo: make the following understandable and speed it up a bit by just using mpz functions
            absprec = min(absprec, val + relprec)
            self.set_precs(mpz_get_ui((<Integer>absprec).value))
            tmp = PY_NEW(Integer)
            mpz_set(tmp.value, self.prime_pow.modulus)
            tmp = <Integer> x % tmp
            self.set_value_from_mpz(tmp.value)
            return
        if isinstance(x, Integer):
            self.set_from_Integers(x, min(absprec, val + relprec))
        else:
            raise TypeError, "unable to create p-adic element from %s of type %s"%(x, type(x))

    def __dealloc__(self):
        mpz_clear(self.value)

    def __reduce__(self):
        """
        sage: a = ZpCA(5)(-3)
        sage: type(a)
        <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
        sage: loads(dumps(a)) == a
        True
        """
        return make_pAdicCappedAbsoluteElement, (self.parent(), self.lift(), self.absprec)

    cdef void set_precs(pAdicCappedAbsoluteElement self, unsigned long absprec):
        """
        Sets self.absprec
        """
        self.absprec = absprec

    cdef void set_value_from_mpz(pAdicCappedAbsoluteElement self, mpz_t value):
        """
        Assuming that self.absprec is set, sets self.value
        """
        if mpz_sgn(value) == -1 or mpz_cmp(value, self.prime_pow.dense_list[self.absprec]) >= 0:
            mpz_mod(self.value, value, self.prime_pow.dense_list[self.absprec])
        else:
            mpz_set(self.value, value)

    cdef void set_from_Integers(pAdicCappedAbsoluteElement self, Integer value, Integer absprec):
        """
        Set self.value and self.absprec.
        """
        self.set_precs(mpz_get_ui((<Integer>absprec).value))
        self.set_value_from_mpz(value.value)

    cdef pAdicCappedAbsoluteElement _new_c(self):
        cdef pAdicCappedAbsoluteElement x
        x = PY_NEW(pAdicCappedAbsoluteElement)
        x._parent = self._parent
        x.prime_pow = self.prime_pow
        mpz_init(x.value)
        return x

    def __richcmp__(left, right, op):
        return (<Element>left)._richcmp(right, op)

    def __invert__(self):
        return self._invert_c_impl()

    cdef RingElement _invert_c_impl(self):
        return self.parent().fraction_field()(self).__invert__()

    cdef ModuleElement _neg_c_impl(self):
        """
        Returns -x.

        INPUT:
        x -- a p-adic capped absolute element
        OUTPUT:
        -x

        EXAMPLES:
        sage: R = Zp(5, prec=10, type='capped-abs')
        sage: a = R(1)
        sage: -a
        4 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)
        """
        if mpz_sgn(self.value) == 0:
            return self
        cdef pAdicCappedAbsoluteElement ans
        ans = self._new_c()
        ans.absprec = self.absprec
        mpz_sub(ans.value, self.prime_pow.dense_list[self.absprec], self.value)
        return ans

    def __pow__(pAdicCappedAbsoluteElement self, right, dummy):
        cdef Integer new, val, absprec
        cdef mpz_t tmp
        new = Integer(right) #Need to make sure that this works for p-adic exponents
        val = self.valuation_c()
        if (val > 0) and isinstance(right, pAdicBaseGenericElement):
            raise ValueError, "Can only have p-adic exponent if base is a unit"
        if (new < 0):
            return (~self).__pow__(-new)
        cdef pAdicCappedAbsoluteElement ans
        cdef Integer preccap
        preccap = <Integer>self.parent().precision_cap()
        ans = self._new_c()
        if val > 0:
            if new * val >= preccap:
                ans.set_from_Integers(Integer(0), preccap)
            else:
                absprec = <Integer> min(self.precision_relative() + new * val, preccap)
                ans.set_precs(mpz_get_ui(absprec.value))
                _sig_on
                mpz_powm(ans.value, self.value, new.value, self.prime_pow.dense_list[ans.absprec])
                _sig_off
        else:
            ans.absprec = self.absprec
            _sig_on
            mpz_powm(ans.value, self.value, new.value, self.prime_pow.dense_list[self.absprec])
            _sig_off
        return ans

    cdef ModuleElement _add_c_impl(self, ModuleElement _right):
        cdef pAdicCappedAbsoluteElement ans
        cdef pAdicCappedAbsoluteElement right = <pAdicCappedAbsoluteElement> _right
        ans = self._new_c()
        if self.absprec < right.absprec:
            ans.absprec = self.absprec
        else:
            ans.absprec = right.absprec
        mpz_add(ans.value, self.value, right.value)
        if mpz_cmp(ans.value, self.prime_pow.dense_list[ans.absprec]) >= 0:
            mpz_mod(ans.value, ans.value, self.prime_pow.dense_list[ans.absprec])
        return ans

    cdef RingElement _div_c_impl(self, RingElement right):
        return self * (~right)

    def __lshift__(pAdicCappedAbsoluteElement self, shift):
        """
        EXAMPLES:
        We create a capped relative field:
            sage: R = Zp(5, 20, 'capped-rel'); a = R(1000); a
            3*5^3 + 5^4 + O(5^23)

        Shifting to the right is the same as dividing by a power of
        the uniformizer $p$ of the $p$-adic ring.
            sage: a >> 1
            3*5^2 + 5^3 + O(5^22)

        Shifting to the left is the same as multiplying by a power of $p$:
            sage: a << 2
            3*5^5 + 5^6 + O(5^25)
            sage: a*5^2
            3*5^5 + 5^6 + O(5^25)

        Shifting by a negative integer to the left is the same as right shifting
        by the absolute value:
            sage: a << -3
            3 + 5 + O(5^20)
            sage: a >> 3
            3 + 5 + O(5^20)
        """
        cdef pAdicCappedAbsoluteElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            if mpz_sgn((<Integer>shift).value) == 1:
                ans.set_precs(self.prime_pow._cache_limit)
            else:
                ans.set_precs(0)
            return ans
        return self._lshift_c(mpz_get_si((<Integer>shift).value))

    cdef pAdicCappedAbsoluteElement _lshift_c(pAdicCappedAbsoluteElement self, long shift):
        cdef unsigned long prec_cap, ansprec
        cdef pAdicCappedAbsoluteElement ans
        cdef PowComputer_class powerer
        if shift < 0:
            return self._rshift_c(-shift)
        prec_cap = self.prime_pow._cache_limit
        if shift >= prec_cap: # if we ever start caching valuation, this check should change
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            ans.set_precs(prec_cap)
            return ans
        elif shift > 0:
            ans = self._new_c()
            ansprec = shift + self.absprec
            if ansprec > prec_cap:
                ansprec = prec_cap
            ans.set_precs(ansprec)
            mpz_mul(ans.value, self.value, self.prime_pow.dense_list[shift])
            if mpz_cmp(ans.value, ans.prime_pow.modulus) >= 0:
                mpz_mod(ans.value, ans.value, ans.prime_pow.modulus)
            return ans
        else:
            return self

    def __rshift__(pAdicCappedAbsoluteElement self, shift):
        """
        EXAMPLES:
            sage: R = Zp(997, 7, 'capped-rel'); a = R(123456878908); a
            964*997 + 572*997^2 + 124*997^3 + O(997^8)

        Shifting to the right divides by a power of p, but dropping terms with
        negative valuation:
            sage: a >> 3
            124 + O(997^5)

        Shifting to the left multiplies by that power of p.
            sage: a << 3
            964*997^4 + 572*997^5 + 124*997^6 + O(997^11)
        """
        cdef pAdicCappedAbsoluteElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_slong_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            if mpz_sgn((<Integer>shift).value) == -1:
                ans.set_precs(self.prime_pow._cache_limit)
            else:
                ans.set_precs(0)
            return ans
        return self._rshift_c(mpz_get_si((<Integer>shift).value))

    cdef pAdicCappedAbsoluteElement _rshift_c(pAdicCappedAbsoluteElement self, long shift):
        cdef unsigned long prec_cap, ansprec
        cdef pAdicCappedAbsoluteElement ans
        if shift < 0:
            return self._lshift_c(-shift)
        prec_cap = self.prime_pow._cache_limit
        if shift >= self.absprec:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            ans.set_precs(0)
            return ans
        elif shift > 0:
            ans = self._new_c()
            ansprec = self.absprec - shift
            ans.set_precs(ansprec)
            mpz_fdiv_q(ans.value, self.value, self.prime_pow.dense_list[shift])
            return ans
        else:
            return self

    cdef RingElement _mul_c_impl(self, RingElement _right):
        """
        EXAMPLES:
        sage: R = ZpCA(5)
        sage: a = R(20,5); b = R(75, 4); a * b
        2*5^3 + 2*5^4 + O(5^5)
        """
        cdef pAdicCappedAbsoluteElement ans
        cdef pAdicCappedAbsoluteElement right = <pAdicCappedAbsoluteElement> _right
        cdef unsigned long sval, rval, prec1, prec2, prec_cap
        ans = self._new_c()
        prec_cap = self.prime_pow._cache_limit
        sval = self.valuation_c()
        rval = right.valuation_c()
        # Need to think about the possibility that the following overflows (too big for long).
        prec1 = sval + right.absprec
        prec2 = rval + self.absprec
        if prec1 < prec2:
            if prec_cap < prec1:
                ans.set_precs(prec_cap)
            else:
                ans.set_precs(prec1)
        else:
            if prec_cap < prec2:
                ans.set_precs(prec_cap)
            else:
                ans.set_precs(prec2)
        mpz_mul(ans.value, self.value, right.value)
        mpz_mod(ans.value, ans.value, self.prime_pow.dense_list[ans.absprec])
        return ans

    cdef ModuleElement _sub_c_impl(self, ModuleElement _right):
        cdef pAdicCappedAbsoluteElement ans
        cdef pAdicCappedAbsoluteElement right = <pAdicCappedAbsoluteElement> _right
        ans = self._new_c()
        if self.absprec < right.absprec:
            ans.absprec = self.absprec
        else:
            ans.absprec = right.absprec
        mpz_sub(ans.value, self.value, right.value)
        if mpz_sgn(ans.value) == -1 or mpz_cmp(ans.value, self.prime_pow.dense_list[ans.absprec]) >= 0:
            mpz_mod(ans.value, ans.value, self.prime_pow.dense_list[ans.absprec])
        return ans

    def add_bigoh(pAdicCappedAbsoluteElement self, absprec):
        """
        Returns a new element with absolute precision decreased to prec

        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            element -- self with precision set to the minimum of self's precision and prec

        EXAMPLES:
            sage: R = Zp(7,4,'capped-abs','series'); a = R(8); a.add_bigoh(1)
            1 + O(7)
        """
        cdef pAdicCappedAbsoluteElement ans
        cdef unsigned long newprec
        cdef Integer _absprec
        if not PY_TYPE_CHECK(absprec, Integer):
            _absprec = Integer(absprec)
        else:
            _absprec = <Integer>absprec
        if mpz_fits_ulong_p(_absprec.value) == 0:
            if mpz_sgn((<Integer>absprec).value) == -1:
                ans = self._new_c()
                ans.set_precs(0)
                mpz_set_ui(ans.value, 0)
                return ans
            else:
                return self
        newprec = mpz_get_ui(_absprec.value)
        if newprec >= self.absprec:
            return self
        elif newprec <= 0:
            ans = self._new_c()
            ans.set_precs(0)
            mpz_set_ui(ans.value, 0)
            return ans
        else:
            ans = self._new_c()
            ans.set_precs(newprec)
            mpz_set(ans.value, self.value)
            if mpz_cmp(ans.value, self.prime_pow.dense_list[ans.absprec]) >= 0:
                mpz_mod(ans.value, ans.value, self.prime_pow.dense_list[ans.absprec])
            return ans

    def copy(pAdicCappedAbsoluteElement self):
        cdef pAdicCappedAbsoluteElement ans
        ans = self._new_c()
        ans.absprec = self.absprec
        mpz_set(ans.value, self.value)
        return ans

    def exp_artin_hasse(self):
        raise NotImplementedError

    def gamma(self):
        raise NotImplementedError

    def is_zero(self, absprec = None):
        r"""
        Returns whether self is zero modulo $p^{\mbox{absprec}}$.

        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            boolean -- whether self is zero

        """
        if absprec is None:
            return mpz_sgn(self.value) == 0
        cdef Integer _absprec
        if not PY_TYPE_CHECK(absprec, Integer):
            _absprec = Integer(absprec)
        else:
            _absprec = <Integer>absprec
        if mpz_cmp_ui(_absprec.value, self.absprec) > 0:
            raise PrecisionError, "Not enough precision to determine if element is zero"
        cdef mpz_t tmp
        mpz_init(tmp)
        cdef unsigned long aprec
        aprec = mpz_get_ui(_absprec.value)
        mpz_mod(tmp, self.value, self.prime_pow.dense_list[aprec])
        if mpz_sgn(tmp) == 0:
            mpz_clear(tmp)
            return True
        else:
            mpz_clear(tmp)
            return False

    def is_equal_to(self, right, absprec = None): #assumes they have the same parent
        r"""
        Returns whether self is equal to right modulo $p^{\mbox{absprec}}$.

        INPUT:
            self -- a p-adic element
            right -- a p-addic element with the same parent
            absprec -- an integer
        OUTPUT:
            boolean -- whether self is equal to right

        """
        cdef unsigned long aprec
        if absprec is None:
            if self.absprec < right.absprec:
                aprec = self.absprec
            else:
                aprec = right.absprec
        else:
            if not PY_TYPE_CHECK(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_slong_p((<Integer>absprec).value) == 0:
                if mpz_sgn((<Integer>absprec).value) == 1:
                    raise PrecisionError, "Elements are not known to specified precision"
                else:
                    return True
            aprec = mpz_get_ui((<Integer>absprec).value)
            if aprec > self.absprec or aprec > right.absprec:
                raise PrecisionError, "Elements are not known to specified precision"
        cdef mpz_t tmp1, tmp2
        mpz_init(tmp1)
        mpz_init(tmp2)
        mpz_mod(tmp1, self.value, self.prime_pow.dense_list[aprec])
        mpz_mod(tmp2, (<pAdicCappedAbsoluteElement>right).value, self.prime_pow.dense_list[aprec])
        if mpz_cmp(tmp1, tmp2) == 0:
            mpz_clear(tmp1)
            mpz_clear(tmp2)
            return True
        else:
            mpz_clear(tmp1)
            mpz_clear(tmp2)
            return False

    def lift(self):
        """
        Returns an integer congruent to this p-adic element modulo p^self.absprec().
        """
        return self.lift_c()

    cdef Integer lift_c(pAdicCappedAbsoluteElement self):
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set(ans.value, self.value)
        return ans

    def lift_to_precision(self, absprec):
        cdef pAdicCappedAbsoluteElement ans
        cdef unsigned long prec_cap
        cdef unsigned long dest_prec
        prec_cap = self.prime_pow._cache_limit
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_fits_ulong_p((<Integer>absprec).value) == 0:
            ans = self._new_c()
            ans.set_precs(prec_cap)
            mpz_set(ans.value, self.value)
        else:
            dest_prec = mpz_get_ui((<Integer>absprec).value)
            if prec_cap < dest_prec:
                dest_prec = prec_cap
            if dest_prec <= self.absprec:
                return self
            else:
                ans = self._new_c()
                ans.set_precs(dest_prec)
                mpz_set(ans.value, self.value)
        return ans

    def list(pAdicCappedAbsoluteElement self, lift_mode = 'simple'):
        """
        Returns a list of coefficients of p starting with $p^0$
        INPUT:
            self -- a p-adic element
            lift_mode -- 'simple', 'smallest' or 'teichmuller' (default 'simple')
        OUTPUT:
            list -- the list of coeficients of self

        NOTES:
        Returns a list [a_0, a_1, \ldots, a_n] so that each a_i is an integer ('simple' and 'smallest')
        or p-adic ('teichmuller') and \sum_{i = 0}^n a_i * p^i = self, modulo the precision of self.
        If lift_mode = 'simple', 0 <= a_i < p.
        If lift_mode = 'smallest', -p/2 < a_i <= p/2.
        If lift_mode = 'teichmuller', a_i^p = a_i (mod p^(self.precision_absolute() - i))

        EXAMPLES:
        sage: R = Zp(7,4,'capped-abs'); a = R(2*7+7**2); a.list()
        [0, 2, 1]
        """
        if lift_mode == 'teichmuller':
            return self.teichmuller_list()
        else:
            return self.base_p_list(self.value, lift_mode)

    cdef object teichmuller_list(pAdicCappedAbsoluteElement self):
        # May eventually want to add a dict to store teichmuller lifts already seen, if p small enough
        cdef unsigned long curpower
        cdef mpz_t tmp
        cdef pAdicCappedAbsoluteElement list_elt
        ans = PyList_New(0)
        curpower = self.absprec
        mpz_init_set(tmp, self.value)
        while mpz_sgn(tmp) != 0:
            curpower -= 1
            list_elt = self._new_c()
            list_elt.set_precs(curpower)
            mpz_mod(list_elt.value, tmp, self.prime_pow.prime.value)
            sage.rings.padics.padic_generic_element.teichmuller_set_c(list_elt.value, self.prime_pow.prime.value, self.prime_pow.dense_list[curpower])
            mpz_sub(tmp, tmp, list_elt.value)
            mpz_divexact(tmp, tmp, self.prime_pow.prime.value)
            mpz_mod(tmp, tmp, self.prime_pow.dense_list[curpower])
            PyList_Append(ans, list_elt)
        mpz_clear(tmp)
        return ans

    def _teichmuller_set(self, Integer n, Integer absprec):
        mpz_set(self.value, n.value)
        cdef unsigned long aprec
        if mpz_fits_ulong_p(absprec.value) == 0:
            aprec = self.prime_pow._cache_limit
        if mpz_sgn(absprec.value) != 1:
            raise ValueError, "can only compute to positive precision"
        aprec = mpz_get_ui(absprec.value)
        if aprec > self.prime_pow._cache_limit:
            aprec = self.prime_pow._cache_limit
        self.set_precs(aprec)
        sage.rings.padics.padic_generic_element.teichmuller_set_c(self.value, self.prime_pow.prime.value, self.prime_pow.dense_list[self.absprec])

    def log_artin_hasse(self):
        raise NotImplementedError

    def multiplicative_order(self):
        r"""
        Returns the minimum possible multiplicative order of self.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the multiplicative order of self.  This is the minimum multiplicative order of all elements of Z_p lifting self to infinite precision.
        """
        cdef mpz_t tmp
        cdef Integer ans
        if mpz_divisible_p(self.value, self.prime_pow.prime.value):
            return infinity
        if mpz_cmp_ui(self.value, 1) == 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 1)
            return ans
        mpz_init(tmp)
        mpz_sub_ui(tmp, self.prime_pow.dense_list[self.absprec], 1)
        if mpz_cmp(self.value, tmp) == 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 2)
            return ans
        # check if self is an approximation to a teichmuller lift:
        mpz_powm(tmp, self.value, self.prime_pow.prime.value, self.prime_pow.dense_list[self.absprec])
        if mpz_cmp(tmp, self.value) == 0:
            mpz_clear(tmp)
            return self.residue(1).multiplicative_order()
        else:
            mpz_clear(tmp)
            return infinity

    def padded_list(self, n, list_mode = 'simple'):
        """
        Returns a list of coefficients of p starting with $p^0$ up to $p^n$ exclusive (padded with zeros if needed)
        INPUT:
            self -- a p-adic element
            n - an integer
        OUTPUT:
            list -- the list of coeficients of self
        EXAMPLES:
            sage: R = Zp(7,4,'capped-abs'); a = R(2*7+7**2); a.padded_list(5)
                [0, 2, 1, 0, 0]

        NOTE:
            this differs from the padded_list method of padic_field_element
            the slice operators throw an error if asked for a slice above the precision, while this function works
        """
        if list_mode == 'simple' or list_mode == 'smallest':
            zero = Integer(0)
        else:
            zero = self.parent()(0, 0)
        L = self.list()
        return L[:n] + [zero] * (n - len(L))

    def precision_absolute(self):
        """
        Returns the absolute precision of self
         INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the absolute precision of self
        EXAMPLES:
            sage: R = Zp(7,4,'capped-abs'); a = R(7); a.precision_absolute()
                4
       """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.absprec)
        return ans

    def precision_relative(self):
        """
        Returns the relative precision of self
         INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the relative precision of self
        EXAMPLES:
            sage: R = Zp(7,4,'capped-abs'); a = R(7); a.precision_relative()
                3
       """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.absprec - self.valuation_c())
        return ans

    def residue(self, absprec):
        r"""
        Reduces self modulo $p^absprec$

        INPUT:
            self -- a p-adic element
            absprec - an integer

        OUTPUT:
            element of Z/(p^prec Z) -- self reduced mod p^prec

        EXAMPLES:
            sage: R = Zp(7,4,'capped-abs'); a = R(8); a.residue(1)
            1
        """
        if absprec > self.precision_absolute():
            raise PrecisionError, "Not enough precision known in order to compute residue."
        elif absprec < 0:
            raise ValueError, "cannot reduce modulo a negative power of p"
        cdef Integer selfvalue
        selfvalue = PY_NEW(Integer)
        mpz_set(selfvalue.value, self.value)
        return Mod(selfvalue, self.parent().prime_pow(absprec))


    #def square_root(self):
    #    r"""
    #    Returns the square root of this p-adic number

    #    INPUT:
    #        self -- a p-adic element
    #    OUTPUT:
    #        p-adic element -- the square root of this p-adic number
    #    EXAMPLES:
    #        sage: R = Zp(3,20,'capped-abs')
    #        sage: R(0).square_root()
    #            O(3^10)
    #        sage: R(1).square_root()
    #            1 + O(3^20)
    #        sage: R(4).square_root() == R(-2)
    #            True
    #        sage: R(9).square_root()
    #            3 + O(3^19)
    #        sage: R2 = Zp(2,20,'capped-abs')
    #        sage: R2(0).square_root()
    #            O(2^10)
    #        sage: R2(1).square_root()
    #            1 + O(2^19)
    #        sage: R2(4).square_root()
    #            2 + O(2^18)
    #        sage: R2(9).square_root() == R2(3) or R2(9).square_root() == R2(-3)
    #            True
    #        sage: R2(17).square_root()
    #            1 + 2^3 + 2^5 + 2^6 + 2^7 + 2^9 + 2^10 + 2^13 + 2^16 + 2^17 + O(2^19)
    #        sage: R3 = Zp(5,20,'capped-abs')
    #        sage: R3(0).square_root()
    #            O(5^10)
    #        sage: R3(1).square_root()
    #            1 + O(5^20)
    #        sage: R3(-1).square_root() == R3.teichmuller(2) or R3(-1).square_root() == R3.teichmuller(3)
    #            True
    #    """
    #    #todo: make more efficient
    #    try:
    #        # use pari
    #        return self.parent()(pari(self).sqrt())
    #    except PariError:
    #        raise ValueError, "element is not a square" # should eventually change to return an element of an extension field

    def unit_part(self):
        r"""
        Returns the unit part of self.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            p-adic element -- the unit part of self
        EXAMPLES:
            sage: R = Zp(17,4,'capped-abs', 'val-unit')
            sage: a = R(18*17)
            sage: a.unit_part()
                18 + O(17^3)
            sage: type(a)
                <type 'sage.rings.padics.padic_capped_absolute_element.pAdicCappedAbsoluteElement'>
        """
        return self.unit_part_c()

    cdef pAdicCappedAbsoluteElement unit_part_c(pAdicCappedAbsoluteElement self):
        cdef pAdicCappedAbsoluteElement ans
        cdef unsigned long v
        if mpz_sgn(self.value) == 0:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            ans.set_precs(0)
            return ans
        elif mpz_divisible_p(self.value, self.prime_pow.prime.value):
            ans = self._new_c()
            _sig_on
            v = mpz_remove(ans.value, self.value, self.prime_pow.prime.value)
            _sig_off
            ans.set_precs(self.absprec - v)
            return ans
        else:
            return self

    def valuation(self):
        """
        Returns the valuation of self, ie the largest power of p dividing self.

        EXAMPLES:
        sage: R = ZpCA(5)
        sage: a = R(5^5*1827); a.valuation()
        5
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_ui(ans.value, self.valuation_c())
        return ans

    cdef unsigned long valuation_c(self):
        if mpz_sgn(self.value) == 0:
            return self.absprec
        cdef mpz_t tmp
        cdef unsigned long ans
        mpz_init(tmp)
        ans = mpz_remove(tmp, self.value, self.prime_pow.prime.value)
        mpz_clear(tmp)
        return ans

    def val_unit(self):
        """
        Returns a 2-tuple, the first element set to the valuation of self, and the second to the unit part of self.

        If self == 0, then the unit part is O(p^0).

        EXAMPLES:
        sage: R = ZpCA(5)
        sage: a = R(75, 6); b = a - a
        sage: a.val_unit()
        (2, 3 + O(5^4))
        sage: b.val_unit()
        (6, O(5^0))
        """
        return self.val_unit_c()

    cdef val_unit_c(self):
        cdef pAdicCappedAbsoluteElement unit
        cdef Integer val
        cdef unsigned long v
        val = PY_NEW(Integer)
        if mpz_sgn(self.value) == 0:
            unit = self._new_c()
            mpz_set_ui(unit.value, 0)
            unit.set_precs(0)
            mpz_set_ui(val.value, self.absprec)
            return (val, unit)
        elif mpz_divisible_p(self.value, self.prime_pow.prime.value):
            unit = self._new_c()
            v = mpz_remove(unit.value, self.value, self.prime_pow.prime.value)
            unit.set_precs(self.absprec - v)
            mpz_set_ui(val.value, v)
            return (val, unit)
        else:
            mpz_set_ui(val.value, 0)
            return (val, self)

    def __hash__(self):
        return hash(self.lift_c())

def make_pAdicCappedAbsoluteElement(parent, x, absprec):
    return parent(x, absprec=absprec)
