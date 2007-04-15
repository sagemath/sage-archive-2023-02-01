"""
Elements of p-Adic Rings with Fixed Modulus

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

cimport sage.rings.integer
cimport sage.rings.rational
cimport sage.rings.padics.pow_computer
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.padics.pow_computer cimport PowComputer_class

import sage.rings.padics.padic_ring_generic_element
import sage.rings.padics.padic_field_generic_element
import sage.rings.padics.padic_lazy_element
import sage.rings.integer_mod
import sage.libs.pari.gen
import sage.rings.integer
import sage.rings.rational

from sage.rings.infinity import infinity
from sage.rings.integer_mod import Mod
from sage.rings.padics.precision_error import PrecisionError

#pAdicLazyElement = sage.rings.padics.padic_lazy_element.pAdicLazyElement
#pAdicGenericElement = sage.rings.padics.padic_generic_element.pAdicGenericElement
pari = sage.libs.pari.gen.pari
pari_gen = sage.libs.pari.gen.gen
PariError = sage.libs.pari.gen.PariError

cdef class pAdicRingFixedModElement(pAdicBaseGenericElement):
    def __init__(pAdicRingFixedModElement self, parent, x, absprec = None, relprec = None, empty = False):
        r"""
        INPUT:
            parent -- a pAdicRingFixedMod object.

        Types currently supported:
            Integers
            Rationals -- denominator must be relatively prime to p
            FixedMod p-adics

        Types that should be supported:
            Finite precision p-adics
            Lazy p-adics
            Elements of local extensions of THIS p-adic ring that actually lie in Zp
            Elements of IntegerModRing(p^k) for k less than or equal to the modulus

        EXAMPLES:
            sage: R = Zp(5, 20, 'fixed-mod', 'terse')

        Construct from integers:
            sage: R(3)
            3 + O(5^20)
            sage: R(75)
            75 + O(5^20)
            sage: R(0)
            0 + O(5^20)

            sage: R(-1)
            95367431640624 + O(5^20)
            sage: R(-5)
            95367431640620 + O(5^20)

        Construct from rationals:
            sage: R(1/2)
            47683715820313 + O(5^20)
            sage: R(-7875/874)
            9493096742250 + O(5^20)
            sage: R(15/425)
            Traceback (most recent call last):
            ...
            ValueError: p divides the denominator

        # todo: the above error message does not agree with the error message
        # in the corresponding capped-relative constructor

        Construct from IntegerMod:
            sage: R(Integers(125)(3))
            3 + O(5^20)
            sage: R(Integers(5)(3))
            3 + O(5^20)
            sage: R(Integers(5^30)(3))
            3 + O(5^20)
            sage: R(Integers(5^30)(1+5^23))
            1 + O(5^20)
            sage: R(Integers(49)(3))
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce from the given integer mod ring (not a power of the same prime)

            sage: R(Integers(48)(3))
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce from the given integer mod ring (not a power of the same prime)

        Some other conversions:
            sage: R(R(5))
            5 + O(5^20)

        # todo: doctests for converting from other types of p-adic rings

        """
        mpz_init(self.value)

        # what I want to do is kinda sketchy: we're direcly copying the mpz_t's, which point
        # to the actual data.  I think this is okay: parent's
        # prime_pow instance will cache p^parent.precision_cap(),
        # holding a reference and thus preventing Python from garbage
        # collecting it, and parent's _p field will store the prime.

        #self.modulus = (<Integer>parent.prime_pow(parent.precision_cap())).value
        #self.p = (<Integer>parent.prime()).value

        #Pyrex won't let me, so I'm doing it the less sketchy way.
        mpz_init(self.modulus)
        mpz_set(self.modulus, (<Integer>parent.prime_pow(parent.precision_cap())).value)
        mpz_init(self.p)
        mpz_set(self.p, (<Integer>parent.prime()).value)

        CommutativeRingElement.__init__(self,parent)
        if empty:
            return

        #if PY_TYPE_CHECK(x, pAdicGenericElement) and x.valuation() < 0:
        #    raise ValueError, "element has negative valuation"
        cdef Integer tmp
        #if PY_TYPE_CHECK(x, pAdicLazyElement):
        #    try:
        #        x.set_precision_absolute(absprec)
        #    except PrecisionError:
        #        pass
        #    if mpz_cmp((<pAdicLazyElement>x).value, self.modulus) >= 0:
        #        mpz_mod(self.value, (<pAdicLazyElement>x).value, self.modulus)
        #    else:
        #        mpz_set(self.value, (<pAdicLazyElement>x).value)
        #    return
        if PY_TYPE_CHECK(x, pAdicGenericElement):
            if x.valuation() < 0:
                raise ValueError, "element has negative valuation"
            if parent.prime() != x.parent().prime():
                raise TypeError, "Cannot coerce between p-adic parents with different primes."
        if PY_TYPE_CHECK(x, pAdicBaseGenericElement):
            tmp = <Integer> x.lift()
            if mpz_cmp(tmp.value, self.modulus) >= 0:
                mpz_mod(self.value, tmp.value, self.modulus)
            else:
                mpz_set(self.value, tmp.value)
            return

        if isinstance(x, pari_gen):
            if x.type() == "t_PADIC":
                x = x.lift()
            if x.type() == 't_INT':
                x = Integer(x)
            elif x.type() == 't_FRAC':
                x = Rational(x)
            else:
                raise TypeError, "unsupported coercion from pari: only p-adics, integers and rationals allowed"

        if sage.rings.integer_mod.is_IntegerMod(x):
            if (<Integer>x.modulus())._is_power_of(<Integer>parent.prime()):
                x = x.lift()
            else:
                raise TypeError, "cannot coerce from the given integer mod ring (not a power of the same prime)"

        #if sage.rings.finite_field_element.is_FiniteFieldElement(x):
        #    if x.parent().order() != parent.prime():
        #        raise TypeError, "can only create p-adic element out of finite field when order of field is p"
        #    x = x.lift()

        #Now use the code below to convert from integer or rational, so don't make the next line elif

        if PY_TYPE_CHECK(x, Integer):
            self.set_from_mpz((<Integer>x).value)
        elif isinstance(x, Rational):
            if parent.prime().divides(x.denominator()):
                raise ValueError, "p divides the denominator"
            else:
                tmp = <Integer> x % parent.prime_pow(parent.precision_cap())
                self.set_from_mpz(tmp.value)
        elif isinstance(x, (int, long)):
            tmp = <Integer> Integer(x)
            self.set_from_mpz(tmp.value)
        else:
            raise TypeError, "unable to create p-adic element"

    def __dealloc__(self):
        # If we switch to using the pointer copying method for setting self.modulus and self.p, the corresponding mpz_clear's here should be removed
        mpz_clear(self.modulus)
        mpz_clear(self.p)
        mpz_clear(self.value)

    cdef void set_from_mpz(pAdicRingFixedModElement self, mpz_t value):
        if mpz_sgn(value) == -1 or mpz_cmp(value, self.modulus) >= 0:
            mpz_mod(self.value, value, self.modulus)
        else:
            mpz_set(self.value, value)

    cdef pAdicRingFixedModElement _new_c(self):
        cdef pAdicRingFixedModElement x
        x = PY_NEW(pAdicRingFixedModElement)
        x._parent = self._parent
        mpz_init(x.value)
        # I want to just copy the pointers:
        #x.modulus = self.modulus
        #x.p = self.p
        # ... but Pyrex apparently won't let me.
        mpz_init_set(x.modulus, self.modulus)
        mpz_init_set(x.p, self.p)
        return x

    def __invert__(self):
        r"""
        Returns multiplicative inverse of this element. Its valuation
        must be zero.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: ~R(2)
            4 + 3*7 + 3*7^2 + 3*7^3 + O(7^4)
            sage: ~R(0)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
            sage: ~R(7)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
        """
        return self._invert_c_impl()

    cdef RingElement _invert_c_impl(self):
        cdef int t
        cdef pAdicRingFixedModElement ans
        t = mpz_divisible_p(self.value, self.p)
        if t:
            raise ValueError, "cannot invert non-unit"
        else:
            ans = self._new_c()
            _sig_on
            mpz_invert(ans.value, self.value, self.modulus)
            _sig_off
            return ans

    cdef pAdicRingFixedModElement _lshift_c(pAdicRingFixedModElement self, int shift):
        cdef pAdicRingFixedModElement ans
        cdef int prec_cap
        cdef mpz_t ppow, low_mod
        cdef PowComputer_class powerer
        if shift < 0:
            return self._rshift_c(-shift)
        prec_cap = mpz_get_si((<Integer>self.parent().precision_cap()).value)
        if shift >= prec_cap:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            return ans
        elif shift > 0:
            ans = self._new_c()
            powerer = <PowComputer_class>self.parent().prime_pow
            # We could revise PowComputer_class to provide direct
            # access to the entries stored in dense_list; this would
            # prevent having to call mpz_init and clear
            mpz_init(ppow)
            mpz_init(low_mod)
            powerer.pow_mpz_ui(low_mod, prec_cap - shift)
            if mpz_cmp(low_mod, self.value) <= 0:
                mpz_mod(ans.value, self.value, low_mod)
            else:
                mpz_set(ans.value, self.value)
            powerer.pow_mpz_ui(ppow, shift)
            mpz_mul(ans.value, ans.value, ppow)
            mpz_clear(low_mod)
            mpz_clear(ppow)
            return ans
        else:
            return self

    def __lshift__(pAdicRingFixedModElement self, shift):
        cdef pAdicRingFixedModElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_sint_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            return ans
        return self._lshift_c(mpz_get_si((<Integer>shift).value))

    cdef pAdicRingFixedModElement _rshift_c(pAdicRingFixedModElement self, int shift):
        cdef pAdicRingFixedModElement ans
        cdef int prec_cap
        cdef mpz_t ppow, low_mod
        cdef PowComputer_class powerer
        if shift < 0:
            return self._lshift_c(-shift)
        prec_cap = mpz_get_si((<Integer>self.parent().precision_cap()).value)
        if shift >= prec_cap:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            return ans
        elif shift > 0:
            ans = self._new_c()
            powerer = <PowComputer_class>self.parent().prime_pow
            # We could revise PowComputer_class to provide direct
            # access to the entries stored in dense_list; this would
            # prevent having to call mpz_init and clear
            mpz_init(ppow)
            powerer.pow_mpz_ui(ppow, shift)
            mpz_fdiv_q(ans.value, self.value, ppow)
            mpz_clear(ppow)
            return ans
        else:
            return self

    def __rshift__(pAdicRingFixedModElement self, shift):
        cdef pAdicRingFixedModElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_sint_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            return ans
        return self._rshift_c(mpz_get_si((<Integer>shift).value))

    cdef ModuleElement _neg_c_impl(self):
        r"""
        Returns negative of self.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: -R(7)
            6*7 + 6*7^2 + 6*7^3 + O(7^4)
        """
        if mpz_sgn(self.value) == 0:
            return self
        cdef pAdicRingFixedModElement ans
        ans = self._new_c()
        mpz_sub(ans.value, self.modulus, self.value)
        return ans

    def __pow__(pAdicRingFixedModElement self, right, m): # NOTE: m ignored, always use self.modulus
        if not PY_TYPE_CHECK(right, Integer):
            right = Integer(right) #Need to make sure that this works for p-adic exponents
        cdef pAdicRingFixedModElement ans
        ans = self._new_c()
        _sig_on
        mpz_powm(ans.value, self.value, (<Integer>right).value, self.modulus)
        _sig_off
        return ans

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        r"""
        Returns sum of self and right.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: x = R(1721); x
            6 + 5*7^3 + O(7^4)
            sage: y = R(1373); y
            1 + 4*7^3 + O(7^4)
            sage: x + y
            7 + 2*7^3 + O(7^4)
        """
        cdef pAdicRingFixedModElement ans
        ans = self._new_c()
        mpz_add(ans.value, self.value, (<pAdicRingFixedModElement>right).value)
        if mpz_cmp(ans.value, self.modulus) >= 0:
            mpz_sub(ans.value, ans.value, self.modulus)
        return ans

    cdef RingElement _mul_c_impl(self, RingElement right):
        r"""
        Returns product of self and right.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: R(3) * R(2)
            6 + O(7^4)
            sage: R(1/2) * R(2)
            1 + O(7^4)
        """
        cdef pAdicRingFixedModElement ans
        ans = self._new_c()
        mpz_mul(ans.value, self.value, (<pAdicRingFixedModElement>right).value)
        mpz_fdiv_r(ans.value, ans.value, self.modulus)
        return ans

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        r"""
        Returns difference of self and right.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: x = R(1721); x
            6 + 5*7^3 + O(7^4)
            sage: y = R(1373); y
            1 + 4*7^3 + O(7^4)
            sage: x - y
            5 + 7^3 + O(7^4)
        """
        cdef pAdicRingFixedModElement ans
        ans = self._new_c()
        mpz_sub(ans.value, self.value, (<pAdicRingFixedModElement>right).value)
        if mpz_sgn(ans.value) == -1:
            mpz_add(ans.value, ans.value, self.modulus)
        return ans

    cdef RingElement _div_c_impl(self, RingElement right):
        r"""
        Returns quotient of self and right. The latter must have
        valuation zero.

        EXAMPLES:
            sage: R = Zp(7, 4, 'fixed-mod', 'series')
            sage: R(3) / R(2)
            5 + 3*7 + 3*7^2 + 3*7^3 + O(7^4)
            sage: R(5) / R(0)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
            sage: R(7) / R(49)
            Traceback (most recent call last):
            ...
            ValueError: cannot invert non-unit
        """
        cdef int t
        cdef pAdicRingFixedModElement ans
        t = mpz_divisible_p((<pAdicRingFixedModElement>right).value, self.p)
        if t:
            raise ValueError, "cannot invert non-unit"
        else:
            ans = self._new_c()
            # Can the cast be inside _sig_on?
            _sig_on
            mpz_invert(ans.value, (<pAdicRingFixedModElement>right).value, self.modulus)
            mpz_mul(ans.value, ans.value, self.value)
            mpz_fdiv_r(ans.value, ans.value, self.modulus)
            _sig_off
            return ans

    def add_bigoh(self, prec):
        """
        Returns a new element with precision decreased to prec
        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            element -- self with precision set to the minimum of  self's precision and prec

        EXAMPLES:
            sage: R = Zp(7,4,'fixed-mod','series'); a = R(8); a.add_bigoh(1)
            1 + O(7^4)
        """
        cdef pAdicRingFixedModElement ans
        cdef mpz_t ppow
        cdef PowComputer_class powerer
        if not PY_TYPE_CHECK(prec, Integer):
            prec = Integer(prec)
        if mpz_cmp((<Integer>prec).value, (<Integer>self.parent().precision_cap()).value) >= 0:
            return self
        ans = self._new_c()
        mpz_init(ppow)
        powerer = <PowComputer_class> self.parent().prime_pow
        powerer.pow_mpz_mpz(ppow, (<Integer>prec).value)
        mpz_mod(ans.value, self.value, ppow)
        mpz_clear(ppow)
        return ans

    def copy(self):
        cdef pAdicRingFixedModElement ans
        ans = self._new_c()
        mpz_set(ans.value, self.value)
        return ans

    def exp_artin_hasse(self):
        raise NotImplementedError

    def gamma(self):
        raise NotImplementedError

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
            return bool(mpz_sgn(self.value) == 0)
        if not PY_TYPE_CHECK(prec, Integer):
            prec = Integer(prec)
        if mpz_cmp((<Integer>prec).value, (<Integer>self.parent().precision_cap()).value) >= 0:
            return bool(mpz_sgn(self.value) == 0)
        cdef mpz_t ppow, tmp
        cdef PowComputer_class powerer
        mpz_init(ppow)
        mpz_init(tmp)
        powerer = <PowComputer_class> self.parent().prime_pow
        powerer.pow_mpz_mpz(ppow, (<Integer>prec).value)
        mpz_mod(tmp, self.value, ppow)
        mpz_clear(ppow)
        if mpz_sgn(tmp) == 0:
            mpz_clear(tmp)
            return True
        else:
            mpz_clear(tmp)
            return False

    def is_equal_to(self, right, prec = None): #assumes they have the same parent
        r"""
        Returns whether self is equal to right modulo $p^{\mbox{prec}}$.

        INPUT:
            self -- a p-adic element
            right -- a p-addic element
            prec -- an integer
        OUTPUT:
            boolean -- whether self is equal to right

        """
        if prec is None:
            return bool(mpz_cmp(self.value, (<pAdicRingFixedModElement>right).value) == 0)
        if not PY_TYPE_CHECK(prec, Integer):
            prec = Integer(prec)
        if mpz_cmp((<Integer>prec).value, (<Integer>self.parent().precision_cap()).value) >= 0:
            return bool(mpz_cmp(self.value, (<pAdicRingFixedModElement>right).value) == 0)
        cdef mpz_t ppow, tmp1, tmp2
        cdef PowComputer_class powerer
        mpz_init(ppow)
        mpz_init(tmp1)
        mpz_init(tmp2)
        powerer = <PowComputer_class> self.parent().prime_pow
        powerer.pow_mpz_mpz(ppow, (<Integer>prec).value)
        mpz_mod(tmp1, self.value, ppow)
        mpz_mod(tmp2, (<pAdicRingFixedModElement>right).value, ppow)
        mpz_clear(ppow)
        if mpz_cmp(tmp1, tmp2) == 0:
            mpz_clear(tmp1)
            mpz_clear(tmp2)
            return True
        else:
            mpz_clear(tmp1)
            mpz_clear(tmp2)
            return False

    def lift(self):
        r"""
        Return an integer congruent to self modulo self's precision.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- a integer congruent to self mod $p^{\mbox{prec}}$
        EXAMPLES:
            sage: R = Zp(7,4,'fixed-mod'); a = R(8); a.lift()
            8
            sage: type(a.lift())
            <type 'sage.rings.integer.Integer'>
        """
        return self.lift_c()

    cdef Integer lift_c(pAdicRingFixedModElement self):
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set(ans.value, self.value)
        return ans

    def lift_to_precision(self, absprec):
        return self

    def list(self, lift_mode = None):
        r"""
        Returns a list of coefficients of p starting with $p^0$.

        INPUT:
            self -- a p-adic element
            lift_mode -- 'simple' (default), 'smallest', 'teichmuller'
        OUTPUT:

            list -- a list [a_0, ..., a_k] such that self is congruent
            to a_0 + a_1*p + ... +a_k*p^k.  If lift_mode is 'simple'
            or 'smallest' the a_i will be Integers; if it is 'teichmuller'
            they will be pAdicRingFixedModElements

        EXAMPLES:
            sage: R = Zp(7,4,'fixed-mod'); a = R(2*7+7**2); a.list()
            [0, 2, 1]

        NOTE:
            this differs from the list method of padic_field_element
        """
        if lift_mode is None:
            lift_mode = 'simple'
        elif lift_mode == 'teichmuller':
            return self.teichmuller_list()
        return self.base_p_list(self.value, self.p, lift_mode, self.parent().prime_pow, mpz_get_si((<Integer>self.parent().precision_cap()).value))

    cdef object teichmuller_list(pAdicRingFixedModElement self):
        # May eventually want to add a dict to store teichmuller lifts already seen, if p small enough
        cdef int curpower, preccap
        cdef mpz_t ppow, tmp
        cdef pAdicRingFixedModElement list_elt
        cdef PowComputer_class powerer
        powerer = <PowComputer_class>self.parent().prime_pow
        ans = PyList_New(0)
        preccap = mpz_get_si((<Integer>self.parent().precision_cap()).value)
        curpower = preccap
        mpz_init(ppow)
        mpz_init_set(tmp, self.value)
        while mpz_sgn(tmp) != 0:
            curpower -= 1
            list_elt = self._new_c()
            mpz_mod(list_elt.value, tmp, self.p)
            powerer.pow_mpz_ui(ppow, preccap)
            sage.rings.padics.padic_generic_element.teichmuller_set_c(list_elt.value, self.p, ppow)
            mpz_sub(tmp, tmp, list_elt.value)
            mpz_divexact(tmp, tmp, self.p)
            powerer.pow_mpz_ui(ppow, curpower)
            mpz_mod(tmp, tmp, ppow)
            PyList_Append(ans, list_elt)
        mpz_clear(ppow)
        mpz_clear(tmp)
        return ans

    def _teichmuller_set(self, Integer n, Integer prec):
        cdef mpz_t ppow
        mpz_init(ppow)
        mpz_set(self.value, n.value)
        if mpz_fits_slong_p(prec.value) == 0:
            raise ValueError, "cannot computer teichmuller lift to that high precision"
        if mpz_sgn(prec.value) != 1:
            raise ValueError, "can only compute to positive precision"
        (<PowComputer_class>self.parent().prime_pow).pow_mpz_mpz(ppow, prec.value)
        sage.rings.padics.padic_generic_element.teichmuller_set_c(self.value, self.p, ppow)
        mpz_clear(ppow)

    def log_artin_hasse(self):
        raise NotImplementedError

    def multiplicative_order(self):
        r"""
        Returns the multiplicative order of self, where self is considered to
        be 1 if it is 1 modulo $p^{\mbox{prec}}$.

        INPUT:
            self -- a p-adic element
            prec -- an integer
        OUTPUT:
            integer -- the multiplicative order of self
        """
        cdef mpz_t tmp
        cdef Integer ans
        if mpz_divisible_p(self.value, self.p):
            return infinity
        if mpz_cmp_ui(self.value, 1):
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 1)
            return ans
        mpz_init(tmp)
        mpz_sub_ui(tmp, self.modulus, 1)
        if mpz_cmp(self.value, tmp) == 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 2)
            return ans
        # check if self is an approximation to a teichmuller lift:
        mpz_powm(tmp, self.value, self.p, self.modulus)
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
            sage: R = Zp(7,4,'fixed-mod'); a = R(2*7+7**2); a.padded_list(5)
                [0, 2, 1, 0, 0]

        NOTE:
            this differs from the padded_list method of padic_field_element
            the slice operators throw an error if asked for a slice above the precision, while this function works
        """
        if list_mode == 'simple' or list_mode == 'smallest':
            zero = Integer(0)
        else:
            zero = self.parent()(0)
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
            sage: R = Zp(7,4,'fixed-mod'); a = R(7); a.precision_absolute()
            4
        """
        return self.parent().precision_cap()

    def precision_relative(self):
        r"""
        Returns the relative precision of self
         INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the relative precision of self
        EXAMPLES:
            sage: R = Zp(7,4,'fixed-mod'); a = R(7); a.precision_relative()
            3
            sage: a = R(0); a.precision_relative()
            0
        """
        cdef Integer ans, preccap
        cdef mpz_t val
        mpz_init(val)
        preccap = <Integer> self.parent().precision_cap()
        ans = PY_NEW(Integer)
        mpz_set_si(val, self.valuation_c())
        mpz_sub(ans.value, preccap.value, val)
        mpz_clear(val)
        return ans

    def residue(self, prec):
        r"""
        Reduces this mod $p^prec$

        INPUT:
            self -- a p-adic element
            prec - an integer

        OUTPUT:
            element of Z/(p^prec Z) -- self reduced mod p^prec

        EXAMPLES:
            sage: R = Zp(7,4,'fixed-mod'); a = R(8); a.residue(1)
            1
        """
        cdef Integer selfvalue
        selfvalue = PY_NEW(Integer)
        mpz_set(selfvalue.value, self.value)
        return Mod(selfvalue, self.parent().prime_pow(prec))

    def square_root(self):
        r"""
        Returns the square root of this p-adic number

        INPUT:
            self -- a p-adic element
        OUTPUT:
            p-adic element -- the square root of this p-adic number

            The square root chosen is the one whose reduction mod p is in
            the range [0, p/2).

            Note that because this is a fixed modulus ring, garbage digits
            may be introduced, if either
            (a) the valuation of the input is positive, or
            (b) p = 2.

            If no square root exists, a ValueError is raised.
            (This may be changed later to return an element of an extension
            field.)

        EXAMPLES:
            sage: R = Zp(3,20,'fixed-mod')
            sage: R(0).square_root()
                O(3^20)
            sage: R(1).square_root()
                1 + O(3^20)
            sage: R(2).square_root()
            Traceback (most recent call last):
            ...
            ValueError: element is not a square
            sage: R(4).square_root() == R(-2)
                True
            sage: R(9).square_root()
                3 + O(3^20)
            sage: R2 = Zp(2,20,'fixed-mod')
            sage: R2(0).square_root()
                O(2^20)
            sage: R2(1).square_root()
                1 + O(2^20)
            sage: R2(4).square_root()
                2 + O(2^20)
            sage: R2(9).square_root() == R2(3) or R2(9).square_root() == R2(-3)
                True
            sage: R2(17).square_root()
                1 + 2^3 + 2^5 + 2^6 + 2^7 + 2^9 + 2^10 + 2^13 + 2^16 + 2^17 + O(2^20)
            sage: R3 = Zp(5,20,'fixed-mod', 'terse')
            sage: R3(0).square_root()
                0 + O(5^20)
            sage: R3(1).square_root()
                1 + O(5^20)
            sage: R3(-1).square_root() == R3.teichmuller(2) or R3(-1).square_root() == R3.teichmuller(3)
                True
        """
        #todo: make more efficient
        try:
            # use pari
            return self.parent()(pari(self).sqrt())
        except PariError:
            # todo: should eventually change to return an element of
            # an extension field
            raise ValueError, "element is not a square"

    def unit_part(self):
        r"""
        Returns the unit part of self.

        If the valuation of self is positive, then the high digits of the
        result will be zero.

        INPUT:
            self -- a p-adic element

        OUTPUT:
            p-adic element -- the unit part of self

        EXAMPLES:
            sage: R = Zp(17, 4, 'fixed-mod')
            sage: R(5).unit_part()
            5 + O(17^4)
            sage: R(18*17).unit_part()
            1 + 17 + O(17^4)
            sage: R(0).unit_part()
            O(17^4)
            sage: type(R(5).unit_part())
            <type 'sage.rings.padics.padic_ring_fixed_mod_element.pAdicRingFixedModElement'>
        """
        return self.unit_part_c()

    cdef pAdicRingFixedModElement unit_part_c(pAdicRingFixedModElement self):
        cdef pAdicRingFixedModElement ans
        if mpz_divisible_p(self.value, self.p):
            ans = self._new_c()
            mpz_remove(ans.value, self.value, self.p)
            return ans
        else:
            return self

    def valuation(self):
        """
        Returns the valuation of self.

        If self is zero, the valuation returned is the precision of the ring.

        INPUT:
            self -- a p-adic element
        OUTPUT:
            integer -- the valuation of self.

        EXAMPLES:
            sage: R = Zp(17, 4,'fixed-mod')
            sage: a = R(2*17^2)
            sage: a.valuation()
            2
            sage: R = Zp(5, 4,'fixed-mod')
            sage: R(0).valuation()
            4
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
        """
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.valuation_c())
        return ans

    cdef int valuation_c(self):
        if mpz_sgn(self.value) == 0:
            return mpz_get_si((<Integer>self.parent().precision_cap()).value)
        cdef mpz_t tmp
        cdef int ans
        mpz_init(tmp)
        ans = mpz_remove(tmp, self.value, self.p)
        mpz_clear(tmp)
        return ans

