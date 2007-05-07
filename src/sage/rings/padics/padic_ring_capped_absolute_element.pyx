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

cdef class pAdicRingCappedAbsoluteElement(pAdicBaseGenericElement):
    def __init__(pAdicRingCappedAbsoluteElement self, parent, x, absprec=infinity, relprec = infinity, empty=False):
        mpz_init(self.value)
        mpz_init_set(self.p, (<Integer>parent.prime()).value)
        mpz_init(self.modulus)
        CommutativeRingElement.__init__(self,parent)
        if empty:
            return
        cdef PowComputer_class powerer
        powerer = <PowComputer_class> parent.prime_pow
        if not absprec is infinity and not relprec is infinity:
            raise ValueError, "can only specify one of absprec and relprec"
        if relprec is infinity:
            if absprec > parent.precision_cap():
                absprec = parent.precision_cap()

        if isinstance(x, pAdicGenericElement):
            if x.valuation() < 0:
                raise ValueError, "element valuation cannot be negative."
            if parent.prime() != x.parent().prime():
                raise TypeError, "Cannot coerce between p-adic parents with different primes"
        if isinstance(x, pAdicLazyElement):
            if relprec is infinity:
                absprec = min(absprec, parent.precision_cap())
                try:
                    x.set_precision_absolute(absprec)
                except PrecisionError:
                    pass
            else:
                try:
                    x.set_precision_relative(relprec)
                except PrecisionError:
                    try:
                        x.set_precision_absolute(parent.precision_cap())
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
        cdef int k
        if sage.rings.integer_mod.is_IntegerMod(x):
            mpz_init_set(modulus, (<Integer>x.modulus()).value)
            k = mpz_remove(modulus, modulus, self.p)
            if mpz_cmp_ui(modulus, 1) == 0:
                # There's a bug here if the user tries to lift to absprec = parent.precision_cap()
                if absprec == parent.precision_cap() and relprec is infinity: #this allows you to lift integer_mod elements to higher precision than they are actually defined.
                    absprec = min(absprec, k)
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
            absprec = min(absprec, val + relprec, parent.precision_cap())
            self.set_precs(mpz_get_ui((<Integer>absprec).value))
            tmp = PY_NEW(Integer)
            mpz_set(tmp.value, self.modulus)
            tmp = <Integer> x % tmp
            self.set_value_from_mpz(tmp.value)
            return
        if isinstance(x, Integer):
            self.set_from_Integers(x, min(absprec, val + relprec, parent.precision_cap()))
        else:
            raise TypeError, "unable to create p-adic element from %s of type %s"%(x, type(x))

    def __dealloc__(self):
        # If we switch to using the pointer copying method for setting self.p, the corresponding mpz_clear here should be removed
        mpz_clear(self.p)
        mpz_clear(self.value)

    cdef void set_precs(pAdicRingCappedAbsoluteElement self, unsigned int absprec):
        """
        Sets self.absprec and self.modulus
        """
        cdef PowComputer_class powerer
        powerer = <PowComputer_class>self.parent().prime_pow
        self.absprec = absprec
        powerer.pow_mpz_ui(self.modulus, absprec)

    cdef void set_value_from_mpz(pAdicRingCappedAbsoluteElement self, mpz_t value):
        """
        Assuming that self.modulus is set, sets self.value
        """
        if mpz_sgn(value) == -1 or mpz_cmp(value, self.modulus) >= 0:
            mpz_mod(self.value, value, self.modulus)
        else:
            mpz_set(self.value, value)

    cdef void set_from_Integers(pAdicRingCappedAbsoluteElement self, Integer value, Integer absprec):
        """
        Set self.value, self.absprec and self.modulus.
        """
        self.set_precs(mpz_get_ui(absprec.value))
        self.set_value_from_mpz(value.value)

    cdef pAdicRingCappedAbsoluteElement _new_c(self):
        cdef pAdicRingCappedAbsoluteElement x
        x = PY_NEW(pAdicRingCappedAbsoluteElement)
        x._parent = self._parent
        mpz_init(x.value)
        # I wanted to just copy pointers:
        # x.p = self.p
        # ... but Pyrex apparently won't let me.
        mpz_init_set(x.p, self.p)
        mpz_init(x.modulus)
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
        cdef pAdicRingCappedAbsoluteElement ans
        ans = self._new_c()
        mpz_set(ans.modulus, self.modulus)
        ans.absprec = self.absprec
        mpz_sub(ans.value, self.modulus, self.value)
        return ans

    def __pow__(pAdicRingCappedAbsoluteElement self, right, dummy):
        cdef Integer new, val, absprec
        cdef mpz_t tmp
        new = Integer(right) #Need to make sure that this works for p-adic exponents
        val = self.valuation()
        if (val > 0) and isinstance(right, pAdicBaseGenericElement):
            raise ValueError, "Can only have p-adic exponent if base is a unit"
        if (new < 0):
            return (~self).__pow__(-new)
        cdef pAdicRingCappedAbsoluteElement ans
        cdef Integer preccap
        preccap = <Integer>self.parent().precision_cap()
        ans = self._new_c()
        if val > 0:
            if new * val >= self.parent().precision_cap():
                ans.set_from_Integers(Integer(0), self.parent().precision_cap())
            else:
                absprec = <Integer> min(self.precision_relative() + new * val, self.parent().precision_cap())
                ans.set_precs(mpz_get_ui(absprec.value))
                _sig_on
                mpz_powm(ans.value, self.value, new.value, ans.modulus)
                _sig_off
        else:
            ans.absprec = self.absprec
            mpz_set(ans.modulus, self.modulus)
            _sig_on
            mpz_powm(ans.value, self.value, new.value, ans.modulus)
            _sig_off
        return ans

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        cdef pAdicRingCappedAbsoluteElement ans
        ans = self._new_c()
        if self.absprec < (<pAdicRingCappedAbsoluteElement>right).absprec:
            ans.absprec = self.absprec
            mpz_set(ans.modulus, self.modulus)
        else:
            ans.absprec = (<pAdicRingCappedAbsoluteElement>right).absprec
            mpz_set(ans.modulus, (<pAdicRingCappedAbsoluteElement>right).modulus)
        mpz_add(ans.value, self.value, (<pAdicRingCappedAbsoluteElement>right).value)
        if mpz_cmp(ans.value, ans.modulus) >= 0:
            mpz_mod(ans.value, ans.value, ans.modulus)
        return ans

    cdef RingElement _div_c_impl(self, RingElement right):
        return self * (<pAdicRingCappedAbsoluteElement>right)._invert_c_impl()

    def __lshift__(pAdicRingCappedAbsoluteElement self, shift):
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
        cdef pAdicRingCappedAbsoluteElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_sint_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            if mpz_sgn((<Integer>shift).value) == 1:
                ans.set_precs(mpz_get_ui((<Integer>self.parent().precision_cap()).value))
            else:
                ans.set_precs(0)
            return ans
        return self._lshift_c(mpz_get_si((<Integer>shift).value))

    cdef pAdicRingCappedAbsoluteElement _lshift_c(pAdicRingCappedAbsoluteElement self, int shift):
        cdef int prec_cap, ansprec
        cdef pAdicRingCappedAbsoluteElement ans
        cdef PowComputer_class powerer
        if shift < 0:
            return self._rshift_c(-shift)
        prec_cap = mpz_get_si((<Integer>self.parent().precision_cap()).value)
        if shift >= prec_cap: # if we ever start caching valuation, this check should change
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            ans.set_precs(prec_cap)
            return ans
        elif shift > 0:
            powerer = self.parent().prime_pow
            ans = self._new_c()
            ansprec = shift + self.absprec
            if ansprec > prec_cap:
                ansprec = prec_cap
            ans.set_precs(ansprec)
            powerer.pow_mpz_ui(ans.value, shift)
            mpz_mul(ans.value, ans.value, self.value)
            if mpz_cmp(ans.value, ans.modulus) >= 0:
                mpz_mod(ans.value, ans.value, ans.modulus)
            return ans
        else:
            return self

    def __rshift__(pAdicRingCappedAbsoluteElement self, shift):
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
        cdef pAdicRingCappedAbsoluteElement ans
        if not PY_TYPE_CHECK(shift, Integer):
            shift = Integer(shift)
        if mpz_fits_sint_p((<Integer>shift).value) == 0:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            if mpz_sgn((<Integer>shift).value) == -1:
                ans.set_precs(mpz_get_ui((<Integer>self.parent().precision_cap()).value))
            else:
                ans.set_precs(0)
            return ans
        return self._rshift_c(mpz_get_si((<Integer>shift).value))

    cdef pAdicRingCappedAbsoluteElement _rshift_c(pAdicRingCappedAbsoluteElement self, int shift):
        cdef int prec_cap, ansprec
        cdef PowComputer_class powerer
        cdef pAdicRingCappedAbsoluteElement ans
        if shift < 0:
            return self._lshift_c(-shift)
        prec_cap = mpz_get_si((<Integer>self.parent().precision_cap()).value)
        if shift >= self.absprec:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            ans.set_precs(0)
            return ans
        elif shift > 0:
            powerer = self.parent().prime_pow
            ans = self._new_c()
            ansprec = self.absprec - shift
            ans.set_precs(ansprec)
            powerer.pow_mpz_ui(ans.value, shift)
            mpz_fdiv_q(ans.value, self.value, ans.value)
            return ans
        else:
            return self

    cdef RingElement _mul_c_impl(self, RingElement right):
        """
        EXAMPLES:

        """
        cdef pAdicRingCappedAbsoluteElement ans
        cdef unsigned int sval, rval, prec1, prec2, prec_cap
        ans = (<pAdicRingCappedAbsoluteElement>self)._new_c()
        prec_cap = mpz_get_ui((<Integer>self.parent().precision_cap()).value)
        sval = (<pAdicRingCappedAbsoluteElement>self).valuation_c()
        rval = (<pAdicRingCappedAbsoluteElement>right).valuation_c()
        # Need to think about the possibility that the following overflows (too big for int).  I think we're okay since they're unsigned.
        prec1 = sval + (<pAdicRingCappedAbsoluteElement>right).absprec
        prec2 = rval + (<pAdicRingCappedAbsoluteElement>self).absprec
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
        mpz_mul(ans.value, self.value, (<pAdicRingCappedAbsoluteElement>right).value)
        mpz_mod(ans.value, ans.value, ans.modulus)
        return ans

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        cdef pAdicRingCappedAbsoluteElement ans
        ans = self._new_c()
        if self.absprec < (<pAdicRingCappedAbsoluteElement>right).absprec:
            ans.absprec = self.absprec
            mpz_set(ans.modulus, self.modulus)
        else:
            ans.absprec = (<pAdicRingCappedAbsoluteElement>right).absprec
            mpz_set(ans.modulus, (<pAdicRingCappedAbsoluteElement>right).modulus)
        mpz_sub(ans.value, self.value, (<pAdicRingCappedAbsoluteElement>right).value)
        if mpz_sgn(ans.value) == -1:
            mpz_mod(ans.value, ans.value, ans.modulus)
        return ans

    def add_bigoh(pAdicRingCappedAbsoluteElement self, prec):
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
        cdef pAdicRingCappedAbsoluteElement ans
        cdef int newprec
        if not PY_TYPE_CHECK(prec, Integer):
            prec = Integer(prec)
        if mpz_fits_sint_p((<Integer>prec).value) == 0:
            if mpz_sgn((<Integer>prec).value) == -1:
                ans = self._new_c()
                ans.set_precs(0)
                mpz_set_ui(ans.value, 0)
                return ans
            else:
                return self
        newprec = mpz_get_si((<Integer>prec).value)
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
            if mpz_cmp(ans.value, ans.modulus) >= 0:
                mpz_mod(ans.value, ans.value, ans.modulus)
            return ans

    def copy(pAdicRingCappedAbsoluteElement self):
        cdef pAdicRingCappedAbsoluteElement ans
        ans = self._new_c()
        ans.absprec = self.absprec
        mpz_set(ans.modulus, self.modulus)
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
            return bool(mpz_sgn(self.value) == 0)
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_cmp_ui((<Integer>absprec).value, self.absprec) > 0:
            raise PrecisionError, "Not enough precision to determine if element is zero"
        cdef mpz_t ppow, tmp
        cdef PowComputer_class powerer
        mpz_init(ppow)
        mpz_init(tmp)
        powerer = <PowComputer_class> self.parent().prime_pow
        powerer.pow_mpz_mpz(ppow, (<Integer>absprec).value)
        mpz_mod(tmp, self.value, ppow)
        mpz_clear(ppow)
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
            right -- a p-addic element
            absprec -- an integer
        OUTPUT:
            boolean -- whether self is equal to right

        """
        cdef int intprec
        if absprec is None:
            if self.absprec < right.absprec:
                intprec = self.absprec
            else:
                intprec = right.absprec
        else:
            if not PY_TYPE_CHECK(absprec, Integer):
                absprec = Integer(absprec)
            if mpz_fits_sint_p((<Integer>absprec).value) == 0:
                if mpz_sgn((<Integer>absprec).value) == 1:
                    raise PrecisionError, "Elements are not known to specified precision"
                else:
                    return True
            intprec = mpz_get_si((<Integer>absprec).value)
            if intprec > self.absprec or intprec > right.absprec:
                raise PrecisionError, "Elements are not known to specified precision"
        cdef mpz_t ppow, tmp1, tmp2
        cdef PowComputer_class powerer
        mpz_init(ppow)
        mpz_init(tmp1)
        mpz_init(tmp2)
        powerer = <PowComputer_class> self.parent().prime_pow
        powerer.pow_mpz_ui(ppow, intprec)
        mpz_mod(tmp1, self.value, ppow)
        mpz_mod(tmp2, (<pAdicRingCappedAbsoluteElement>right).value, ppow)
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
        return self.lift_c()

    cdef Integer lift_c(pAdicRingCappedAbsoluteElement self):
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set(ans.value, self.value)
        return ans

    def lift_to_precision(self, absprec):
        cdef pAdicRingCappedAbsoluteElement ans
        cdef int prec_cap
        cdef int dest_prec
        prec_cap = mpz_get_ui((<Integer>self.parent().precision_cap()).value)
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_fits_sint_p((<Integer>absprec).value) == 0:
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

    def list(pAdicRingCappedAbsoluteElement self, lift_mode = 'simple'):
        """
        Returns a list of coeficiants of p starting with $p^0$
        INPUT:
            self -- a p-adic element
        OUTPUT:
            list -- the list of coeficients of self
        EXAMPLES:
            sage: R = Zp(7,4,'capped-abs'); a = R(2*7+7**2); a.list()
                [0, 2, 1]

        NOTE:
            this differs from the list method of padic_field_element
            use slice operators to get a particular range

        """
        if lift_mode == 'teichmuller':
            return self.teichmuller_list()
        else:
            return self.base_p_list(self.value, self.p, lift_mode, self.parent().prime_pow, mpz_get_si((<Integer>self.parent().precision_cap()).value))

    cdef object teichmuller_list(pAdicRingCappedAbsoluteElement self):
        # May eventually want to add a dict to store teichmuller lifts already seen, if p small enough
        cdef int curpower, preccap
        cdef mpz_t ppow, tmp
        cdef pAdicRingCappedAbsoluteElement list_elt
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
            list_elt.set_precs(preccap)
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

    def _teichmuller_set(self, Integer n, Integer absprec):
        cdef mpz_t ppow
        mpz_set(self.value, n.value)
        if mpz_fits_sint_p(absprec.value) == 0:
            raise ValueError, "cannot compute teichmuller lift to that high precision"
        if mpz_sgn(absprec.value) != 1:
            raise ValueError, "can only compute to positive precision"
        self.set_precs(mpz_get_si(absprec.value))
        sage.rings.padics.padic_generic_element.teichmuller_set_c(self.value, self.p, self.modulus)

    def log_artin_hasse(self):
        raise NotImplementedError

    def multiplicative_order(self):
        r"""
        Returns the multiplicative order of self, where self is considered to be one if it is one modulo $p^{\mbox{prec}}$.

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
        if mpz_cmp_ui(self.value, 1) == 0:
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
            sage: R = Zp(7,4,'capped-abs'); a = R(2*7+7**2); a.padded_list(5)
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

    def residue(self, prec):
        r"""
        Reduces this mod $p^prec$

        INPUT:
            self -- a p-adic element
            prec - an integer

        OUTPUT:
            element of Z/(p^prec Z) -- self reduced mod p^prec

        EXAMPLES:
            sage: R = Zp(7,4,'capped-abs'); a = R(8); a.residue(1)
            1
        """
        if prec > self.precision_absolute():
            raise PrecisionError, "Not enough precision known in order to compute residue."
        elif prec < 0:
            raise ValueError, "cannot reduce modulo a negative power of p"
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
        EXAMPLES:
            sage: R = Zp(3,20,'capped-abs')
            sage: R(0).square_root()
                O(3^10)
            sage: R(1).square_root()
                1 + O(3^20)
            sage: R(4).square_root() == R(-2)
                True
            sage: R(9).square_root()
                3 + O(3^19)
            sage: R2 = Zp(2,20,'capped-abs')
            sage: R2(0).square_root()
                O(2^10)
            sage: R2(1).square_root()
                1 + O(2^19)
            sage: R2(4).square_root()
                2 + O(2^18)
            sage: R2(9).square_root() == R2(3) or R2(9).square_root() == R2(-3)
                True
            sage: R2(17).square_root()
                1 + 2^3 + 2^5 + 2^6 + 2^7 + 2^9 + 2^10 + 2^13 + 2^16 + 2^17 + O(2^19)
            sage: R3 = Zp(5,20,'capped-abs')
            sage: R3(0).square_root()
                O(5^10)
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
            raise ValueError, "element is not a square" # should eventually change to return an element of an extension field

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
                <type 'sage.rings.padics.padic_ring_capped_absolute_element.pAdicRingCappedAbsoluteElement'>
        """
        return self.unit_part_c()

    cdef pAdicRingCappedAbsoluteElement unit_part_c(pAdicRingCappedAbsoluteElement self):
        cdef pAdicRingCappedAbsoluteElement ans
        cdef int v
        if mpz_sgn(self.value) == 0:
            ans = self._new_c()
            mpz_set_ui(ans.value, 0)
            ans.set_precs(0)
            return ans
        elif mpz_divisible_p(self.value, self.p):
            ans = self._new_c()
            v = mpz_remove(ans.value, self.value, self.p)
            ans.set_precs(self.absprec - v)
            return ans
        else:
            return self

    def valuation(self):
        cdef Integer ans
        ans = PY_NEW(Integer)
        mpz_set_si(ans.value, self.valuation_c())
        return ans

    cdef int valuation_c(self):
        if mpz_sgn(self.value) == 0:
            return self.absprec
        cdef mpz_t tmp
        cdef int ans
        mpz_init(tmp)
        ans = mpz_remove(tmp, self.value, self.p)
        mpz_clear(tmp)
        return ans

    def val_unit(self):
        return self.val_unit_c()

    cdef val_unit_c(self):
        cdef pAdicRingCappedAbsoluteElement unit
        cdef Integer val
        cdef int v
        val = PY_NEW(Integer)
        if mpz_sgn(self.value) == 0:
            unit = self._new_c()
            mpz_set_ui(unit.value, 0)
            unit.set_precs(0)
            mpz_set_ui(val.value, self.absprec)
            return (val, unit)
        elif mpz_divisible_p(self.value, self.p):
            unit = self._new_c()
            v = mpz_remove(unit.value, self.value, self.p)
            unit.set_precs(self.absprec - v)
            mpz_set_ui(val.value, v)
            return (val, unit)
        else:
            mpz_set_ui(val.value, 0)
            return (val, self)

    def __hash__(self):
        return self._hash()

    cdef long _hash(self) except -1:
        cdef Integer ans
        if self.absprec == 0:
            # This is so that different primes are distinguished.  If someone else can think of a better idea, go for it.
            return hash(self.parent().prime_pow(self.parent().precision_cap()))
        else:
            ans = PY_NEW(Integer)
            mpz_xor(ans.value, self.modulus, self.value)
            return hash(ans)
