"""
`p`-Adic Capped Relative Elements

Elements of `p`-Adic Rings with Capped Relative Precision

AUTHORS:

- David Roe: initial version, rewriting to use templates (2012-3-1)
- Genya Zaytman: documentation
- David Harvey: doctests
"""

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/libs/linkages/padics/mpz.pxi"
include "CR_template.pxi"

from sage.libs.pari.gen cimport PariInstance
cdef PariInstance P = sage.libs.pari.all.pari
from sage.rings.finite_rings.integer_mod import Mod

cdef class pAdicCappedRelativeElement(CRElement):
    """
    Constructs new element with given parent and value.

    INPUT:

    - ``x`` -- value to coerce into a capped relative ring or field

    - ``absprec`` -- maximum number of digits of absolute precision

    - ``relprec`` -- maximum number of digits of relative precision

    EXAMPLES::

        sage: R = Zp(5, 10, 'capped-rel')

    Construct from integers::

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

    Construct from rationals::

        sage: R(1/2)
        3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 2*5^8 + 2*5^9 + O(5^10)
        sage: R(-7875/874)
        3*5^3 + 2*5^4 + 2*5^5 + 5^6 + 3*5^7 + 2*5^8 + 3*5^10 + 3*5^11 + 3*5^12 + O(5^13)
        sage: R(15/425)
        Traceback (most recent call last):
        ...
        ValueError: p divides the denominator

    Construct from IntegerMod::

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

    ::

        sage: R(Integers(48)(3))
        Traceback (most recent call last):
        ...
        TypeError: cannot coerce from the given integer mod ring (not a power of the same prime)

    Some other conversions::

        sage: R(R(5))
        5 + O(5^11)

    Construct from Pari objects::

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
        sage: R(pari(R(0,5)))
        O(5^5)

    # todo: doctests for converting from other types of p-adic rings

    """
    def lift(self):
        """
        Return an integer or rational congruent to ``self`` modulo ``self``'s
        precision.  If a rational is returned, its denominator will equal
        ``p^ordp(self)``.

        EXAMPLES::

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

    def __copy__(self):
        """
        Returns a copy of self.

        EXAMPLES::

            sage: a = Zp(5,6)(17); b = copy(a)
            sage: a == b
            True
            sage: a is b
            False
        """
        cdef pAdicCappedRelativeElement ans
        ans = self._new_c()
        ans.relprec = self.relprec
        ans.ordp = self.ordp
        ans._normalized = self._normalized
        mpz_set(ans.unit, self.unit)
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
        return mpz_sgn(self.unit) == -1

    cpdef bint _is_inexact_zero(self) except -1:
        """
        Returns True if this element is indistinguishable from zero but has finite precision.

        EXAMPLES::

            sage: R = Zp(5)
            sage: R(0)._is_inexact_zero()
            False
            sage: R(0,5)._is_inexact_zero()
            True
            sage: R(17)._is_inexact_zero()
            False
        """
        self._normalize()
        return self.relprec == 0 and not self._is_exact_zero()

    def is_zero(self, absprec = None):
        r"""
        Returns whether self is zero modulo $p^{\mbox{absprec}}$.

        If absprec is None, returns True if this element is indistinguishable from zero.

        INPUT:

        - self -- a p-adic element
        - absprec -- (default: None) an integer or None

        OUTPUT:

        - boolean -- whether self is zero

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

        TESTS:

        Check that :trac:`12549` is fixed::

            sage: a = Zp(5)(1) + Zp(5)(-1)
            sage: a.is_zero()
            True
        """
        self._normalize()
        if absprec is None:
            return mpz_sgn(self.unit) <= 0
        if mpz_sgn(self.unit) == -1:
            return True
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        elif mpz_sgn(self.unit) == 0:
            if mpz_fits_slong_p((<Integer>absprec).value) == 0 or self.ordp < mpz_get_si((<Integer>absprec).value):
                raise PrecisionError, "Not enough precision to determine if element is zero"
            else:
                return True
        return self.ordp >= mpz_get_si((<Integer>absprec).value)

    def __nonzero__(self):
        """
        Returns True if self is distinguishable from zero.

        For most applications, explicitly specifying the power of p modulo which the element is supposed to be nonzero is preferable.

        EXAMPLES::

            sage: R = Zp(5); a = R(0); b = R(0,5); c = R(75)
            sage: bool(a), bool(b), bool(c)
            (False, False, True)
        """
        self._normalize()
        return mpz_sgn(self.unit) > 0

    def is_equal_to(self, right, absprec=None):
        r"""
        Returns whether self is equal to right modulo $p^{\mbox{absprec}}$.

        if absprec is None, returns True if self and right are equal to the minimum of their precisions.

        INPUT:

        - self -- a p-adic element
        - right -- a p-addic element
        - absprec -- an integer or None

        OUTPUT:

        - boolean -- whether self is equal to right (modulo $p^{\mbox{absprec}}$)

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
            elif (<pAdicCappedRelativeElement>right).ordp >= aprec:
                return True
            else:
                return False
        if aprec > (<pAdicCappedRelativeElement>self).ordp + (<pAdicCappedRelativeElement>self).relprec:
            raise PrecisionError, "Elements not known to enough precision"
        if mpz_sgn((<pAdicCappedRelativeElement>right).unit) == -1:
            if self.ordp >= aprec:
                return True
            else:
                return False
        if aprec > (<pAdicCappedRelativeElement>right).ordp + (<pAdicCappedRelativeElement>right).relprec:
            raise PrecisionError, "Elements not known to enough precision"
        # We now know that both self and right have enough precision to determine modulo p^absprec
        if self.ordp >= aprec and (<pAdicCappedRelativeElement>right).ordp >= aprec:
            return True
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
        Return an integer or rational congruent to self modulo self's
        precision.  If a rational is returned, its denominator will
        eqaul p^ordp(self).

        INPUT:

        - self -- a p-adic element

        OUTPUT:

        - integer -- a integer congruent to self mod $p^{\mbox{prec}}$

        EXAMPLES::

            sage: R = Zp(7,4,'capped-rel'); a = R(8); a.lift()
            8
            sage: R = Qp(7,4); a = R(8); a.lift()
            8
            sage: R = Qp(7,4); a = R(8/7); a.lift()
            8/7
        """
        return self.lift_c()

    cdef lift_c(self):
        """
        Implementation of lift.

        TESTS::

            sage: O(5^5).lift() #indirect doctest
            0
            sage: R = Qp(5); R(0).lift()
            0
            sage: R(5/9).lift()
            264909532335070
            sage: R(9/5).lift()
            9/5
        """
        cdef Integer ans
        cdef Rational ansr
        if self.ordp >= 0:
            ans = PY_NEW(Integer)
            if self.relprec == 0:
                mpz_set_ui(ans.value, 0)
            else:
                mpz_set(ans.value, self.unit)
                mpz_mul(ans.value, ans.value, self.prime_pow.pow_mpz_t_tmp(self.ordp)[0])
            return ans
        else:
            ansr = PY_NEW(Rational)
            if self.relprec == 0:
                mpq_set_si(ansr.value, 0, 1)
                return self
            else:
                mpz_set(mpq_numref(ansr.value), self.unit)
                mpz_set(mpq_denref(ansr.value), self.prime_pow.pow_mpz_t_tmp(-self.ordp)[0])
            return ansr

    def _pari_(self):
        """
        Converts this element to an equivalent pari element.

        EXAMPLES::

            sage: R = Zp(17, 10); a = ~R(14); pari(a) #indirect doctest
            11 + 3*17 + 17^2 + 6*17^3 + 13*17^4 + 15*17^5 + 10*17^6 + 3*17^7 + 17^8 + 6*17^9 + O(17^10)
            sage: pari(R(0))
            0
            sage: pari(R(0,5))
            O(17^5)
        """
        return self._to_gen()

    cdef pari_gen _to_gen(self):
        """
        Converts this element to an equivalent pari element.

        EXAMPLES::

           sage: R = Zp(5, 10); a = R(17); pari(a) #indirect doctest
           2 + 3*5 + O(5^10)
           sage: pari(R(0))
           0
           sage: pari(R(0,5))
           O(5^5)
        """
        if exactzero(self.ordp):
            return P.new_gen_from_int(0)
        else:
            return P.new_gen_from_padic(self.ordp, self.relprec,
                                        self.prime_pow.prime.value,
                                        self.prime_pow.pow_mpz_t_tmp(self.relprec)[0],
                                        self.unit)
    def _integer_(self, Z=None):
        """
        Returns an integer congruent to this element modulo
        ``p^self.absolute_precision()``.

         EXAMPLES::

            sage: R = Zp(5); a = R(-1); a._integer_()
            95367431640624
         """
        if self.ordp < 0:
            raise ValueError, "Cannot form an integer out of a p-adic field element with negative valuation"
        return self.lift_c()

    def residue(self, absprec=1):
        """
        Reduces this element modulo `p^{\mbox{absprec}}`.

        INPUT:

        - ``absprec`` - an integer (default: ``1``)

        OUTPUT:

        Element of `\ZZ/(p^{\mbox{absprec}} \ZZ)` -- the reduction modulo
        `p^{\mbox{absprec}}`

        EXAMPLES::

            sage: R = Zp(7,4,'capped-rel'); a = R(8); a.residue(1)
            1
            sage: R = Qp(7,4,'capped-rel'); a = R(8); a.residue(1)
            1
            sage: a.residue(6)
            Traceback (most recent call last):
            ...
            PrecisionError: Not enough precision known in order to compute residue.
            sage: b = a/7
            sage: b.residue(1)
            Traceback (most recent call last):
            ...
            ValueError: Element must have non-negative valuation in order to compute residue.
        """
        cdef Integer selfvalue, modulus
        cdef long aprec
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if absprec > self.precision_absolute():
            raise PrecisionError, "Not enough precision known in order to compute residue."
        elif absprec < 0:
            raise ValueError, "cannot reduce modulo a negative power of p"
        aprec = mpz_get_ui((<Integer>absprec).value)
        if self.ordp < 0:
            raise ValueError, "Element must have non-negative valuation in order to compute residue."
        modulus = PY_NEW(Integer)
        mpz_set(modulus.value, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
        selfvalue = PY_NEW(Integer)
        if self.relprec == 0:
            mpz_set_ui(selfvalue.value, 0)
        else:
            # Need to do this better.
            mpz_mul(selfvalue.value, self.prime_pow.pow_mpz_t_tmp(self.ordp)[0], self.unit)
        return Mod(selfvalue, modulus)

def unpickle_pcre_v1(R, unit, ordp, relprec):
    """
    Unpickles a capped relative element.

    EXAMPLES::

        sage: from sage.rings.padics.padic_capped_relative_element import unpickle_pcre_v1
        sage: R = Zp(5)
        sage: a = unpickle_pcre_v1(R, 17, 2, 5); a
        2*5^2 + 3*5^3 + O(5^7)
    """
    return unpickle_cre_v2(pAdicCappedRelativeElement, R, unit, ordp, relprec)

def base_p_list(Integer n, bint pos, PowComputer_class prime_pow):
    """
    Returns a base-`p` list of digits of ``n``.

    INPUT:

    - ``n`` -- a positive Integer.

    - ``pos`` -- a boolean.  If True, then returns the standard base `p` expansion.
                 Otherwise, the digits lie in the range `-p/2` to `p/2`.

    - ``prime_pow`` -- A PowComputer giving the prime.

    EXAMPLES::

        sage: from sage.rings.padics.padic_capped_relative_element import base_p_list
        sage: base_p_list(192837, True, Zp(5).prime_pow)
        [2, 2, 3, 2, 3, 1, 2, 2]
        sage: 2 + 2*5 + 3*5^2 + 2*5^3 + 3*5^4 + 5^5 + 2*5^6 + 2*5^7
        192837
        sage: base_p_list(192837, False, Zp(5).prime_pow)
        [2, 2, -2, -2, -1, 2, 2, 2]
        sage: 2 + 2*5 - 2*5^2 - 2*5^3 - 5^4 + 2*5^5 + 2*5^6 + 2*5^7
        192837
    """
    if mpz_sgn(n.value) < 0:
        raise ValueError("n must be nonnegative")
    return clist(n.value, prime_pow.prec_cap, pos, prime_pow)
