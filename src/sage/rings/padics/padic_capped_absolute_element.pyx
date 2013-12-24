"""
`p`-Adic Capped Absolute Elements

Elements of `p`-Adic Rings with Absolute Precision Cap

AUTHORS:

- David Roe
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
include "CA_template.pxi"

from sage.libs.pari.pari_instance cimport PariInstance
cdef PariInstance P = sage.libs.pari.pari_instance.pari
from sage.rings.finite_rings.integer_mod import Mod

cdef class pAdicCappedAbsoluteElement(CAElement):
    """
    Constructs new element with given parent and value.

    INPUT:

    - ``x`` -- value to coerce into a capped absolute ring

    - ``absprec`` -- maximum number of digits of absolute precision

    - ``relprec`` -- maximum number of digits of relative precision

    EXAMPLES::

        sage: R = ZpCA(3, 5)
        sage: R(2)
        2 + O(3^5)
        sage: R(2, absprec=2)
        2 + O(3^2)
        sage: R(3, relprec=2)
        3 + O(3^3)
        sage: R(Qp(3)(10))
        1 + 3^2 + O(3^5)
        sage: R(pari(6))
        2*3 + O(3^5)
        sage: R(pari(1/2))
        2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5)
        sage: R(1/2)
        2 + 3 + 3^2 + 3^3 + 3^4 + O(3^5)
        sage: R(mod(-1, 3^7))
        2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + O(3^5)
        sage: R(mod(-1, 3^2))
        2 + 2*3 + O(3^2)
        sage: R(3 + O(3^2))
        3 + O(3^2)
    """
    def lift(self):
        """
            sage: R = ZpCA(3)
            sage: R(10).lift()
            10
            sage: R(-1).lift()
            3486784400
        """
        return self.lift_c()

    cdef lift_c(self):
        """
        Implementation of lift.

        TESTS::

            sage: ZpCA(3,3)(1/4).lift() # indirect doctest
            7
        """
        cdef Integer ans = PY_NEW(Integer)
        mpz_set(ans.value, self.value)
        return ans

    def _pari_(self):
        """
        Conversion to pari.

        EXAMPLES::

            sage: R = ZpCA(5)
            sage: pari(R(1777)) #indirect doctest
            2 + 5^2 + 4*5^3 + 2*5^4 + O(5^20)
            sage: pari(R(0,0))
            O(5^0)
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
        # holder is defined in the linkage file
        cdef long val = mpz_remove(holder.value, self.value, self.prime_pow.prime.value)
        return P.new_gen_from_padic(val, self.absprec - val,
                                    self.prime_pow.prime.value,
                                    self.prime_pow.pow_mpz_t_tmp(self.absprec - val)[0],
                                    holder.value)

    def _integer_(self, Z=None):
        r"""
        Converts this element to an integer.

        TESTS::

            sage: R = ZpCA(5, prec = 4); a = R(642); a
            2 + 3*5 + O(5^4)
            sage: a._integer_()
            17
        """
        return self.lift_c()

    def residue(self, absprec=1):
        r"""
        Reduces this element modulo ``p^absprec``.

        INPUT:

        - ``absprec`` - an integer

        OUTPUT:

        element of `\mathbb{Z}/p^{\mbox{absprec}} \mathbb{Z}` -- ``self``
        reduced modulo ``p^absprec``.

        EXAMPLES::

            sage: R = Zp(7,4,'capped-abs'); a = R(8); a.residue(1)
            1
        """
        cdef Integer selfvalue, modulus
        if not PY_TYPE_CHECK(absprec, Integer):
            absprec = Integer(absprec)
        if mpz_cmp_si((<Integer>absprec).value, self.absprec) > 0:
            raise PrecisionError, "Not enough precision known in order to compute residue."
        elif mpz_sgn((<Integer>absprec).value) < 0:
            raise ValueError, "cannot reduce modulo a negative power of p"
        cdef long aprec = mpz_get_ui((<Integer>absprec).value)
        modulus = PY_NEW(Integer)
        mpz_set(modulus.value, self.prime_pow.pow_mpz_t_tmp(aprec)[0])
        selfvalue = PY_NEW(Integer)
        mpz_set(selfvalue.value, self.value)
        return Mod(selfvalue, modulus)

    def multiplicative_order(self):
        r"""
        Returns the minimum possible multiplicative order of this element.

        OUTPUT:
        the multiplicative order of self.  This is the minimum multiplicative
        order of all elements of `\mathbb{Z}_p` lifting ``self`` to infinite
        precision.

        EXAMPLES::

            sage: R = ZpCA(7, 6)
            sage: R(1/3)
            5 + 4*7 + 4*7^2 + 4*7^3 + 4*7^4 + 4*7^5 + O(7^6)
            sage: R(1/3).multiplicative_order()
            +Infinity
            sage: R(7).multiplicative_order()
            +Infinity
            sage: R(1).multiplicative_order()
            1
            sage: R(-1).multiplicative_order()
            2
            sage: R.teichmuller(3).multiplicative_order()
            6
        """
        cdef mpz_t ppow_minus_one
        cdef Integer ans
        if mpz_divisible_p(self.value, self.prime_pow.prime.value):
            return infinity
        if mpz_cmp_ui(self.value, 1) == 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 1)
            return ans
        mpz_init(ppow_minus_one)
        mpz_sub_ui(ppow_minus_one, self.prime_pow.pow_mpz_t_tmp(self.absprec)[0], 1)
        if mpz_cmp(self.value, ppow_minus_one) == 0:
            ans = PY_NEW(Integer)
            mpz_set_ui(ans.value, 2)
            mpz_clear(ppow_minus_one)
            return ans
        # check if self is an approximation to a teichmuller lift:
        mpz_powm(ppow_minus_one, self.value, self.prime_pow.prime.value, self.prime_pow.pow_mpz_t_tmp(self.absprec)[0])
        if mpz_cmp(ppow_minus_one, self.value) == 0:
            mpz_clear(ppow_minus_one)
            return self.residue(1).multiplicative_order()
        else:
            mpz_clear(ppow_minus_one)
            return infinity

def make_pAdicCappedAbsoluteElement(parent, x, absprec):
    """
    Unpickles a capped absolute element.

    EXAMPLES::

        sage: from sage.rings.padics.padic_capped_absolute_element import make_pAdicCappedAbsoluteElement
        sage: R = ZpCA(5)
        sage: a = make_pAdicCappedAbsoluteElement(R, 17*25, 5); a
        2*5^2 + 3*5^3 + O(5^5)
    """
    return unpickle_cae_v2(pAdicCappedAbsoluteElement, parent, x, absprec)
